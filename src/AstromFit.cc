#include <iostream>
#include <algorithm>
#include "lsst/meas/simastrom/AstromFit.h"
#include "lsst/meas/simastrom/Associations.h"
#include "lsst/meas/simastrom/Mapping.h"
//#include "preferences.h"
#include "lsst/meas/simastrom/Gtransfo.h"
#include "Eigen/Sparse"
//#include "Eigen/CholmodSupport" // to switch to cholmod
#include <time.h> // for clock
#include "lsst/pex/exceptions.h"
#include <fstream>
#include "lsst/meas/simastrom/Matvect.h"
#include "lsst/meas/simastrom/Tripletlist.h"

typedef Eigen::SparseMatrix<double> SpMat;

using namespace std;

static double sqr(const double &x) {return x*x;}

const double posErrorIncrement=0.02;

namespace lsst {
namespace meas {
namespace simastrom {


AstromFit::AstromFit(Associations &A, DistortionModel *D) : 
  _assoc(A),  _distortionModel(D)
{
  _LastNTrip = 0;
  _JDRef = 0;

  _referenceColor = 0; 
  _sigCol = 0;
  unsigned count = 0;
  for (auto i=_assoc.fittedStarList.begin(); 
       i!= _assoc.fittedStarList.end() ; ++i) 
    {
      _referenceColor += (*i)->color;
      _sigCol += sqr((*i)->color);
      count++;
    }
  if (count)
    {
      _referenceColor /= double(count);
      if (_sigCol>0) _sigCol = sqrt(_sigCol/count - sqr(_referenceColor));
    }
  cout << "INFO: reference Color : " << _referenceColor << " sig " << _sigCol << " " << count << endl;

  _nRefrac = _assoc.NBands();
  _refracCoefficient.resize(_nRefrac,0);

  _nMeasuredStars = 0;
  // The various _npar... are initialized in AssignIndices.
  // Although there is no reason to adress them before one might be tempted by
  // evaluating a Chi2 rightaway, .. which uses these counts, so:
  AssignIndices("");

}



#define NPAR_PM 2

/* ! this routine is used in 3 instances: when computing
the derivatives, when computing the Chi2, when filling a tuple.
*/
Point AstromFit::TransformFittedStar(const FittedStar &F,
				     const Gtransfo * Sky2TP,
				     const Point &RefractionVector,
				     const double RefractionCoeff,
				     const double Jd) const
{
  Point fittedStarInTP =  Sky2TP->apply(F);
  if (F.mightMove)
    {
      fittedStarInTP.x += F.pmx*Jd;
      fittedStarInTP.y += F.pmy*Jd;
    }
  // account for atmospheric refraction: does nothing if color 
  // have not been assigned
  // the color definition shouldbe the same when computing derivatives
  double color = F.color - _referenceColor;
  fittedStarInTP.x += RefractionVector.x * color * RefractionCoeff;
  fittedStarInTP.y += RefractionVector.y * color * RefractionCoeff;
  return fittedStarInTP;
}

/*! this is the first implementation of an error "model". 
  We'll certainly have to upgrade it. MeasuredStar provided
in case we need the mag.  */
static void TweakAstromMeasurementErrors(FatPoint &P, const MeasuredStar &Ms)
{
  static bool called=false;
  static double increment = 0;
  if (!called)
    {
      increment = sqr(posErrorIncrement);
      called = true;
    }
  P.vx += increment;
  P.vy += increment;
}

static bool heavyDebug = false;
static unsigned fsIndexDebug = 0;



// we could consider computing the chi2 here.
// (although it is not extremely useful)
void AstromFit::LSDerivatives(const CcdImage &Ccd, 
			      TripletList &TList, Eigen::VectorXd &Rhs) const
{
  /***************************************************************************/
  /**  Changes in this routine should be reflected into AccumulateStatImage  */
  /***************************************************************************/
  /* Setup */ 
  cout << "LSDerivative three " << endl; 
  // 1 : get the Mapping's
  const Mapping *mapping = _distortionModel->GetMapping(Ccd);
  unsigned npar_mapping = (_fittingDistortions) ? mapping->Npar() : 0;
  unsigned npar_pos = (_fittingPos) ? 2 : 0;
  unsigned npar_refrac = (_fittingRefrac) ? 1 : 0;
  unsigned npar_pm = (_fittingPM) ? NPAR_PM : 0;
  unsigned npar_tot =  npar_mapping + npar_pos + npar_refrac + npar_pm;
  // if (npar_tot == 0) this CcdImage does not contribute 
  // any contraint to the fit, so :
  if (npar_tot == 0) return;
  vector<unsigned> indices(npar_tot,-1);
  if (_fittingDistortions)  mapping->GetMappingIndices(indices);

  // proper motion stuff
  double jd = Ccd.JD() - _JDRef;
  // refraction stuff
  Point refractionVector = Ccd.ParallacticVector();
  unsigned iband = Ccd.BandRank();
  double refractionCoefficient = _refracCoefficient.at(iband);
  // transformation from sky to TP
  const Gtransfo* sky2TP = _distortionModel->Sky2TP(mapping, Ccd);
  // reserve matrices once for all measurements
  GtransfoLin dypdy;
  // the shape of h (et al) is required this way in order to able to 
  // separate derivatives along x and y as vectors.
  Eigen::MatrixX2d h(npar_tot,2), halpha(npar_tot,2), hw(npar_tot,2); 
  Eigen::Matrix2d transW(2,2);
  Eigen::Matrix2d alpha(2,2);
  Eigen::VectorXd grad(npar_tot);
  cout << "alors npar_tot " << npar_tot << " " << npar_mapping << std::endl;
  // current position in the Jacobian
  unsigned kTriplets = TList.NextFreeIndex();
  const MeasuredStarList &catalog = Ccd.CatalogForFit();

  cout << "Catalog size " << catalog.size() << endl;
  int icat=0;
  for (auto i = catalog.begin(); i!= catalog.end(); ++i)
    {
      const MeasuredStar& ms = **i;
      if (!ms.IsValid()) continue;
      icat++;
      // tweak the measurement errors
      FatPoint inPos = ms;
      TweakAstromMeasurementErrors(inPos, ms);
      h.setZero(); // we cannot be sure that all entries will be overwritten.
      FatPoint outPos;
      // should *not* fill h if WhatToFit excludes mapping parameters.
      if (_fittingDistortions) 
	  mapping->ComputeTransformAndDerivatives(inPos, outPos, h);
      else mapping->TransformPosAndErrors(inPos,outPos);

      unsigned ipar = npar_mapping;
      double det = outPos.vx*outPos.vy-sqr(outPos.vxy);
      if (det <=0 || outPos.vx <=0 || outPos.vy<=0) {
	cout << " WARNING: inconsistent measurement errors :drop measurement at " << Point(ms) << " in image " << Ccd.Name() << endl;
	continue;
      }	
      transW(0,0) = outPos.vy/det;
      transW(1,1) = outPos.vx/det;
      transW(0,1) = transW(1,0) = -outPos.vxy/det;
      // compute alpha, a triangular square root
      // of transW (i.e. a Cholesky factor)
      alpha(0,0) = sqrt(transW(0,0));
      // checked that  alpha*alphaT = transW
      alpha(1,0) = transW(0,1)/alpha(0,0); 
      alpha(1,1) = 1./sqrt(det*transW(0,0));
      alpha(0,1) = 0;
      
      const FittedStar *fs = ms.GetFittedStar();
      // DEBUG
      if (heavyDebug && fsIndexDebug == fs->IndexInMatrix())
	{
	  double x = 1;
	}
  
      Point fittedStarInTP = TransformFittedStar(*fs, sky2TP,
						 refractionVector, 
						 refractionCoefficient,
						 jd);

      // compute derivative of TP position w.r.t sky position ....
      if (npar_pos>0) // ... if actually fitting FittedStar position
	{
	  sky2TP->Derivative(*fs, dypdy, 1e-3);
	  // sign checked
	  // TODO Still have to check with non trivial non-diagonal terms
	  h(npar_mapping,0) = -dypdy.A11();
	  h(npar_mapping+1,0) = -dypdy.A12();
	  h(npar_mapping,1) = -dypdy.A21();	  
	  h(npar_mapping+1,1) = -dypdy.A22();
	  indices[npar_mapping] = fs->IndexInMatrix();
	  indices.at(npar_mapping+1) = fs->IndexInMatrix()+1;
	  ipar += npar_pos;
	}
      if (_fittingPM) // should add  "&& fs->mightMove" in the test ? 
	{
	  h(ipar,0) = -jd; // Sign unchecked but consistent with above
	  h(ipar+1, 1) = -jd;
	  indices[ipar] = fs->IndexInMatrix()+2;
	  indices[ipar+1] = fs->IndexInMatrix()+3;
	  ipar+= npar_pm;
	}
      if (_fittingRefrac)
	{
	  /* if the definition of color changes, it has to remain
	     consistent with TransformFittedStar */
	  double color = fs->color - _referenceColor;
	  // sign checked 
	  h(ipar,0) = -refractionVector.x*color;
	  h(ipar,1) = -refractionVector.y*color;
	  indices[ipar] = _refracPosInMatrix+iband;
	  ipar += 1;
	}

      // We can now compute the residual
      Eigen::Vector2d res(fittedStarInTP.x-outPos.x, fittedStarInTP.y-outPos.y);

      // do not write grad = h*transW*res to avoid 
      // dynamic allocation of a temporary
      halpha = h*alpha;
      hw = h*transW;
      grad = hw*res;
      // now feed in triplets and Rhs
      cout << "allo? " << npar_tot << endl;
      for (unsigned ipar=0; ipar<npar_tot; ++ipar)
	{
	  for (unsigned  ic=0; ic<2; ++ic)
	    {
	      double val = halpha(ipar,ic);
	      if (val ==0) continue;
#if (TRIPLET_INTERNAL_COORD == COL)
	      cout << "add tripleta " << kTriplets+ic << " " << indices[ipar] << " " << val << " " << ipar << " " << ic << endl;
	      TList.AddTriplet(indices[ipar], kTriplets+ic,val);
#else
	      cout << "add triplet " << kTriplets+ic << " " << indices[ipar] << " " << val << endl;
	      TList.AddTriplet(kTriplets+ic, indices[ipar], val);
#endif
	    }
	  Rhs(indices[ipar]) += grad(ipar); 
	  if (heavyDebug && fsIndexDebug == fs->IndexInMatrix() && (indices[ipar] == fsIndexDebug || indices[ipar] == fsIndexDebug+1))
	    {
	      cout << " DEBUG : res " << res(0) << ' ' << res(1) << endl;
	      cout << " DEBUG : contribs g, j " << grad(ipar) << ' ' <<  halpha(ipar,0) << ' ' <<  halpha(ipar,1) << endl;
	      unsigned ic = (indices[ipar] == fsIndexDebug) ? 0 : 1;
	      cout << " DEBUG : h alpha " << h(ipar,0) << ' ' << h(ipar,1) << ' ' << alpha(0, ic) << ' ' << alpha(1,ic) << endl;
	    }
	  //#warning : flipped sign of Gradient
	  //	  Rhs(indices[ipar]) -= grad(ipar); 
	}
      kTriplets += 2; // each measurement contributes 2 columns in the Jacobian
    }

  cout << "icat " << icat << endl;
  TList.SetNextFreeIndex(kTriplets);
}    

// we could consider computing the chi2 here.
// (although it is not extremely useful)
void AstromFit::LSDerivatives(const CcdImageList  &L, 
			      TripletList &TList, Eigen::VectorXd &Rhs)
{ 
  cout << "LSDerivative two " << endl; 
  
  for (auto im=L.cbegin(); im!=L.end() ; ++im)
    {
      LSDerivatives(**im, TList, Rhs);
    }
}

#define HACK_REF_ERRORS 1. // used to isolate the measurement or ref terms

//! this routine computes the derivatives of all LS terms, including the ones that refer to references stars, if any
void AstromFit::LSDerivatives(TripletList &TList, Eigen::VectorXd &Rhs)
{
  cout << "LSDerivative one " << endl; 
  //the terms involving fittedstars and measurements
  LSDerivatives(_assoc.TheCcdImageList(), TList,Rhs);
  // terms involving fitted stars and reference stars.
  // they only provide derivatives if we are fitting positions:
  if (! _fittingPos) return;
  /* the other case where the accumulation of derivatives stops 
     here is when there are no RefStars */
  if (_assoc.refStarList.size() == 0) return;
  const FittedStarList &fsl = _assoc.fittedStarList;
  Eigen::Matrix2d w(2,2);
  Eigen::Matrix2d alpha(2,2);
  Eigen::Matrix2d h(2,2), halpha(2,2), hw(2,2);
  GtransfoLin der;
  Eigen::Vector2d res,grad;
  unsigned indices[2+NPAR_PM];
  unsigned kTriplets = TList.NextFreeIndex();
  /* We cannot use the spherical coordinates directly to evaluate
     Euclidean distances, we have to use a projector on some plane in
     order to express least squares. No projecting could lead to a
     disaster around the poles or across alpha=0.  So we need a
     projector. We construct a projector and will change its
     projection point at every object */
  TanRaDec2Pix proj(GtransfoLin(), Point(0.,0.));
  for (auto i = fsl.cbegin(); i!= fsl.end(); ++i)
    {
      const FittedStar &fs = **i;
      const RefStar *rs = fs.GetRefStar();
      if (rs == NULL) continue;
      proj.SetTangentPoint(fs);
      // fs projects to (0,0), no need to compute its transform.
      FatPoint rsProj;
      proj.TransformPosAndErrors(*rs, rsProj);
      proj.Derivative(fs, der, 1e-4);
      // sign checked. TODO check that the off-diagonal terms are OK.
      h(0,0) = -der.A11();
      h(1,0) = -der.A12();
      h(1,0) = -der.A21();
      h(1,1) = -der.A22();
      // TO DO : account for proper motions.
      double det = rsProj.vx*rsProj.vy - sqr(rsProj.vxy);
      if (rsProj.vx <=0 || rsProj.vy <=0 || det <= 0) 
	{
	  cout << " WARNING: Ref star error matrix not posdef:  " << endl
	       << *rs << endl;
	  continue;
	}
      w(0,0) = rsProj.vy/det;
      w(0,1) = w(1,0) = -rsProj.vxy/det;
      w(1,1) = rsProj.vx/det;
      w *= HACK_REF_ERRORS;
      det /= sqr(HACK_REF_ERRORS);
      // compute alpha, a triangular square root
      // of w (i.e. a Cholesky factor)
      alpha(0,0) = sqrt(w(0,0));
      // checked that  alpha*alphaT = transW
      alpha(1,0) = w(0,1)/alpha(0,0); 
      alpha(1,1) = 1./sqrt(det*w(0,0));
      alpha(0,1) = 0;
      indices[0] = fs.IndexInMatrix();
      indices[1] = fs.IndexInMatrix()+1;
      unsigned npar_tot = 2;
      /* TODO: account here for proper motions in the reference
      catalog. We can code the effect and set the value to 0. Most
      (all?)  catalogs do not even come with a reference epoch. Gaia
      will change that. When refraction enters into the game, one should 
      pay attention to the orientation of the frame */

      /* The residual should be Proj(fs)-Proj(*rs) in order to be consistent
	 with the measurement terms. Since P(fs) = 0, we have: */
      res[0] = -rsProj.x;
      res[1] = -rsProj.y;
      halpha = h*alpha;
      // grad = h*w*res
      hw = h*w;
      grad = hw*res;
      // now feed in triplets and Rhs
      for (unsigned ipar=0; ipar<npar_tot; ++ipar)
	{
	  for (unsigned  ic=0; ic<2; ++ic)
	    {
	      double val = halpha(ipar,ic);
	      if (val ==0) continue;
#if (TRIPLET_INTERNAL_COORD == COL)
	      TList.AddTriplet(indices[ipar], kTriplets+ic,val);
#else
	      TList.AddTriplet(kTriplets+ic, indices[ipar], val);
#endif
	    }
	  Rhs(indices[ipar]) += grad(ipar); 
	}
      kTriplets += 2; // each measurement contributes 2 columns in the Jacobian      
    }
  TList.SetNextFreeIndex(kTriplets);
}


// This is almost a selection of lines of LSDerivatives(CccdImge ...)
/* This routine (and the following one) is template because it is used
both with its first argument as "const CCdImage &" and "CcdImage &",
and I did not want to replicate it.  The constness of the iterators is
automagically set by declaring them as "auto" */

template <class ImType, class Accum> 
void AstromFit::AccumulateStatImage(ImType &Ccd, Accum &Accu) const
{
  /*********************************************************************/
  /**  Changes in this routine should be reflected into LSDerivatives  */
  /*********************************************************************/
  /* Setup */
  // 1 : get the Mapping's
  const Mapping *mapping = _distortionModel->GetMapping(Ccd);
  // proper motion stuff
  double jd = Ccd.JD() - _JDRef;
  cout << "in accumulateststimage jd " << jd << endl;
  // refraction stuff
  Point refractionVector = Ccd.ParallacticVector();
  double refractionCoefficient = _refracCoefficient.at(Ccd.BandRank());
  // transformation from sky to TP
  const Gtransfo* sky2TP = _distortionModel->Sky2TP(mapping, Ccd);
  // reserve matrix once for all measurements
  Eigen::Matrix2Xd transW(2,2);

  auto &catalog = Ccd.CatalogForFit();
  for (auto i = catalog.begin(); i!= catalog.end(); ++i)
    {
      auto &ms = **i;
      if (!ms.IsValid()) continue;
      // tweak the measurement errors
      FatPoint inPos = ms;
      TweakAstromMeasurementErrors(inPos, ms);

      FatPoint outPos;
      // should *not* fill h if WhatToFit excludes mapping parameters.
      mapping->TransformPosAndErrors(inPos, outPos);
      double det = outPos.vx*outPos.vy-sqr(outPos.vxy);
      if (det <=0 || outPos.vx <=0 || outPos.vy<=0) {
	cout << " WARNING: inconsistent measurement errors :drop measurement at " << Point(ms) << " in image " << Ccd.Name() << endl;
	continue;
      }	
      transW(0,0) = outPos.vy/det;
      transW(1,1) = outPos.vx/det;
      transW(0,1) = transW(1,0) = -outPos.vxy/det;

      const FittedStar *fs = ms.GetFittedStar();
      Point fittedStarInTP = TransformFittedStar(*fs, sky2TP,
						 refractionVector, 
						 refractionCoefficient,
						 jd);

      Eigen::Vector2d res(fittedStarInTP.x-outPos.x, fittedStarInTP.y-outPos.y); 
      double chi2Val = res.transpose()*transW*res;

      cout << "in accumulateststimage chi2Val " << chi2Val << endl;
      Accu.AddEntry(chi2Val, 2, &ms);
    }// end of loop on measurements
}


//! for a list of images.
template <class ListType, class Accum> 
void AstromFit::AccumulateStatImageList(ListType &L, Accum &Accu) const
{
  for (auto im=L.begin(); im!=L.end() ; ++im)
    {
      AccumulateStatImage(**im, Accu);
    }
}

//! for the list of images in the provided  association and the reference stars, if any
AstromFit::Chi2 AstromFit::ComputeChi2() const
{
  Chi2 chi2;
  AccumulateStatImageList(_assoc.TheCcdImageList(), chi2);
  // Now add the ref stars
  const FittedStarList &fsl = _assoc.fittedStarList;
  /* If you wonder why we project here, read comments in 
     AstromFit::LSDerivatives(TripletList &TList, Eigen::VectorXd &Rhs) */
  TanRaDec2Pix proj(GtransfoLin(), Point(0.,0.));
  for (auto i = fsl.cbegin(); i!= fsl.end(); ++i)
    {
      const FittedStar &fs = **i;
      const RefStar *rs = fs.GetRefStar();
      if (rs == NULL) continue;
      proj.SetTangentPoint(fs);
      // fs projects to (0,0), no need to compute its transform.
      FatPoint rsProj;
      proj.TransformPosAndErrors(*rs, rsProj);
      // TO DO : account for proper motions.
      double rx = rsProj.x; // -fsProj.x (which is 0)
      double ry = rsProj.y;
      double det = rsProj.vx*rsProj.vy - sqr(rsProj.vxy);
      double wxx = rsProj.vy/det;
      double wyy = rsProj.vx/det;
      double wxy = -rsProj.vxy/det;
      wxx *= HACK_REF_ERRORS;
      wyy *= HACK_REF_ERRORS;
      wxy *= HACK_REF_ERRORS;
      chi2.AddEntry(wxx*sqr(rx) + 2*wxy*rx*ry+ wyy*sqr(ry), 2, NULL);
    }
  // so far, ndof contains the number of squares.
  // So, subtract here the number of parameters.
  chi2.ndof -= _nParTot;
  cout << "Finally the chi2 " << chi2 << " " << _nParTot << endl;
  return chi2;
}

//! a class to accumulate chi2 contributions together with pointers to the contributors.
/*! This structure allows to compute the chi2 statistics (average and
  variance) and directly point back to the bad guys without
  relooping. The Chi2Entry routine makes it compatible with
  AccumulateStatImage and AccumulateStatImageList. */
struct Chi2Entry
{
  double chi2;
  MeasuredStar *ms;

  Chi2Entry(const double &c, MeasuredStar *s): chi2(c), ms(s) {}
  // for sort
  bool operator < (const Chi2Entry &R) const {return (chi2<R.chi2);}
};

struct Chi2Vect : public vector<Chi2Entry>
{
  void AddEntry(const double &Chi2Val, unsigned ndof, MeasuredStar *ms)
  { push_back(Chi2Entry(Chi2Val,ms));}

};

//! this routine is to be used only in the framework of outlier removal
/*! it fills the array of indices of parameters that a Measured star
    contrains. Not really all of them if you check. */
void AstromFit::GetMeasuredStarIndices(const MeasuredStar &Ms, 
				       std::vector<unsigned> &Indices) const
{
  Indices.clear();
  if (_fittingDistortions)
    {
      const Mapping *mapping = _distortionModel->GetMapping(*Ms.ccdImage);
      mapping->GetMappingIndices(Indices);
    }
  const FittedStar *fs= Ms.GetFittedStar();
  unsigned fsIndex = fs->IndexInMatrix();
  if (_fittingPos)
    {
      Indices.push_back(fsIndex);
      Indices.push_back(fsIndex+1);
    }
  // For securing the outlier removal, the next block is just useless 
  if (_fittingPM)
     {
       for (unsigned k=0; k<NPAR_PM; ++k) Indices.push_back(fsIndex+2+k);
     }
  /* Should not put the index of refaction stuff or we will not be
     able to remove more than 1 star at a time. */
}

//! Discards measurements contributing more than a cut, computed as <chi2>+NSigCut+rms(chi2). Returns the number of removed outliers. No refit done.
/*! After returning form here, there are still measurements that
  contribute above the cut, but their contribution should be
  evaluated after a refit before discarding them . */
unsigned AstromFit::RemoveOutliers(const double &NSigCut)
{
  /* Some reshuffling would be needed if we wish to use the small-rank
     update trick rather than resolving over again. Typically We would
     need to compute the Jacobian and RHS contributions of the
     discarded measurement and update the current factorization and
     solution. */
  CcdImageList &L=_assoc.ccdImageList;
  // collect chi2 contributions
  Chi2Vect chi2s;
  chi2s.reserve(_nMeasuredStars);
  AccumulateStatImageList(_assoc.ccdImageList, chi2s);
  // do some stat
  unsigned nval = chi2s.size();
  if (nval==0) return 0;
  sort(chi2s.begin(), chi2s.end());
  double median = (nval & 1)? chi2s[nval/2].chi2 :  
    0.5*(chi2s[nval/2-1].chi2 + chi2s[nval/2].chi2);
  // some more stats. should go into the class if recycled anywhere else
  double sum=0; double sum2 = 0;
  for (auto i=chi2s.begin(); i!=chi2s.end(); ++i) 
    {sum+= i->chi2;sum2+= sqr(i->chi2);}
  double average = sum/nval;
  double sigma = sqrt(sum2/nval - sqr(average));
  cout << "INFO : RemoveOutliers chi2 stat: mean/median/sigma " 
       << average << '/'<< median << '/' << sigma << endl;
  double cut = average+NSigCut*sigma;
  /* For each of the parameters, we will not remove more than 1
     measurement that contributes to constraining it. Keep track using
     of what we are touching using an integer vector. This is the
     trick that Marc Betoule came up to for outlier removals in "star
     flats" fits. */
  Eigen::VectorXi affectedParams(_nParTot);
  affectedParams.setZero();

  unsigned removed = 0; // returned to the caller
  // start from the strongest outliers.
  for (auto i = chi2s.rbegin(); i != chi2s.rend(); ++i)
    {
      if (i->chi2 < cut) break; // because the array is sorted. 
      vector<unsigned> indices;
      GetMeasuredStarIndices(*(i->ms), indices);
      bool drop_it = true;
      /* find out is a stronger outlier contraining on the parameters
	 this one contrains was already discarded. If yes, we keep this one */
      for (auto i=indices.cbegin(); i!= indices.end(); ++i)
	if (affectedParams(*i) !=0) drop_it = false;
      
      if (drop_it)
	{
	  FittedStar *fs = i->ms->GetFittedStar();
	  i->ms->SetValid(false); removed++;
	  fs->MeasurementCount()--; // could be put in SetValid
	  /* By making sure that we do not remove all MeasuredStars
	     pointing to a FittedStar in a single go,
	     fs->MeasurementCount() should never go to 0. 
	     
	     It seems plausible that the adopted mechanism prevents as
	     well to end up with under-constrained transfos. */
	  for (auto i=indices.cbegin(); i!= indices.end(); ++i)
	    affectedParams(*i)++;
	}
    } // end loop on measurements
  cout << "INFO : RemoveOutliers : found and removed " 
       << removed << " outliers" << endl;
  return removed;
}

/*! WhatToFit is searched for strings : "Distortions", "Positions",
"Refrac", "PM" which define which parameter set is going to be
variable when computing derivatives (LSDerivatives) and minimizing
(Minimize()).  WhatToFit="Positions Distortions" will minimize w.r.t
mappings and objects positions, and not w.r.t proper motions and
refraction modeling.  However if proper motions and/or refraction
parameters have already been set, then they are accounted for when
computing residuals.  The string is forwarded to the DistortionModel,
and it can then be used to turn subsets of distortion parameter on or
off, if the DistortionModel implements such a thing.
*/
void AstromFit::AssignIndices(const std::string &WhatToFit)
{
  _WhatToFit = WhatToFit;
  cout << "INFO: we are going to fit : " << WhatToFit << endl;
  _fittingDistortions = (_WhatToFit.find("Distortions") != string::npos);
  _fittingPos = (_WhatToFit.find("Positions") != string::npos);
  _fittingRefrac = (_WhatToFit.find("Refrac") != string::npos);
  if (_sigCol == 0 && _fittingRefrac) 
    {
      cout << "WARNING: We cannot fit refraction coefficients without a color lever arm. Ignoring refraction" << endl;
      _fittingRefrac = false;
    }	
  _fittingPM = (_WhatToFit.find("PM") != string::npos);
// When entering here, we assume that WhatToFit has already been interpreted.


  _nParDistortions = 0;
  if (_fittingDistortions) 
    {
      _nParDistortions = _distortionModel->AssignIndices(0,_WhatToFit);
      cout << "fitting distorsions " << _nParDistortions << endl;
    }
      
  unsigned ipar = _nParDistortions;

  if (_fittingPos)
    {
      FittedStarList &fsl = _assoc.fittedStarList;
      for (FittedStarIterator i= fsl.begin(); i != fsl.end(); ++i)
	{
	  FittedStar &fs = **i;
	  // the parameter layout here is used also
	  // - when filling the derivatives
	  // - when updating (OffsetParams())
	  // - in GetMeasuredStarIndices
	  fs.SetIndexInMatrix(ipar);
	  ipar+=2;
	  if ((_fittingPM) & fs.mightMove) ipar+= NPAR_PM;
	}
    }
  unsigned _nParPositions = ipar-_nParDistortions;
  if (_fittingRefrac)
    { 
      _refracPosInMatrix = ipar;
      ipar += _nRefrac;
    }
  _nParTot = ipar;

  cout << "ici npartot " << _nParTot << std::endl;
#if (0)  
  //DEBUG
  cout << " INFO: np(d,p, total) = " 
       <<  _nParDistortions << ' '
       << _nParPositions << ' '  
       << _nParTot << ' '
       << WhatToFit << endl;
  const FittedStar &ffs = (**(_assoc.fittedStarList.begin()));
  cout << " INFO : first Star Index : " <<  ffs.IndexInMatrix() << ' ' << Point(ffs) << endl;
#endif
}

void AstromFit::OffsetParams(const Eigen::VectorXd& Delta)
{
  if (Delta.size() != _nParTot) 
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "AstromFit::OffsetParams : the provided vector length is not compatible with the current WhatToFit setting");

  if (_fittingDistortions) 
    _distortionModel->OffsetParams(Delta);

  if (_fittingPos)
    {
      FittedStarList &fsl = _assoc.fittedStarList;
      for (FittedStarIterator i= fsl.begin(); i != fsl.end(); ++i)
	{
	  FittedStar &fs = **i;
	  // the parameter layout here is used also
	  // - when filling the derivatives
	  // - when assigning indices (AssignIndices())
	  unsigned index = fs.IndexInMatrix();
	  fs.x += Delta(index);
	  fs.y += Delta(index+1);
	  if ((_fittingPM) & fs.mightMove)
	    {
	      fs.pmx += Delta(index+2);
	      fs.pmy += Delta(index+3);
	    }
	}
    }
  if (_fittingRefrac)
    {
      for (unsigned k=0; k<_nRefrac; ++k)
	_refracCoefficient[k] += Delta(_refracPosInMatrix+k);
    }
}

// should not be too large !

static void write_sparse_matrix_in_fits(const SpMat &mat, const string &FitsName)
{
  if (mat.rows()*mat.cols() > 2e8)
    {
      cout << "WARNING :  write_sparse_matrix_in_fits : yout matrix is too large. " << FitsName << " not generated"<< endl;
      return;
    }
  Mat m(mat.rows(),mat.cols());
  for (int k=0; k<mat.outerSize(); ++k)
    for (SpMat::InnerIterator it(mat,k); it; ++it)
      {
	m (it.row(), it.col()) = it.value();
      }
  m.writeFits(FitsName);
}

static void write_vect_in_fits(const Eigen::VectorXd &V, const string &FitsName)
{
  Vect v(V.size());
  for (int k=0; k <V.size(); ++k) v(k) = V(k);
  Mat(v).writeFits(FitsName);
}




/*! This is a complete Newton Raphson step. Compute first and 
  second derivatives, solve for the step and apply it, without 
  a line search. */
bool AstromFit::Minimize(const std::string &WhatToFit)
{
  AssignIndices(WhatToFit);
  const CcdImageList &ccdImageList = _assoc.TheCcdImageList();
  // Count measurements. unused at the moment ....
  _nMeasuredStars = 0;
  for (auto i = ccdImageList.begin(); i!= ccdImageList.end() ; ++i)
    _nMeasuredStars += (*i)->CatalogForFit().size();

  // TODO : write a guesser for the number of triplets
  unsigned nTrip = (_LastNTrip) ? _LastNTrip: 1e6; // GuessNTrip(ccdImageList);
  TripletList tList(nTrip);
  Eigen::VectorXd grad(_nParTot);  grad.setZero();

  //Fill the triplets
  clock_t tstart = clock();
  LSDerivatives(tList, grad);
  clock_t tend = clock();
  _LastNTrip = tList.size(); 

  cout << " INFO: End of triplet filling, ntrip = " << tList.size() 
       << " CPU = " << float(tend-tstart)/float(CLOCKS_PER_SEC) 
       << endl;

  SpMat hessian;
  {
#if (TRIPLET_INTERNAL_COORD == COL)
    cout << "filling jacobian " << _nParTot << " " << tList.NextFreeIndex() << endl;
    SpMat jacobian(_nParTot,tList.NextFreeIndex());
    jacobian.setFromTriplets(tList.begin(), tList.end());
    // release memory shrink_to_fit is C++11
    cout << jacobian << endl;
    tList.clear(); 
    //tList.shrink_to_fit();
    clock_t tstart = clock();
    hessian = jacobian*jacobian.transpose();
    clock_t tend = clock();
    std::cout << "INFO: CPU for J*Jt " 
	      << float(tend-tstart)/float(CLOCKS_PER_SEC) << std::endl;

#else
    SpMat jacobian(tList.NextRank(), _nParTot);
    jacobian.setFromTriplets(tList.begin(), tList.end());
    // release memory shrink_to_fit is C++11
    tList.clear(); 
    //tList.shrink_to_fit(); 
    cout << " starting H=JtJ " << endl;
    hessian = jacobian.transpose()*jacobian;
#endif
  }// release the Jacobian

  //  write_sparse_matrix_in_fits(hessian, "h.fits");
  cout << "INFO: hessian : dim=" << hessian.rows() << " " << hessian.cols() 
       << " nnz=" << hessian.nonZeros() 
       << " filling-frac = " << hessian.nonZeros()/sqr(hessian.rows()) << endl;
  cout << "INFO: starting factorization" << endl;

  tstart = clock();
  Eigen::SimplicialLDLT<SpMat> chol(hessian);
  if (chol.info() != Eigen::Success)
    {
      cout << "ERROR: AstromFit::Minimize : factorization failed " << endl;
      return false;
    }

  cout << "grad size " << grad.rows() << " " << grad.cols() << endl;
  Eigen::VectorXd delta = chol.solve(grad);

  //  cout << " offsetting parameters" << endl;
  OffsetParams(delta);
  tend = clock();
  std::cout << "INFO: CPU for factor-solve-update " 
  	    << float(tend-tstart)/float(CLOCKS_PER_SEC) << std::endl;
}


void AstromFit::CheckStuff()
{
#if (0)
  const char *what2fit[] = {"Positions", "Distortions", "Refrac",
		      "Positions Distortions", "Positions Refrac", 
		      "Distortions Refrac",
		      "Positions Distortions Refrac"};
#endif
  const char *what2fit[] = {"Positions", "Distortions",
		      "Positions Distortions"};
  // DEBUG
  for (int k=0; k < sizeof(what2fit)/sizeof(what2fit[0]); ++k)
    {
      AssignIndices(what2fit[k]);
#if (0)
      fsIndexDebug = _nParDistortions;
      heavyDebug = true;
#endif
      TripletList tList(10000);
      Eigen::VectorXd rhs(_nParTot);  rhs.setZero();		    
      LSDerivatives(tList, rhs);
      SpMat jacobian(_nParTot,tList.NextFreeIndex());
      jacobian.setFromTriplets(tList.begin(), tList.end());
      SpMat hessian = jacobian*jacobian.transpose();
      
      char name[24];
      sprintf(name,"h%d.fits", k);
      write_sparse_matrix_in_fits(hessian, name);
      sprintf(name,"g%d.fits", k);
      write_vect_in_fits(rhs, name);
      cout << "npar : " << _nParTot << ' ' << _nParDistortions << ' ' << endl;

    }
}

void AstromFit::MakeResTuple(const std::string &TupleName) const
{
  std::ofstream tuple(TupleName.c_str());
  tuple << "#xccd: coordinate in CCD" << endl
	<< "#yccd: " << endl
	<< "#rx:   residual in degrees in TP" << endl
	<< "#ry:" << endl
	<< "#xtp: transformed coordinate in TP " << endl
	<< "#ytp:" << endl  
	<< "#mag: rough mag" << endl
	<< "#jd: Julian date of the measurement" << endl
    	<< "#rvx: transformed measurement uncertainty " << endl
    	<< "#rvy:" << endl
    	<< "#rvxy:" << endl
	<< "#color : " << endl
	<< "#chi2: contribution to Chi2 (2D dofs)" << endl
	<< "#nm: number of measurements of this FittedStar" << endl
    	<< "#chip: chip number" << endl
    	<< "#shoot: shoot id" << endl
	<< "#end" << endl;
  const CcdImageList &L=_assoc.TheCcdImageList();
  for (auto i=L.cbegin(); i!=L.end() ; ++i)
    {
      const CcdImage &im = **i;
      const MeasuredStarList &cat = im.CatalogForFit();
      const Mapping *mapping = _distortionModel->GetMapping(im);
      const Point &refractionVector = im.ParallacticVector();
      double jd = im.JD() - _JDRef;
      unsigned iband= im.BandIndex();
      for (auto is=cat.cbegin(); is!=cat.end(); ++is)
	{
	  const MeasuredStar &ms = **is;
	  FatPoint tpPos;
	  mapping->TransformPosAndErrors(ms, tpPos);
	  const Gtransfo* sky2TP = _distortionModel->Sky2TP(mapping, im);
	  const FittedStar *fs = ms.GetFittedStar();
	  
	  Point fittedStarInTP = TransformFittedStar(*fs, sky2TP,
						     refractionVector, 
						     _refracCoefficient[iband],
						     jd);
	  Point res=tpPos-fittedStarInTP;
	  double det = tpPos.vx*tpPos.vy-sqr(tpPos.vxy);
	  double wxx = tpPos.vy/det;
	  double wyy = tpPos.vx/det;
	  double wxy = -tpPos.vxy/det;
	  //	  double chi2 = rx*(wxx*rx+wxy*ry)+ry*(wxy*rx+wyy*ry);
	  double chi2 = wxx*res.x*res.x + wyy*res.y*res.y + 2*wxy*res.x*res.y;
	  tuple << ms.x << ' ' << ms.y << ' ' 
		<< res.x << ' ' << res.y << ' '
		<< tpPos.x << ' ' << tpPos.y << ' '
		<< fs->Mag() << ' ' << jd << ' ' 
		<< tpPos.vx << ' ' << tpPos.vy << ' ' << tpPos.vxy << ' ' 
		<< fs->color << ' ' 
		<< chi2 << ' ' 
		<< fs->MeasurementCount() << ' ' 
		<< im.Chip() << ' ' << im.Shoot() << endl;
	}// loop on measurements in image
    }// loop on images

}

}}}