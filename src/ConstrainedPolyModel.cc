#include "lsst/log/Log.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/ConstrainedPolyModel.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/AstrometryModel.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/ProjectionHandler.h"
#include "lsst/jointcal/AstroUtils.h" // ApplyTransfo(Frame)

#include "lsst/pex/exceptions.h"
namespace pexExcept = lsst::pex::exceptions;

#include <string>
#include <iostream>

namespace {
    LOG_LOGGER _log = LOG_GET("jointcal.ConstrainedPolyModel");
}

namespace lsst {
namespace jointcal {

/* This code does not contain anything involved. It just maps the
routines AstrometryFit needs to what is needed for this two-transfo model.
The two-transfo mappings are implemented using two one-transfo
mappings.*/

// TODO : separate the polynomial degrees for chip and visit transfos.
// TODO propagate those into python:
  static int DistortionDegree = 3;

using namespace std;

ConstrainedPolyModel::ConstrainedPolyModel(const CcdImageList &ccdImageList,
                                           const ProjectionHandler* projectionHandler,
                                           bool initFromWCS,
                                           unsigned nNotFit) : _sky2TP(projectionHandler)

{
  // from datacards (or default)
  unsigned degree = DistortionDegree;
  unsigned count = 0;
  VisitIdType refVisit;
  // first loop to initialize all visit  and chip transfos.
  for (auto i=ccdImageList.cbegin(); i!= ccdImageList.cend(); ++i, ++count)
    {
      const CcdImage &im = **i;
      unsigned visit = im.getVisit();
      unsigned chip = im.getCcdId();
      auto visitp = _visitMap.find(visit);
      if (visitp == _visitMap.end())
	{
	  if (_visitMap.size() == 0)
	    {
	      // if one fits all of them, the model is degenerate.
#ifdef ROTATE_T2
# warning : hack in ConstrainedPolyModel::ConstrainedPolyModel : rotated frame
	      _visitMap[visit] = new SimpleGtransfoMapping(GtransfoLinRot(3.141927/2.), /* ToFit = */ false);
#else
	      _visitMap[visit] = std::unique_ptr<SimpleGtransfoMapping>(new SimpleGtransfoMapping(GtransfoIdentity()));
#endif
	      refVisit = visit;
	    }
	    else
#ifdef ROTATE_T2
	      {
		GtransfoPoly poly(degree);
		poly = GtransfoPoly(GtransfoLinRot(3.141927/2.))*poly;
		_visitMap[visit] = new SimplePolyMapping(GtransfoLin(), poly);
	      }
#else
	  _visitMap[visit] = std::unique_ptr<SimplePolyMapping>(new SimplePolyMapping(GtransfoLin(),
										      GtransfoPoly(degree)));
#endif
	}
      auto chipp = _chipMap.find(chip);
      if ((chipp == _chipMap.end()) && visit == refVisit )
	{
	  const Frame &frame = im.getImageFrame();

	  _tpFrame += ApplyTransfo(frame, *im.Pix2CommonTangentPlane(), LargeFrame);
	  GtransfoPoly pol(im.Pix2TangentPlane(),
			   frame,
			   degree);
	  GtransfoLin shiftAndNormalize = NormalizeCoordinatesTransfo(frame);

	  _chipMap[chip] = std::unique_ptr<SimplePolyMapping>(new SimplePolyMapping(shiftAndNormalize, pol*shiftAndNormalize.invert()));
	}
    }
  // now, second loop to set the mappings of the CCdImages
  for (auto i=ccdImageList.cbegin(); i!= ccdImageList.cend(); ++i, ++count)
    {
      const CcdImage &im = **i;
      unsigned visit = im.getVisit();
      unsigned chip = im.getCcdId();
      // check that the chip_indexed part was indeed assigned
      // (i.e. the reference visit was complete)
      if (_chipMap.find(chip) == _chipMap.end())
	{
        LOGLS_WARN(_log, "Chip " << chip << " is missing in the reference exposure, expect troubles.");
	  GtransfoLin norm = NormalizeCoordinatesTransfo(im.getImageFrame());
	  _chipMap[chip] = std::unique_ptr<SimplePolyMapping>( new SimplePolyMapping(norm,
										     GtransfoPoly(degree)));
	}
      _mappings[&im] = std::unique_ptr<TwoTransfoMapping>(new TwoTransfoMapping(_chipMap[chip].get(), _visitMap[visit].get()));

    }
  LOGLS_INFO(_log, "Constructor got " << _chipMap.size() << " chip mappings and "
             << _visitMap.size() << " visit mappings.");
  // DEBUG
  for (auto i=_visitMap.begin(); i != _visitMap.end(); ++i)
    LOGLS_DEBUG(_log, i->first);
}

const Mapping* ConstrainedPolyModel::getMapping(const CcdImage &C) const
{
  mappingMapType::const_iterator i = _mappings.find(&C);
  if  (i==_mappings.end()) return nullptr;
  return (i->second.get());
}

/*! This routine decodes "DistortionsChip" and "DistortionsVisit" in
  WhatToFit. If WhatToFit contains "Distortions" and not
  Distortions<Something>, it is understood as both chips and
  visits. */
unsigned ConstrainedPolyModel::assignIndices(unsigned FirstIndex,
					     std::string &WhatToFit)
{
  unsigned index=FirstIndex;
  if (WhatToFit.find("Distortions") == std::string::npos)
    {
        LOGLS_ERROR(_log, "assignIndices was called and Distortions is *not* in WhatToFit");
        return 0;
    }
  // if we get here "Distortions" is in WhatToFit
  _fittingChips = (WhatToFit.find("DistortionsChip") != std::string::npos);
  _fittingVisits = (WhatToFit.find("DistortionsVisit") != std::string::npos);
  // If nothing more than "Distortions" is specified, it means all:
  if ((!_fittingChips)&&(!_fittingVisits))
    {_fittingChips = _fittingVisits = true;}
  if (_fittingChips)
    for (auto i = _chipMap.begin(); i!=_chipMap.end(); ++i)
      {
	SimplePolyMapping *p = dynamic_cast<SimplePolyMapping *>(&*(i->second));
	if (!p)
	  throw LSST_EXCEPT(pexExcept::InvalidParameterError,"ERROR: in ConstrainedPolyModel, all chip \
mappings should be SimplePolyMappings");
	p->SetIndex(index);
	index+= p->Npar();
      }
  if (_fittingVisits)
    for (auto i = _visitMap.begin(); i!=_visitMap.end(); ++i)
      {
	SimplePolyMapping *p = dynamic_cast<SimplePolyMapping *>(&*(i->second));
	if (!p) continue; // it should be GtransfoIdentity
	p->SetIndex(index);
	index+= p->Npar();
      }
  // Tell the mappings which derivatives they will have to fill:
  for (auto i = _mappings.begin(); i != _mappings.end() ; ++i)
    {
      i->second->SetWhatToFit(_fittingChips, _fittingVisits);
    }
  return index;
}

void ConstrainedPolyModel::offsetParams(const Eigen::VectorXd &Delta)
{
  if (_fittingChips)
    for (auto i = _chipMap.begin(); i!=_chipMap.end(); ++i)
      {
        auto *p = (&*(i->second));
        if (p->Npar()) // probably useless test
	  p->OffsetParams(&Delta(p->Index()));
      }
  if (_fittingVisits)
    for (auto i = _visitMap.begin(); i!=_visitMap.end(); ++i)
      {
        auto *p = (&*(i->second));
        if (p->Npar()) // probably useless test
	  p->OffsetParams(&Delta(p->Index()));
      }
}

void ConstrainedPolyModel::freezeErrorScales()
{
  for (auto i = _visitMap.begin(); i!=_visitMap.end(); ++i)
    i->second->FreezeErrorScales();
  for (auto i = _chipMap.begin(); i!=_chipMap.end(); ++i)
    i->second->FreezeErrorScales();
}


const Gtransfo& ConstrainedPolyModel::getChipTransfo(const unsigned Chip) const
{
  auto chipp = _chipMap.find(Chip);
  if (chipp == _chipMap.end()) {
    std::stringstream errMsg;
    errMsg << "No such chipId: '" << Chip << "' found in chipMap of:  " << this;
    throw pexExcept::InvalidParameterError(errMsg.str());
  }
  return chipp->second->Transfo();
}

// Array of visits involved in the solution.
std::vector<VisitIdType> ConstrainedPolyModel::getVisits() const
{
  std::vector<VisitIdType> res;
  res.reserve(_visitMap.size());
  for (auto i = _visitMap.begin(); i!=_visitMap.end(); ++i)
    res.push_back(i->first);
  return res;
}

const Gtransfo& ConstrainedPolyModel::getVisitTransfo(const VisitIdType &Visit) const
{
  auto visitp = _visitMap.find(Visit);
  if (visitp == _visitMap.end()) {
    std::stringstream errMsg;
    errMsg << "No such visitId: '" << Visit << "' found in visitMap of: " << this;
    throw pexExcept::InvalidParameterError(errMsg.str());
  }
  return visitp->second->Transfo();
}


std::shared_ptr<TanSipPix2RaDec> ConstrainedPolyModel::produceSipWcs(const CcdImage &ccdImage) const
{
  mappingMapType::const_iterator i = _mappings.find(&ccdImage);
  if  (i==_mappings.end()) return nullptr;
  const TwoTransfoMapping *m = i->second.get();

  const GtransfoPoly &t1=dynamic_cast<const GtransfoPoly&>(m->T1());
  const GtransfoPoly &t2=dynamic_cast<const GtransfoPoly&>(m->T2());
  const TanRaDec2Pix *proj=dynamic_cast<const TanRaDec2Pix*>(sky2TP(ccdImage));
  if (!(&t1)  || !(&t2) || !proj) return nullptr;

  GtransfoPoly pix2Tp = t2*t1;

  const GtransfoLin &projLinPart = proj->LinPart(); // should be the identity, but who knows? So, let us incorporate it into the pix2TP part.
  GtransfoPoly wcsPix2Tp = GtransfoPoly(projLinPart.invert())*pix2Tp;

  // compute a decent approximation, if higher order corrections get ignored
  GtransfoLin cdStuff = wcsPix2Tp.LinearApproximation(ccdImage.getImageFrame().Center());

  // wcsPix2TP = cdStuff*sip , so
  GtransfoPoly sip = GtransfoPoly(cdStuff.invert())*wcsPix2Tp;
  Point tangentPoint( proj->TangentPoint());
  return std::shared_ptr<TanSipPix2RaDec>(new TanSipPix2RaDec(cdStuff, tangentPoint, &sip));
}

}} // end of namespaces
