// -*- LSST-C++ -*-
/*
 * This file is part of jointcal.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <iomanip>
#include <iterator> /* for ostream_iterator */
#include <limits>
#include <cmath>
#include <math.h>
#include <fstream>
#include "assert.h"
#include <sstream>

#include "Eigen/Core"

#include "lsst/log/Log.h"
#include "lsst/geom/Point.h"
#include "lsst/jointcal/AstrometryTransform.h"
#include "lsst/jointcal/Frame.h"
#include "lsst/jointcal/StarMatch.h"
#include "lsst/pex/exceptions.h"
#include "Eigen/Cholesky"

using long_double=double;

namespace pexExcept = lsst::pex::exceptions;

using namespace std;

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.AstrometryTransform");
}

namespace lsst {
namespace jointcal {

bool isIntegerShift(const AstrometryTransform *transform) {
    const AstrometryTransformPolynomial *shift =
            dynamic_cast<const AstrometryTransformPolynomial *>(transform);
    if (shift == nullptr) return false;

    static const double eps = 1e-5;

    double dx = shift->getCoefficient(0, 0, 0);
    double dy = shift->getCoefficient(0, 0, 1);

    static Point dumb(4000, 4000);
    if (fabs(dx - int(floor(dx + 0.5))) < eps && fabs(dy - int(floor(dy + 0.5))) < eps &&
        fabs(dumb.x + dx - shift->apply(dumb).x) < eps && fabs(dumb.y + dy - shift->apply(dumb).y) < eps)
        return true;

    return false;
}

/********* AstrometryTransform ***********************/

Frame AstrometryTransform::apply(Frame const &inputframe, bool inscribed) const {
    // 2 opposite corners
    double xtmin1, xtmax1, ytmin1, ytmax1;
    apply(inputframe.xMin, inputframe.yMin, xtmin1, ytmin1);
    apply(inputframe.xMax, inputframe.yMax, xtmax1, ytmax1);
    Frame fr1(std::min(xtmin1, xtmax1), std::min(ytmin1, ytmax1), std::max(xtmin1, xtmax1),
              std::max(ytmin1, ytmax1));
    // 2 other corners
    double xtmin2, xtmax2, ytmin2, ytmax2;
    apply(inputframe.xMin, inputframe.yMax, xtmin2, ytmax2);
    apply(inputframe.xMax, inputframe.yMin, xtmax2, ytmin2);
    Frame fr2(std::min(xtmin2, xtmax2), std::min(ytmin2, ytmax2), std::max(xtmin2, xtmax2),
              std::max(ytmin2, ytmax2));

    if (inscribed) return fr1 * fr2;
    return fr1 + fr2;
}

std::unique_ptr<AstrometryTransform> AstrometryTransform::composeAndReduce(
        AstrometryTransform const &) const {  // by default no way to compose
    return std::unique_ptr<AstrometryTransform>(nullptr);
}

double AstrometryTransform::getJacobian(const double x, const double y) const {
    double x2, y2;
    double eps = x * 0.01;
    if (eps == 0) eps = 0.01;
    apply(x, y, x2, y2);
    double dxdx, dydx;
    apply(x + eps, y, dxdx, dydx);
    dxdx -= x2;
    dydx -= y2;
    double dxdy, dydy;
    apply(x, y + eps, dxdy, dydy);
    dxdy -= x2;
    dydy -= y2;
    return ((dxdx * dydy - dxdy * dydx) / (eps * eps));
}

/*! the Derivative is represented by a AstrometryTransformLinear, in which
  (hopefully), the offset terms are zero. Derivative should
  transform a vector of offsets into a vector of offsets. */
void AstrometryTransform::computeDerivative(Point const &where, AstrometryTransformLinear &derivative,
                                            const double step) const {
    double x = where.x;
    double y = where.y;
    double xp0, yp0;
    apply(x, y, xp0, yp0);

    double xp, yp;
    apply(x + step, y, xp, yp);
    derivative.a11() = (xp - xp0) / step;
    derivative.a21() = (yp - yp0) / step;
    apply(x, y + step, xp, yp);
    derivative.a12() = (xp - xp0) / step;
    derivative.a22() = (yp - yp0) / step;
    derivative.dx() = 0;
    derivative.dy() = 0;
}

AstrometryTransformLinear AstrometryTransform::linearApproximation(Point const &where,
                                                                   const double step) const {
    Point outwhere = apply(where);
    AstrometryTransformLinear der;
    computeDerivative(where, der, step);
    return AstrometryTransformLinearShift(outwhere.x, outwhere.y) * der *
           AstrometryTransformLinearShift(-where.x, -where.y);
}

void AstrometryTransform::transformPosAndErrors(FatPoint const &in, FatPoint &out) const {
    FatPoint res;  // in case in and out are the same address...
    res = apply(in);
    AstrometryTransformLinear der;
    // could save a call here, since Derivative needs the transform of where that we already have
    // 0.01 may not be a very good idea in all cases. May be we should provide a way of altering that.
    computeDerivative(in, der, 0.01);
    double a11 = der.A11();
    double a22 = der.A22();
    double a21 = der.A21();
    double a12 = der.A12();
    res.vx = a11 * (a11 * in.vx + 2 * a12 * in.vxy) + a12 * a12 * in.vy;
    res.vy = a21 * a21 * in.vx + a22 * a22 * in.vy + 2. * a21 * a22 * in.vxy;
    res.vxy = a21 * a11 * in.vx + a22 * a12 * in.vy + (a21 * a12 + a11 * a22) * in.vxy;
    out = res;
}

void AstrometryTransform::transformErrors(Point const &where, const double *vIn, double *vOut) const {
    AstrometryTransformLinear der;
    computeDerivative(where, der, 0.01);
    double a11 = der.A11();
    double a22 = der.A22();
    double a21 = der.A21();
    double a12 = der.A12();

    /*      (a11 a12)          (vxx  vxy)
       M =  (       )  and V = (        )
            (a21 a22)          (xvy  vyy)

       Vxx = Vin[0], vyy = Vin[1], Vxy = Vin[2];
       we want to compute M*V*tp(M)
       A lin alg light package would be perfect...
    */
    int xx = 0;
    int yy = 1;
    int xy = 2;
    // M*V :

    double b11 = a11 * vIn[xx] + a12 * vIn[xy];
    double b22 = a21 * vIn[xy] + a22 * vIn[yy];
    double b12 = a11 * vIn[xy] + a12 * vIn[yy];
    double b21 = a21 * vIn[xx] + a22 * vIn[xy];

    // (M*V) * tp(M)

    vOut[xx] = b11 * a11 + b12 * a12;
    vOut[xy] = b11 * a21 + b12 * a22;
    vOut[yy] = b21 * a21 + b22 * a22;
}

std::unique_ptr<AstrometryTransform> AstrometryTransform::roughInverse(const Frame &region) const {
    // "in" and "out" refer to the inverse direction.
    Point centerOut = region.getCenter();
    Point centerIn = apply(centerOut);
    AstrometryTransformLinear der;
    computeDerivative(centerOut, der, std::sqrt(region.getArea()) / 5.);
    der = der.inverted();
    der = AstrometryTransformLinearShift(centerOut.x, centerOut.y) * der *
          AstrometryTransformLinearShift(-centerIn.x, -centerIn.y);
    return std::unique_ptr<AstrometryTransform>(new AstrometryTransformLinear(der));
}

/* implement one in AstrometryTransform, so that all derived
   classes do not need to provide one... */

/* the routines that follow are used for ea generic parameter
   transformation serialization, used e.g. for fits. Enables
   to manipulate transformation parameters as vectors.
*/

// not dummy : what it does is virtual because paramRef is virtual.
void AstrometryTransform::getParams(double *params) const {
    std::size_t npar = getNpar();
    for (std::size_t i = 0; i < npar; ++i) params[i] = paramRef(i);
}

void AstrometryTransform::offsetParams(Eigen::VectorXd const &delta) {
    std::size_t npar = getNpar();
    for (std::size_t i = 0; i < npar; ++i) paramRef(i) += delta[i];
}

double AstrometryTransform::paramRef(Eigen::Index const) const {
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                      std::string("AstrometryTransform::paramRef should never be called "));
}

double &AstrometryTransform::paramRef(Eigen::Index const) {
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                      "AstrometryTransform::paramRef should never be called ");
}

void AstrometryTransform::paramDerivatives(Point const &, double *, double *) const {
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                      "AstrometryTransform::paramDerivatives() should never be called ");
}

ostream &operator<<(ostream &stream, AstrometryTransform const &transform) {
    transform.print(stream);
    return stream;
}

void AstrometryTransform::write(const std::string &fileName) const {
    ofstream s(fileName.c_str());
    write(s);
    bool ok = !s.fail();
    s.close();
    if (!ok)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "AstrometryTransform::write, something went wrong for file " + fileName);
}

void AstrometryTransform::write(ostream &stream) const {
    throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            "AstrometryTransform::write(ostream), should never be called. MEans that it is missing in some "
            "derived class ");
}

/******************* GTransformInverse ****************/
/* inverse transformation, solved by iterations. Before using
   it (probably via AstrometryTransform::inverseTransform), consider
   seriously StarMatchList::inverseTransform */
class AstrometryTransformInverse : public AstrometryTransform {
private:
    std::unique_ptr<AstrometryTransform> _direct;
    std::unique_ptr<AstrometryTransform> _roughInverse;
    double precision2;

public:
    AstrometryTransformInverse(const AstrometryTransform *direct, const double precision,
                               const Frame &region);

    //! implements an iterative (Gauss-Newton) solver. It resorts to the Derivative function: 4 calls to the
    //! direct transform per iteration.
    void apply(const double xIn, const double yIn, double &xOut, double &yOut) const;

    void print(ostream &out) const;

    double fit(StarMatchList const &starMatchList);

    virtual std::unique_ptr<AstrometryTransform> clone() const;

    AstrometryTransformInverse(AstrometryTransformInverse const &);

    //! Overload the "generic routine"
    std::unique_ptr<AstrometryTransform> roughInverse(const Frame &) const { return _direct->clone(); }

    //! Inverse transform: returns the direct one!
    std::unique_ptr<AstrometryTransform> inverseTransform(double, const Frame &) const {
        return _direct->clone();
    }

    ~AstrometryTransformInverse();

private:
    void operator=(AstrometryTransformInverse const &);
};

std::unique_ptr<AstrometryTransform> AstrometryTransform::inverseTransform(const double precision,
                                                                           const Frame &region) const {
    return std::unique_ptr<AstrometryTransform>(new AstrometryTransformInverse(this, precision, region));
}

AstrometryTransformInverse::AstrometryTransformInverse(const AstrometryTransform *direct,
                                                       const double precision, const Frame &region) {
    _direct = direct->clone();
    _roughInverse = _direct->roughInverse(region);
    precision2 = precision * precision;
}

AstrometryTransformInverse::AstrometryTransformInverse(AstrometryTransformInverse const &model)
        : AstrometryTransform() {
    _direct = model._direct->clone();
    _roughInverse = model._roughInverse->clone();
    precision2 = model.precision2;
}

AstrometryTransformInverse::~AstrometryTransformInverse() {}

void AstrometryTransformInverse::operator=(AstrometryTransformInverse const &model) {
    _direct = model._direct->clone();
    _roughInverse = model._roughInverse->clone();
    precision2 = model.precision2;
}

void AstrometryTransformInverse::apply(const double xIn, const double yIn, double &xOut, double &yOut) const {
    Point in(xIn, yIn);
    Point outGuess = _roughInverse->apply(in);
    AstrometryTransformLinear directDer, reverseDer;
    int loop = 0;
    int maxloop = 20;
    double move2;
    do {
        loop++;
        Point inGuess = _direct->apply(outGuess);
        _direct->computeDerivative(outGuess, directDer);
        reverseDer = directDer.inverted();
        double xShift, yShift;
        reverseDer.apply(xIn - inGuess.x, yIn - inGuess.y, xShift, yShift);
        outGuess.x += xShift;
        outGuess.y += yShift;
        move2 = xShift * xShift + yShift * yShift;
    } while ((move2 > precision2) && (loop < maxloop));
    if (loop == maxloop) LOGLS_WARN(_log, "Problems applying AstrometryTransformInverse at " << in);
    xOut = outGuess.x;
    yOut = outGuess.y;
}

void AstrometryTransformInverse::print(ostream &stream) const {
    stream << " AstrometryTransformInverse of  :" << endl << *_direct << endl;
}

double AstrometryTransformInverse::fit(StarMatchList const &) {
    throw pexExcept::RuntimeError(
            "Cannot fit a AstrometryTransformInverse. Use StarMatchList::inverseTransform instead.");
}

std::unique_ptr<AstrometryTransform> AstrometryTransformInverse::clone() const {
    return std::unique_ptr<AstrometryTransform>(new AstrometryTransformInverse(*this));
}

/************* AstrometryTransformComposition **************/

// This class was done to allow composition of AstrometryTransform's, without specifications of their types.
// does not need to be public. Invoked  by compose(left,right)

//! Private class to handle AstrometryTransform compositions (i.e. piping). Use the routine compose if
//! you need this functionnality.
class AstrometryTransformComposition : public AstrometryTransform {
private:
    std::unique_ptr<AstrometryTransform> _first, _second;

public:
    AstrometryTransformComposition(AstrometryTransform const &second, AstrometryTransform const &first);

    //! return second(first(xIn,yIn))
    void apply(const double xIn, const double yIn, double &xOut, double &yOut) const;
    void print(ostream &stream) const;

    //!
    double fit(StarMatchList const &starMatchList);

    std::unique_ptr<AstrometryTransform> clone() const;
    ~AstrometryTransformComposition();
};

AstrometryTransformComposition::AstrometryTransformComposition(AstrometryTransform const &second,
                                                               AstrometryTransform const &first) {
    _first = first.clone();
    _second = second.clone();
}

void AstrometryTransformComposition::apply(const double xIn, const double yIn, double &xOut,
                                           double &yOut) const {
    double xout, yout;
    _first->apply(xIn, yIn, xout, yout);
    _second->apply(xout, yout, xOut, yOut);
}

void AstrometryTransformComposition::print(ostream &stream) const {
    stream << "Composed AstrometryTransform consisting of:" << std::endl;
    _first->print(stream);
    _second->print(stream);
}

double AstrometryTransformComposition::fit(StarMatchList const &starMatchList) {
    /* fits only one of them. could check that first can actually be fitted... */
    return _first->fit(starMatchList);
}

std::unique_ptr<AstrometryTransform> AstrometryTransformComposition::clone() const {
    return std::make_unique<AstrometryTransformComposition>(*_second, *_first);
}

AstrometryTransformComposition::~AstrometryTransformComposition() {}

std::unique_ptr<AstrometryTransform> compose(AstrometryTransform const &left,
                                             AstrometryTransformIdentity const &right) {
    return left.clone();
}

std::unique_ptr<AstrometryTransform> compose(AstrometryTransform const &left,
                                             AstrometryTransform const &right) {
    // Try to use the composeAndReduce method from left. If absent, AstrometryTransform::composeAndReduce
    // returns NULL. composeAndReduce is non trivial for polynomials.
    std::unique_ptr<AstrometryTransform> composition(left.composeAndReduce(right));
    // composition == NULL means no reduction: just build a Composition that pipelines "left" and "right".
    if (composition == nullptr)
        return std::make_unique<AstrometryTransformComposition>(left, right);
    else
        return composition;
}

// just a speed up, to avoid useless numerical derivation.
void AstrometryTransformIdentity::computeDerivative(Point const &, AstrometryTransformLinear &derivative,
                                                    const double) const {
    derivative = AstrometryTransformLinear();
}

AstrometryTransformLinear AstrometryTransformIdentity::linearApproximation(Point const &,
                                                                           const double) const {
    AstrometryTransformLinear result;
    return result;  // rely on default AstrometryTransformlin constructor;
}

std::shared_ptr<ast::Mapping> AstrometryTransformIdentity::toAstMap(jointcal::Frame const &) const {
    return std::make_shared<ast::UnitMap>(2);  // a AstrometryTransformIdentity is identically ast::UnitMap(2)
}

void AstrometryTransformIdentity::write(ostream &stream) const {
    stream << "AstrometryTransformIdentity 1" << endl;
}

void AstrometryTransformIdentity::read(istream &stream) {
    int format;
    stream >> format;
    if (format != 1)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          " AstrometryTransformIdentity::read : format is not 1 ");
}

/***************  AstrometryTransformPolynomial **************************************/

//! Default transform : identity for all orders (>=1 )

AstrometryTransformPolynomial::AstrometryTransformPolynomial(std::size_t order) : _order(order) {
    _nterms = (order + 1) * (order + 2) / 2;

    // allocate and fill coefficients
    _coeffs.resize(2 * _nterms, 0.);
    // the default is supposed to be the identity, (for order>=1).
    if (_order >= 1) {
        getCoefficient(1, 0, 0) = 1;
        getCoefficient(0, 1, 1) = 1;
    }
}

//#ifdef TO_BE_FIXED
AstrometryTransformPolynomial::AstrometryTransformPolynomial(const AstrometryTransform *transform,
                                                             const Frame &frame, std::size_t order,
                                                             std::size_t nPoint) {
    StarMatchList sm;

    double step = std::sqrt(fabs(frame.getArea()) / double(nPoint));
    for (double x = frame.xMin + step / 2; x <= frame.xMax; x += step)
        for (double y = frame.yMin + step / 2; y <= frame.yMax; y += step) {
            auto pix = std::make_shared<BaseStar>(x, y, 0, 0);
            double xtr, ytr;
            transform->apply(x, y, xtr, ytr);
            auto tp = std::make_shared<BaseStar>(xtr, ytr, 0, 0);
            /* These are fake stars so no need to transform fake errors.
               all errors (and weights) will be equal : */
            sm.push_back(StarMatch(*pix, *tp, pix, tp));
        }
    AstrometryTransformPolynomial ret(order);
    ret.fit(sm);
    *this = ret;
}
//#endif

AstrometryTransformPolynomial::AstrometryTransformPolynomial(
        std::shared_ptr<afw::geom::TransformPoint2ToPoint2> transform, jointcal::Frame const &domain,
        std::size_t order, std::size_t nSteps) {
    jointcal::StarMatchList starMatchList;
    double xStart = domain.xMin;
    double yStart = domain.yMin;
    double xStep = domain.getWidth() / (nSteps + 1);
    double yStep = domain.getHeight() / (nSteps + 1);
    for (std::size_t i = 0; i < nSteps; ++i) {
        for (std::size_t j = 0; j < nSteps; ++j) {
            // TODO: once DM-4044 is done, we can remove the redundancy in `Point`/`Point2D` here
            jointcal::Point in(xStart + i * xStep, yStart + j * yStep);
            geom::Point2D inAfw(in.x, in.y);
            geom::Point2D outAfw = transform->applyForward(inAfw);
            jointcal::Point out(outAfw.getX(), outAfw.getY());
            starMatchList.emplace_back(in, out, nullptr, nullptr);
        }
    }
    AstrometryTransformPolynomial poly(order);
    poly.fit(starMatchList);
    *this = poly;
}

void AstrometryTransformPolynomial::computeMonomials(double xIn, double yIn, double *monomial) const {
    /* The ordering of monomials is implemented here.
       You may not change it without updating the "mapping" routines
      getCoefficient(std::size_t, std::size_t, std::size_t).
      I (P.A.) did not find a clever way to loop over monomials.
      Improvements welcome.
      This routine is used also by the fit to fill monomials.
      We could certainly be more elegant.
    */

    double xx = 1;
    for (std::size_t ix = 0; ix <= _order; ++ix) {
        double yy = 1;
        std::size_t k = ix * (ix + 1) / 2;
        for (std::size_t iy = 0; iy <= _order - ix; ++iy) {
            monomial[k] = xx * yy;
            yy *= yIn;
            k += ix + iy + 2;
        }
        xx *= xIn;
    }
}

void AstrometryTransformPolynomial::setOrder(std::size_t order) {
    _order = order;
    std::size_t old_nterms = _nterms;
    _nterms = (_order + 1) * (_order + 2) / 2;

    // temporarily save coefficients
    vector<double> old_coeffs = _coeffs;
    // reallocate enough size
    _coeffs.resize(2 * _nterms);
    // reassign to zero (this is necessary because ycoeffs
    // are after xcoeffs and so their meaning changes
    for (std::size_t k = 0; k < _nterms; ++k) _coeffs[k] = 0;
    // put back what we had before
    std::size_t kmax = min(old_nterms, _nterms);
    for (std::size_t k = 0; k < kmax; ++k) {
        _coeffs[k] = old_coeffs[k];                         // x terms
        _coeffs[k + _nterms] = old_coeffs[k + old_nterms];  // y terms
    }
}

/* this is reasonably fast, when optimized */
void AstrometryTransformPolynomial::apply(const double xIn, const double yIn, double &xOut,
                                          double &yOut) const {
    /*
      This routine computes the monomials only once for both
      polynomials.  This is why AstrometryTransformPolynomial does not use an auxilary
      class (such as PolyXY) to handle each polynomial.

      The code works even if &xIn == &xOut (or &yIn == &yOut)
      It uses Variable Length Allocation (VLA) rather than a vector<double>
      because allocating the later costs about 50 ns. All VLA uses are tagged.
    */
    double monomials[_nterms];  // this is VLA, which is (perhaps) not casher C++
    computeMonomials(xIn, yIn, monomials);

    xOut = 0;
    yOut = 0;
    const double *c = &_coeffs[0];
    const double *pm = &monomials[0];
    // the ordering of the coefficients and the monomials are identical.
    for (int k = _nterms; k--;) xOut += (*(pm++)) * (*(c++));
    pm = &monomials[0];
    for (int k = _nterms; k--;) yOut += (*(pm++)) * (*(c++));
}

void AstrometryTransformPolynomial::computeDerivative(Point const &where,
                                                      AstrometryTransformLinear &derivative,
                                                      const double step)
        const { /* routine checked against numerical derivatives from AstrometryTransform::Derivative */
    if (_order == 1) {
        derivative = AstrometryTransformLinear(*this);
        derivative.dx() = derivative.dy() = 0;
        return;
    }

    double dermx[2 * _nterms];  // VLA
    double *dermy = dermx + _nterms;
    double xin = where.x;
    double yin = where.y;

    double xx = 1;
    double xxm1 = 1;  // xx^(ix-1)
    for (std::size_t ix = 0; ix <= _order; ++ix) {
        std::size_t k = (ix) * (ix + 1) / 2;
        // iy = 0
        dermx[k] = ix * xxm1;
        dermy[k] = 0;
        k += ix + 2;
        double yym1 = 1;  // yy^(iy-1)
        for (std::size_t iy = 1; iy <= _order - ix; ++iy) {
            dermx[k] = ix * xxm1 * yym1 * yin;
            dermy[k] = iy * xx * yym1;
            yym1 *= yin;
            k += ix + iy + 2;
        }
        xx *= xin;
        if (ix >= 1) xxm1 *= xin;
    }

    derivative.dx() = 0;
    derivative.dy() = 0;

    const double *mx = &dermx[0];
    const double *my = &dermy[0];
    const double *c = &_coeffs[0];
    // dx'
    double a11 = 0, a12 = 0;
    for (int k = _nterms; k--;) {
        a11 += (*(mx++)) * (*c);
        a12 += (*(my++)) * (*(c++));
    }
    derivative.a11() = a11;
    derivative.a12() = a12;
    // dy'
    double a21 = 0, a22 = 0;
    mx = &dermx[0];
    my = &dermy[0];
    for (int k = _nterms; k--;) {
        a21 += (*(mx++)) * (*c);
        a22 += (*(my++)) * (*(c++));
    }
    derivative.a21() = a21;
    derivative.a22() = a22;
}

void AstrometryTransformPolynomial::transformPosAndErrors(FatPoint const &in, FatPoint &out) const {
    /*
       The results from this routine were compared to what comes out
       from apply and transformErrors. The Derivative routine was
       checked against numerical derivatives from
       AstrometryTransform::Derivative. (P.A dec 2009).

       This routine could be made much simpler by calling apply and
       Derivative (i.e. you just suppress it, and the fallback is the
       generic version in AstrometryTransform).  BTW, I checked that both routines
       provide the same result. This version is however faster
       (monomials get recycled).
    */
    double monomials[_nterms];  // VLA

    FatPoint res;  // to store the result, because nothing forbids &in == &out.

    double dermx[2 * _nterms];        // monomials for derivative w.r.t. x (VLA)
    double *dermy = dermx + _nterms;  // same for y
    double xin = in.x;
    double yin = in.y;

    double xx = 1;
    double xxm1 = 1;  // xx^(ix-1)
    for (std::size_t ix = 0; ix <= _order; ++ix) {
        std::size_t k = (ix) * (ix + 1) / 2;
        // iy = 0
        dermx[k] = ix * xxm1;
        dermy[k] = 0;
        monomials[k] = xx;
        k += ix + 2;
        double yy = yin;
        double yym1 = 1;  // yy^(iy-1)
        for (std::size_t iy = 1; iy <= _order - ix; ++iy) {
            monomials[k] = xx * yy;
            dermx[k] = ix * xxm1 * yy;
            dermy[k] = iy * xx * yym1;
            yym1 *= yin;
            yy *= yin;
            k += ix + iy + 2;
        }
        xx *= xin;
        if (ix >= 1) xxm1 *= xin;
    }

    // output position
    double xout = 0, yout = 0;
    const double *c = &_coeffs[0];
    const double *pm = &monomials[0];
    for (int k = _nterms; k--;) xout += (*(pm++)) * (*(c++));
    pm = &monomials[0];
    for (int k = _nterms; k--;) yout += (*(pm++)) * (*(c++));
    res.x = xout;
    res.y = yout;

    // derivatives
    c = &_coeffs[0];
    const double *mx = &dermx[0];
    const double *my = &dermy[0];
    double a11 = 0, a12 = 0;
    for (int k = _nterms; k--;) {
        a11 += (*(mx++)) * (*c);
        a12 += (*(my++)) * (*(c++));
    }

    double a21 = 0, a22 = 0;
    mx = &dermx[0];
    my = &dermy[0];
    for (int k = _nterms; k--;) {
        a21 += (*(mx++)) * (*c);
        a22 += (*(my++)) * (*(c++));
    }

    // output co-variance
    res.vx = a11 * (a11 * in.vx + 2 * a12 * in.vxy) + a12 * a12 * in.vy;
    res.vy = a21 * a21 * in.vx + a22 * a22 * in.vy + 2. * a21 * a22 * in.vxy;
    res.vxy = a21 * a11 * in.vx + a22 * a12 * in.vy + (a21 * a12 + a11 * a22) * in.vxy;
    out = res;
}

/* The coefficient ordering is defined both here *AND* in the
   AstrometryTransformPolynomial::apply, AstrometryTransformPolynomial::Derivative, ... routines
   Change all or none ! */

double AstrometryTransformPolynomial::getCoefficient(std::size_t degX, std::size_t degY,
                                                     std::size_t whichCoord) const {
    assert((degX + degY <= _order) && whichCoord < 2);
    /* this assertion above is enough to ensure that the index used just
       below is within bounds since the reserved length is
       2*_nterms=(order+1)*(order+2) */
    return _coeffs[(degX + degY) * (degX + degY + 1) / 2 + degY + whichCoord * _nterms];
}

double &AstrometryTransformPolynomial::getCoefficient(std::size_t degX, std::size_t degY,
                                                      std::size_t whichCoord) {
    assert((degX + degY <= _order) && whichCoord < 2);
    return _coeffs[(degX + degY) * (degX + degY + 1) / 2 + degY + whichCoord * _nterms];
}

double AstrometryTransformPolynomial::coeffOrZero(std::size_t degX, std::size_t degY,
                                                  std::size_t whichCoord) const {
    //  assert((degX+degY<=order) && whichCoord<2);
    assert(whichCoord < 2);
    if (degX + degY <= _order)
        return _coeffs[(degX + degY) * (degX + degY + 1) / 2 + degY + whichCoord * _nterms];
    return 0;
}

/* parameter serialization for "virtual" fits */
double AstrometryTransformPolynomial::paramRef(Eigen::Index const i) const {
    assert(i < 2 * Eigen::Index(_nterms));
    return _coeffs[i];
}

double &AstrometryTransformPolynomial::paramRef(Eigen::Index const i) {
    assert(i < 2 * Eigen::Index(_nterms));
    return _coeffs[i];
}

void AstrometryTransformPolynomial::paramDerivatives(Point const &where, double *dx, double *dy)
        const { /* first half : dxout/dpar, second half : dyout/dpar */
    computeMonomials(where.x, where.y, dx);
    for (std::size_t k = 0; k < _nterms; ++k) {
        dy[_nterms + k] = dx[k];
        dx[_nterms + k] = dy[k] = 0;
    }
}

/* utility for the print(ostream&) routine */
static string monomialString(std::size_t powX, std::size_t powY) {
    stringstream ss;
    if (powX + powY) ss << "*";
    if (powX > 0) ss << "x";
    if (powX > 1) ss << "^" << powX;
    if (powY > 0) ss << "y";
    if (powY > 1) ss << "^" << powY;
    return ss.str();
}

void AstrometryTransformPolynomial::print(ostream &stream) const {
    auto oldPrecision = stream.precision();
    stream.precision(12);
    stream << "AstrometryTransformPolynomial: order=" << getOrder() << std::endl;
    for (std::size_t ic = 0; ic < 2; ++ic) {
        if (ic == 0)
            stream << "newx = ";
        else
            stream << "newy = ";
        bool printed = false;  // whether we've printed any monomials
        for (std::size_t p = 0; p <= _order; ++p)
            for (std::size_t py = 0; py <= p; ++py) {
                double coefficient = getCoefficient(p - py, py, ic);
                if (coefficient != 0 || isnan(coefficient)) {
                    // Only print "interesting" coefficients.
                    if (printed && p + py != 0) stream << " + ";
                    printed = true;
                    stream << coefficient << monomialString(p - py, py);
                }
            }
        // if all components are zero, nothing is output by the loops above
        if (!printed) {
            stream << 0;
        }
        stream << endl;
    }
    if (_order > 0) stream << "Linear determinant = " << determinant();
    stream.precision(oldPrecision);
}

double AstrometryTransformPolynomial::determinant() const {
    if (_order >= 1)
        return getCoefficient(1, 0, 0) * getCoefficient(0, 1, 1) -
               getCoefficient(0, 1, 0) * getCoefficient(1, 0, 1);
    return 0;
}

//! Returns the transformation that maps the input frame along both axes to [-1,1]
AstrometryTransformLinear normalizeCoordinatesTransform(const Frame &frame) {
    Point center = frame.getCenter();
    return AstrometryTransformLinearScale(2. / frame.getWidth(), 2. / frame.getHeight()) *
           AstrometryTransformLinearShift(-center.x, -center.y);
}

/*utility for the AstrometryTransformPolynomial::fit() routine */
static AstrometryTransformLinear shiftAndNormalize(StarMatchList const &starMatchList) {
    double xav = 0;
    double x2 = 0;
    double yav = 0;
    double y2 = 0;
    double count = 0;
    for (auto it = starMatchList.begin(); it != starMatchList.end(); ++it) {
        const StarMatch &a_match = *it;
        Point const &point1 = a_match.point1;
        xav += point1.x;
        yav += point1.y;
        x2 += std::pow(point1.x, 2);
        y2 += std::pow(point1.y, 2);
        count++;
    }
    if (count == 0) return AstrometryTransformLinear();
    xav /= count;
    yav /= count;
    // 3.5 stands for sqrt(12).
    double xspan = 3.5 * std::sqrt(x2 / count - std::pow(xav, 2));
    double yspan = 3.5 * std::sqrt(y2 / count - std::pow(yav, 2));
    return AstrometryTransformLinearScale(2. / xspan, 2. / yspan) *
           AstrometryTransformLinearShift(-xav, -yav);
}

static double sq(double x) { return x * x; }

double AstrometryTransformPolynomial::computeFit(StarMatchList const &starMatchList,
                                                 AstrometryTransform const &shiftToCenter,
                                                 const bool useErrors) {
    Eigen::MatrixXd A(2 * _nterms, 2 * _nterms);
    A.setZero();
    Eigen::VectorXd B(2 * _nterms);
    B.setZero();
    double sumr2 = 0;
    double monomials[_nterms];
    for (auto it = starMatchList.begin(); it != starMatchList.end(); ++it) {
        const StarMatch &a_match = *it;
        Point tmp = shiftToCenter.apply(a_match.point1);
        FatPoint point1(tmp, a_match.point1.vx, a_match.point1.vy, a_match.point1.vxy);
        FatPoint const &point2 = a_match.point2;
        double wxx, wyy, wxy;
        FatPoint tr1;
        computeMonomials(point1.x, point1.y, monomials);
        if (useErrors) {
            transformPosAndErrors(point1, tr1);  // we might consider recycling the monomials
            double vxx = (tr1.vx + point2.vx);
            double vyy = (tr1.vy + point2.vy);
            double vxy = (tr1.vxy + point2.vxy);
            double det = vxx * vyy - vxy * vxy;
            wxx = vyy / det;
            wyy = vxx / det;
            wxy = -vxy / det;
        } else {
            wxx = wyy = 1;
            wxy = 0;
            apply(point1.x, point1.y, tr1.x, tr1.y);
        }
        double resx = point2.x - tr1.x;
        double resy = point2.y - tr1.y;
        sumr2 += wxx * sq(resx) + wyy * sq(resy) + 2 * wxy * resx * resy;

        double bxcoeff = wxx * resx + wxy * resy;
        double bycoeff = wyy * resy + wxy * resx;
        for (std::size_t j = 0; j < _nterms; ++j) {
            for (std::size_t i = j; i < _nterms; ++i) {
                A(i, j) += wxx * monomials[i] * monomials[j];
                A(i + _nterms, j + _nterms) += wyy * monomials[i] * monomials[j];
                A(j, i + _nterms) = A(i, j + _nterms) += wxy * monomials[i] * monomials[j];
            }
            B(j) += bxcoeff * monomials[j];
            B(j + _nterms) += bycoeff * monomials[j];
        }
    }  // end loop on points
    Eigen::LDLT<Eigen::MatrixXd, Eigen::Lower> factor(A);
    // should probably throw
    if (factor.info() != Eigen::Success) {
        LOGL_ERROR(_log, "AstrometryTransformPolynomial::fit could not factorize");
        return -1;
    }

    Eigen::VectorXd sol = factor.solve(B);
    for (std::size_t k = 0; k < 2 * _nterms; ++k) _coeffs[k] += sol(k);
    if (starMatchList.size() == _nterms) return 0;
    return (sumr2 - B.dot(sol));
}

double AstrometryTransformPolynomial::fit(StarMatchList const &starMatchList) {
    if (starMatchList.size() < _nterms) {
        LOGLS_FATAL(_log, "AstrometryTransformPolynomial::fit trying to fit a polynomial transform of order "
                                  << _order << " with only " << starMatchList.size() << " matches.");
        return -1;
    }

    AstrometryTransformPolynomial conditionner = shiftAndNormalize(starMatchList);

    computeFit(starMatchList, conditionner, false);               // get a rough solution
    computeFit(starMatchList, conditionner, true);                // weight with it
    double chi2 = computeFit(starMatchList, conditionner, true);  // once more

    (*this) = (*this) * conditionner;
    if (starMatchList.size() == _nterms) return 0;
    return chi2;
}

std::unique_ptr<AstrometryTransform> AstrometryTransformPolynomial::composeAndReduce(
        AstrometryTransformPolynomial const &right) const {
    if (getOrder() == 1 && right.getOrder() == 1)
        return std::make_unique<AstrometryTransformLinear>((*this) * (right));  // does the composition
    else
        return std::make_unique<AstrometryTransformPolynomial>((*this) * (right));  // does the composition
}

/*  PolyXY the class used to perform polynomial algebra (and in
    particular composition) at the coefficient level. This class
    handles a single polynomial, while a AstrometryTransformPolynomial is a couple of
    polynomials. This class does not have any routine to evaluate
    polynomials. Efficiency is not a concern since these routines are
    seldom used.  There is no need to expose this tool class to
    AstrometryTransform users.
*/

class PolyXY {
    std::size_t order;
    std::size_t nterms;
    vector<long_double> coeffs;

public:
    PolyXY(const int order) : order(order), nterms((order + 1) * (order + 2) / 2) {
        coeffs.reserve(nterms);
        coeffs.insert(coeffs.begin(), nterms, 0L);  // fill & initialize to 0.
    }

    std::size_t getOrder() const { return order; }

    PolyXY(AstrometryTransformPolynomial const &transform, std::size_t whichCoord)
            : order(transform.getOrder()), nterms((order + 1) * (order + 2) / 2), coeffs(nterms, 0L) {
        for (std::size_t px = 0; px <= order; ++px)
            for (std::size_t py = 0; py <= order - px; ++py) {
                getCoefficient(px, py) = transform.getCoefficient(px, py, whichCoord);
            }
    }

    long_double getCoefficient(std::size_t powX, std::size_t powY) const {
        assert(powX + powY <= order);
        return coeffs.at((powX + powY) * (powX + powY + 1) / 2 + powY);
    }

    long_double &getCoefficient(std::size_t powX, std::size_t powY) {
        assert(powX + powY <= order);
        return coeffs.at((powX + powY) * (powX + powY + 1) / 2 + powY);
    }
};

/* =====================  PolyXY Algebra routines ================== */

static void operator+=(PolyXY &left, const PolyXY &right) {
    std::size_t rdeg = right.getOrder();
    assert(left.getOrder() >= rdeg);
    for (std::size_t i = 0; i <= rdeg; ++i)
        for (std::size_t j = 0; j <= rdeg - i; ++j) left.getCoefficient(i, j) += right.getCoefficient(i, j);
}

/* multiplication by a scalar */
static PolyXY operator*(const long_double &a, const PolyXY &polyXY) {
    PolyXY result(polyXY);
    // no direct access to coefficients: do it the soft way
    std::size_t order = polyXY.getOrder();
    for (std::size_t i = 0; i <= order; ++i)
        for (std::size_t j = 0; j <= order - i; ++j) result.getCoefficient(i, j) *= a;
    return result;
}

/*! result(x,y) = p1(x,y)*p2(x,y) */
static PolyXY product(const PolyXY &p1, const PolyXY &p2) {
    std::size_t deg1 = p1.getOrder();
    std::size_t deg2 = p2.getOrder();
    PolyXY result(deg1 + deg2);
    for (std::size_t i1 = 0; i1 <= deg1; ++i1)
        for (std::size_t j1 = 0; j1 <= deg1 - i1; ++j1)
            for (std::size_t i2 = 0; i2 <= deg2; ++i2)
                for (std::size_t j2 = 0; j2 <= deg2 - i2; ++j2)
                    result.getCoefficient(i1 + i2, j1 + j2) +=
                            p1.getCoefficient(i1, j1) * p2.getCoefficient(i2, j2);
    return result;
}

/* powers[k](x,y) = polyXY(x,y)**k, 0 <= k <= maxP */
static void computePowers(const PolyXY &polyXY, std::size_t maxP, vector<PolyXY> &powers) {
    powers.reserve(maxP + 1);
    powers.push_back(PolyXY(0));
    powers[0].getCoefficient(0, 0) = 1L;
    for (std::size_t k = 1; k <= maxP; ++k) powers.push_back(product(powers[k - 1], polyXY));
}

/*! result(x,y) = polyXY(polyX(x,y),polyY(x,y)) */
static PolyXY composition(const PolyXY &polyXY, const PolyXY &polyX, const PolyXY &polyY) {
    std::size_t pdeg = polyXY.getOrder();
    PolyXY result(pdeg * max(polyX.getOrder(), polyY.getOrder()));
    vector<PolyXY> pXPowers;
    vector<PolyXY> pYPowers;
    computePowers(polyX, pdeg, pXPowers);
    computePowers(polyY, pdeg, pYPowers);
    for (std::size_t px = 0; px <= pdeg; ++px)
        for (std::size_t py = 0; py <= pdeg - px; ++py)
            result += polyXY.getCoefficient(px, py) * product(pXPowers.at(px), pYPowers.at(py));
    return result;
}

/* ===================== end of  PolyXY Algebra routines ============= */

/* reducing polynomial composition is the reason for PolyXY stuff : */

AstrometryTransformPolynomial AstrometryTransformPolynomial::operator*(
        AstrometryTransformPolynomial const &right) const {
    // split each transform into 2d polynomials
    PolyXY plx(*this, 0);
    PolyXY ply(*this, 1);
    PolyXY prx(right, 0);
    PolyXY pry(right, 1);

    // compute the compositions
    PolyXY rx(composition(plx, prx, pry));
    PolyXY ry(composition(ply, prx, pry));

    // copy the results the hard way.
    AstrometryTransformPolynomial result(_order * right._order);
    for (std::size_t px = 0; px <= result._order; ++px)
        for (std::size_t py = 0; py <= result._order - px; ++py) {
            result.getCoefficient(px, py, 0) = rx.getCoefficient(px, py);
            result.getCoefficient(px, py, 1) = ry.getCoefficient(px, py);
        }
    return result;
}

AstrometryTransformPolynomial AstrometryTransformPolynomial::operator+(
        AstrometryTransformPolynomial const &right) const {
    if (_order >= right._order) {
        AstrometryTransformPolynomial res(*this);
        for (std::size_t i = 0; i <= right._order; ++i)
            for (std::size_t j = 0; j <= right._order - i; ++j) {
                res.getCoefficient(i, j, 0) += right.getCoefficient(i, j, 0);
                res.getCoefficient(i, j, 1) += right.getCoefficient(i, j, 1);
            }
        return res;
    } else
        return (right + (*this));
}

AstrometryTransformPolynomial AstrometryTransformPolynomial::operator-(
        AstrometryTransformPolynomial const &right) const {
    AstrometryTransformPolynomial res(std::max(_order, right._order));
    for (std::size_t i = 0; i <= res._order; ++i)
        for (std::size_t j = 0; j <= res._order - i; ++j) {
            res.getCoefficient(i, j, 0) = coeffOrZero(i, j, 0) - right.coeffOrZero(i, j, 0);
            res.getCoefficient(i, j, 1) = coeffOrZero(i, j, 1) - right.coeffOrZero(i, j, 1);
        }
    return res;
}

std::shared_ptr<ast::Mapping> AstrometryTransformPolynomial::toAstMap(jointcal::Frame const &domain) const {
    auto inverse = inversePolyTransform(*this, domain, 1e-7, _order + 2, 100);
    return std::make_shared<ast::PolyMap>(toAstPolyMapCoefficients(), inverse->toAstPolyMapCoefficients());
}

void AstrometryTransformPolynomial::write(ostream &s) const {
    s << " AstrometryTransformPolynomial 1" << endl;
    s << "order " << _order << endl;
    int oldprec = s.precision();
    s << setprecision(12);
    for (std::size_t k = 0; k < 2 * _nterms; ++k) s << _coeffs[k] << ' ';
    s << endl;
    s << setprecision(oldprec);
}

void AstrometryTransformPolynomial::read(istream &s) {
    int format;
    s >> format;
    if (format != 1)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          " AstrometryTransformPolynomial::read : format is not 1 ");

    string order;
    s >> order >> _order;
    if (order != "order")
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          " AstrometryTransformPolynomial::read : expecting \"order\" and found " + order);
    setOrder(_order);
    for (std::size_t k = 0; k < 2 * _nterms; ++k) s >> _coeffs[k];
}

ndarray::Array<double, 2, 2> AstrometryTransformPolynomial::toAstPolyMapCoefficients() const {
    std::size_t nCoeffs = _coeffs.size();
    ndarray::Array<double, 2, 2> result = ndarray::allocate(ndarray::makeVector(nCoeffs, std::size_t(4)));

    ndarray::Size k = 0;
    for (std::size_t iCoord = 0; iCoord < 2; ++iCoord) {
        for (std::size_t p = 0; p <= _order; ++p) {
            for (std::size_t py = 0; py <= p; ++py, ++k) {
                result[k][0] = getCoefficient(p - py, py, iCoord);
                result[k][1] = iCoord + 1;
                result[k][2] = p - py;
                result[k][3] = py;
            }
        }
    }

    return result;
}

std::shared_ptr<AstrometryTransformPolynomial> inversePolyTransform(AstrometryTransform const &forward,
                                                                    Frame const &domain,
                                                                    double const precision,
                                                                    std::size_t maxOrder,
                                                                    std::size_t nSteps) {
    StarMatchList sm;
    double xStart = domain.xMin;
    double yStart = domain.yMin;
    double xStep = domain.getWidth() / (nSteps - 1);
    double yStep = domain.getHeight() / (nSteps - 1);
    for (std::size_t i = 0; i < nSteps; ++i) {
        for (std::size_t j = 0; j < nSteps; ++j) {
            Point in(xStart + i * xStep, yStart + j * yStep);
            Point out(forward.apply(in));
            sm.push_back(StarMatch(out, in, nullptr, nullptr));
        }
    }
    std::size_t npairs = sm.size();
    std::size_t order;
    std::shared_ptr<AstrometryTransformPolynomial> poly;
    std::shared_ptr<AstrometryTransformPolynomial> oldPoly;
    double chi2 = 0;
    double oldChi2 = std::numeric_limits<double>::infinity();
    for (order = 1; order <= maxOrder; ++order) {
        poly.reset(new AstrometryTransformPolynomial(order));
        auto success = poly->fit(sm);
        if (success == -1) {
            std::stringstream errMsg;
            errMsg << "Cannot fit a polynomial of order " << order << " with " << nSteps << "^2 points";
            throw pexExcept::RuntimeError(errMsg.str());
        }
        // compute the chi2 ignoring errors:
        chi2 = 0;
        for (auto const &i : sm) chi2 += i.point2.computeDist2(poly->apply((i.point1)));
        LOGLS_TRACE(_log, "inversePoly order " << order << ": " << chi2 << " / " << npairs << " = "
                                               << chi2 / npairs << " < " << precision * precision);

        if (chi2 / npairs < precision * precision) break;

        // If this triggers, we know we did not reach the required precision.
        if (chi2 > oldChi2) {
            LOGLS_WARN(_log, "inversePolyTransform: chi2 increases ("
                                     << chi2 << " > " << oldChi2 << "); ending fit with order: " << order);
            LOGLS_WARN(_log, "inversePolyTransform: requested precision not reached: "
                                     << chi2 << " / " << npairs << " = " << chi2 / npairs << " < "
                                     << precision * precision);
            poly = std::move(oldPoly);
            order--;
            break;
        } else {
            oldChi2 = chi2;
            // Clone it so we don't lose it in the next iteration.
            oldPoly = dynamic_pointer_cast<AstrometryTransformPolynomial>(
                    std::shared_ptr<AstrometryTransform>(poly->clone()));
        }
    }
    if (order > maxOrder)
        LOGLS_WARN(_log, "inversePolyTransform: Reached max order without reaching requested precision: "
                                 << chi2 << " / " << npairs << " = " << chi2 / npairs << " < "
                                 << precision * precision);
    return poly;
}

/**************** AstrometryTransformLinear ***************************************/
/* AstrometryTransformLinear is a specialized constructor of AstrometryTransformPolynomial
   May be it could just disappear ??
*/

AstrometryTransformLinear::AstrometryTransformLinear(const double Dx, const double Dy, const double A11,
                                                     const double A12, const double A21, const double A22)
        : AstrometryTransformPolynomial(1) {
    dx() = Dx;
    a11() = A11;
    a12() = A12;
    dy() = Dy;
    a21() = A21;
    a22() = A22;
}

AstrometryTransformLinear::AstrometryTransformLinear(AstrometryTransformPolynomial const &transform)
        : AstrometryTransformPolynomial(1) {
    if (transform.getOrder() != 1)
        throw pexExcept::InvalidParameterError(
                "Trying to build a AstrometryTransformLinear from a higher order transform. Aborting. ");
    (AstrometryTransformPolynomial &)(*this) = transform;
}

AstrometryTransformLinear AstrometryTransformLinear::operator*(AstrometryTransformLinear const &right) const {
    // There is a general routine in AstrometryTransformPolynomial that would do the job:
    //  return AstrometryTransformLinear(AstrometryTransformPolynomial::operator*(right));
    // however, we are using this composition of linear stuff heavily in
    // AstrometryTransform::linearApproximation, itself used in inverseTransform::apply.
    // So, for sake of efficiency, and since it is easy, we take a shortcut:
    AstrometryTransformLinear result;
    apply(right.Dx(), right.Dy(), result.dx(), result.dy());
    result.a11() = A11() * right.A11() + A12() * right.A21();
    result.a12() = A11() * right.A12() + A12() * right.A22();
    result.a21() = A21() * right.A11() + A22() * right.A21();
    result.a22() = A21() * right.A12() + A22() * right.A22();
    return result;
}

void AstrometryTransformLinear::computeDerivative(Point const &, AstrometryTransformLinear &derivative,
                                                  const double) const {
    derivative = *this;
    derivative.getCoefficient(0, 0, 0) = 0;
    derivative.getCoefficient(0, 0, 1) = 0;
}

AstrometryTransformLinear AstrometryTransformLinear::linearApproximation(Point const &, const double) const {
    return *this;
}

AstrometryTransformLinear AstrometryTransformLinear::inverted() const {
    //
    //   (T1,M1) * (T2,M2) = 1 i.e (0,1) implies
    //   T1 = -M1 * T2
    //   M1 = M2^(-1)
    //

    double a11 = A11();
    double a12 = A12();
    double a21 = A21();
    double a22 = A22();
    double d = (a11 * a22 - a12 * a21);
    if (d == 0) {
        LOGL_FATAL(_log,
                   "AstrometryTransformLinear::invert singular transformation: transform contents will be "
                   "dumped to stderr.");
        print(cerr);
    }

    AstrometryTransformLinear result(0, 0, a22 / d, -a12 / d, -a21 / d, a11 / d);
    double rdx, rdy;
    result.apply(Dx(), Dy(), rdx, rdy);
    result.dx() = -rdx;
    result.dy() = -rdy;
    return result;
}

std::unique_ptr<AstrometryTransform> AstrometryTransformLinear::inverseTransform(const double,
                                                                                 const Frame &) const {
    return std::unique_ptr<AstrometryTransform>(new AstrometryTransformLinear(inverted()));
}

void AstrometryTransformLinear::print(ostream &stream) const {
    auto oldPrecision = stream.precision();
    stream.precision(12);
    stream << "Linear AstrometryTransform:" << std::endl;
    stream << A11() << " " << A12() << " + " << Dx() << std::endl;
    stream << A21() << " " << A22() << " + " << Dy() << std::endl;
    stream.precision(oldPrecision);
    stream << "determinant = " << determinant();
    stream.precision(oldPrecision);
}

double AstrometryTransformLinearRot::fit(StarMatchList const &) {
    throw pexExcept::NotFoundError("AstrometryTransformLinearRot::fit not implemented! aborting");
}

double AstrometryTransformLinearShift::fit(StarMatchList const &starMatchList) {
    std::size_t npairs = starMatchList.size();
    if (npairs < 3) {
        LOGLS_FATAL(_log, "AstrometryTransformLinearShift::fit trying to fit a linear transform with only "
                                  << npairs << " matches.");
        return -1;
    }

    double sumr2 = 0; /* used to compute chi2 without relooping */
    /* loop on pairs  and fill */
    Eigen::VectorXd B(2);
    B.setZero();
    Eigen::MatrixXd A(2, 2);
    A.setZero();

    for (auto const &it : starMatchList) {
        FatPoint const &point1 = it.point1;
        FatPoint const &point2 = it.point2;
        double deltax = point2.x - point1.x;
        double deltay = point2.y - point1.y;
        double vxx = point1.vx + point2.vx;
        double vyy = point1.vy + point2.vy;
        double vxy = point1.vxy + point2.vxy;
        double det = vxx * vyy - vxy * vxy;
        double wxx = vyy / det;
        double wyy = vxx / det;
        double wxy = -vxy / det;
        B(0) += deltax * wxx + wxy * deltay;
        B(1) += deltay * wyy + wxy * deltax;
        A(0, 0) += wxx;
        A(1, 1) += wyy;
        A(0, 1) += wxy;
        sumr2 += deltax * deltax * wxx + deltay * deltay * wyy + 2. * wxy * deltax * deltay;
    }
    double det = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
    if (det <= 0) return -1;
    double tmp = A(0, 0);
    A(0, 0) = A(1, 1) / det;
    A(1, 1) = tmp / det;
    A(0, 1) = A(1, 0) = -A(0, 1) / det;
    Eigen::VectorXd sol = A * B;
    (*this) = AstrometryTransformLinearShift(sol(0), sol(1));
    return (sumr2 - sol.dot(B));  // chi2 compact form
}

AstrometryTransformLinearRot::AstrometryTransformLinearRot(const double angleRad, const Point *center,
                                                           const double scaleFactor) {
    double c = scaleFactor * std::cos(angleRad);
    double s = scaleFactor * std::sin(angleRad);
    a11() = a22() = c;
    a21() = s;
    a12() = -s;

    // we want that the center does not move : transform+M*C = C ==> transform = C - M*C
    Point a_point(0., 0.);
    if (center) a_point = *center;

    dx() = dy() = 0;
    AstrometryTransformPolynomial::apply(a_point.x, a_point.y, dx(), dy());  // compute M*C
    dx() = a_point.x - Dx();
    dy() = a_point.y - dy();
}

static double deg2rad(double degree) { return degree * M_PI / 180.; }

static double rad2deg(double rad) { return rad * 180. / M_PI; }

/*************  WCS transform ******************/
/************** LinPixelToTan *******************/

/* Implementation note : it seemed wise to incorporate
   the radians to degreess convertion into the linPixelToTan
   part (right in the constructor), and to do the
   opposite operation in the LinPart routine.
   When I was coding the fit, I realized that it was a
   bad idea. then I realized that the fitting routine
   itself was probably useless. Finaly, for a persistor,
   it seems bizarre that the stored data is different
   from what convention (angles in degrees for WCS's)
   would expect.
   So, no more "automatic" degrees to radians and
   radians to degrees conversion. They are explicitely
   done in apply (for TanPixelToRaDec and TanRaDecToPixel).
   This is a minor concern though....
*/
BaseTanWcs::BaseTanWcs(AstrometryTransformLinear const &pixToTan, Point const &tangentPoint,
                       const AstrometryTransformPolynomial *corrections) {
    // the angles returned by linPixelToTan should be in degrees.
    linPixelToTan = pixToTan;
    ra0 = deg2rad(tangentPoint.x);
    dec0 = deg2rad(tangentPoint.y);
    cos0 = std::cos(dec0);
    sin0 = std::sin(dec0);
    corr = nullptr;
    if (corrections) corr.reset(new AstrometryTransformPolynomial(*corrections));
}

/* with some sort of smart pointer ro handle "corr", we could remove the
   copy constructor, the operator = and the destructor */

// ": AstrometryTransform" suppresses a warning
BaseTanWcs::BaseTanWcs(const BaseTanWcs &original) : AstrometryTransform() {
    corr = nullptr;
    *this = original;
}

void BaseTanWcs::operator=(const BaseTanWcs &original) {
    linPixelToTan = original.linPixelToTan;
    ra0 = original.ra0;
    dec0 = original.dec0;
    cos0 = std::cos(dec0);
    sin0 = std::sin(dec0);
    corr = nullptr;
    if (original.corr) corr.reset(new AstrometryTransformPolynomial(*original.corr));
}

void BaseTanWcs::apply(const double xIn, const double yIn, double &xOut, double &yOut) const {
    double l, m;                        // radians in the tangent plane
    pixToTangentPlane(xIn, yIn, l, m);  // l, m in degrees.
    l = deg2rad(l);
    m = deg2rad(m);  // now in radians
                     // Code inspired from worldpos.c in wcssubs (ancestor of the wcslib)
                     /* At variance with wcslib, it collapses the projection to a plane
                        and expression of sidereal cooordinates into a single set of
                        operations. */
    double dect = cos0 - m * sin0;
    if (dect == 0) {
        LOGL_WARN(_log, "No sidereal coordinates at pole!");
        xOut = 0;
        yOut = 0;
        return;
    }
    double rat = ra0 + atan2(l, dect);
    dect = atan(std::cos(rat - ra0) * (m * cos0 + sin0) / dect);
    if (rat - ra0 > M_PI) rat -= (2. * M_PI);
    if (rat - ra0 < -M_PI) rat += (2. * M_PI);
    if (rat < 0.0) rat += (2. * M_PI);
    // convert to degree
    xOut = rad2deg(rat);
    yOut = rad2deg(dect);
}

Point BaseTanWcs::getTangentPoint() const { return Point(rad2deg(ra0), rad2deg(dec0)); }

AstrometryTransformLinear BaseTanWcs::getLinPart() const { return linPixelToTan; }

void BaseTanWcs::setCorrections(std::unique_ptr<AstrometryTransformPolynomial> corrections) {
    corr = std::move(corrections);
}

Point BaseTanWcs::getCrPix() const {
    /* CRPIX's are defined by:
                      ( CD1_1  CD1_2 )   (x - crpix1)
       transformed =  (              ) * (          )
                      ( CD2_1  CD2_2 )   (y - crpix2)

       so that CrPix is the point which transforms to (0,0)
    */
    const AstrometryTransformLinear inverse = linPixelToTan.inverted();
    return Point(inverse.Dx(), inverse.Dy());
}

BaseTanWcs::~BaseTanWcs() {}

AstrometryTransformSkyWcs::AstrometryTransformSkyWcs(std::shared_ptr<afw::geom::SkyWcs> skyWcs)
        : _skyWcs(skyWcs) {}

void AstrometryTransformSkyWcs::apply(const double xIn, const double yIn, double &xOut, double &yOut) const {
    auto const outCoord = _skyWcs->pixelToSky(geom::Point2D(xIn, yIn));
    xOut = outCoord[0].asDegrees();
    yOut = outCoord[1].asDegrees();
}

void AstrometryTransformSkyWcs::print(std::ostream &out) const {
    out << "AstrometryTransformSkyWcs(" << *_skyWcs << ")";
}

double AstrometryTransformSkyWcs::fit(const StarMatchList &starMatchList) {
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not implemented");
}

std::unique_ptr<AstrometryTransform> AstrometryTransformSkyWcs::clone() const {
    return std::unique_ptr<AstrometryTransformSkyWcs>(new AstrometryTransformSkyWcs(getSkyWcs()));
}

/*************************** TanPixelToRaDec ***************/

TanPixelToRaDec::TanPixelToRaDec(AstrometryTransformLinear const &pixToTan, Point const &tangentPoint,
                                 const AstrometryTransformPolynomial *corrections)
        : BaseTanWcs(pixToTan, tangentPoint, corrections) {}

// ": AstrometryTransform" suppresses a warning
TanPixelToRaDec::TanPixelToRaDec() : BaseTanWcs(AstrometryTransformLinear(), Point(0, 0), nullptr) {}

std::unique_ptr<AstrometryTransform> TanPixelToRaDec::composeAndReduce(
        AstrometryTransformLinear const &right) const {
    if (right.getOrder() == 1) {
        return std::make_unique<TanPixelToRaDec>((*this) * (right));
    } else {
        return std::unique_ptr<AstrometryTransform>(nullptr);
    }
}

TanPixelToRaDec TanPixelToRaDec::operator*(AstrometryTransformLinear const &right) const {
    TanPixelToRaDec result(*this);
    result.linPixelToTan = result.linPixelToTan * right;
    return result;
}

TanRaDecToPixel TanPixelToRaDec::inverted() const {
    if (corr != nullptr) {
        LOGL_WARN(_log, "You are inverting a TanPixelToRaDec with corrections.");
        LOGL_WARN(_log, "The inverse you get ignores the corrections!");
    }
    return TanRaDecToPixel(getLinPart().inverted(), getTangentPoint());
}

std::unique_ptr<AstrometryTransform> TanPixelToRaDec::roughInverse(const Frame &) const {
    return std::unique_ptr<AstrometryTransform>(
            new TanRaDecToPixel(getLinPart().inverted(), getTangentPoint()));
}

std::unique_ptr<AstrometryTransform> TanPixelToRaDec::inverseTransform(const double precision,
                                                                       const Frame &region) const {
    if (!corr)
        return std::unique_ptr<AstrometryTransform>(
                new TanRaDecToPixel(getLinPart().inverted(), getTangentPoint()));
    else
        return std::unique_ptr<AstrometryTransform>(new AstrometryTransformInverse(this, precision, region));
}

AstrometryTransformPolynomial TanPixelToRaDec::getPixelToTangentPlane() const {
    if (corr)
        return (*corr) * linPixelToTan;
    else
        return linPixelToTan;
}

void TanPixelToRaDec::pixToTangentPlane(double xPixel, double yPixel, double &xTangentPlane,
                                        double &yTangentPlane) const {
    // xTangentPlane, yTangentPlane in degrees.
    linPixelToTan.apply(xPixel, yPixel, xTangentPlane, yTangentPlane);
    if (corr) {
        double xtmp = xTangentPlane;
        double ytmp = yTangentPlane;
        corr->apply(xtmp, ytmp, xTangentPlane, yTangentPlane);  // still in degrees.
    }
}

std::unique_ptr<AstrometryTransform> TanPixelToRaDec::clone() const {
    return std::unique_ptr<AstrometryTransform>(
            new TanPixelToRaDec(getLinPart(), getTangentPoint(), corr.get()));
}

void TanPixelToRaDec::print(ostream &stream) const {
    stream << " TanPixelToRaDec, lin part :" << endl << linPixelToTan;
    Point tp = getTangentPoint();
    stream << " tangent point " << tp.x << ' ' << tp.y << endl;
    Point crpix = getCrPix();
    stream << " crpix " << crpix.x << ' ' << crpix.y << endl;
    if (corr) stream << "PV correction: " << endl << *corr;
}

double TanPixelToRaDec::fit(StarMatchList const &) {
    /* OK we could implement this routine, but it is
       probably useless since to do the match, we have to
       project from sky to tangent plane. When a match is
       found, it is easier to carry out the fit in the
       tangent plane, rather than going back to the celestial
       sphere (and reproject to fit...). Anyway if this
       message shows up, we'll think about it.
    */
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                      "TanPixelToRaDec::fit is NOT implemented (although it is doable)) ");
    return -1;
}

/*************************** TanSipPixelToRaDec ***************/

TanSipPixelToRaDec::TanSipPixelToRaDec(AstrometryTransformLinear const &pixToTan, Point const &tangentPoint,
                                       const AstrometryTransformPolynomial *corrections)
        : BaseTanWcs(pixToTan, tangentPoint, corrections) {}

// ": AstrometryTransform" suppresses a warning
TanSipPixelToRaDec::TanSipPixelToRaDec() : BaseTanWcs(AstrometryTransformLinear(), Point(0, 0), nullptr) {}

/* Would require some checks before cooking up something more efficient
   than just a linear approximation */
#if 0
std::unique_ptr<AstrometryTransform> TanPixelToRaDec::roughInverse(const Frame &region) const
{
  if (&region) {}
  return std::unique_ptr<AstrometryTransform>(new TanRaDecToPixel(getLinPart().inverted(),getTangentPoint()));
}
#endif

std::unique_ptr<AstrometryTransform> TanSipPixelToRaDec::inverseTransform(const double precision,
                                                                          const Frame &region)
        const { /* We have not implemented (yet) the reverse corrections available in SIP */
    return std::unique_ptr<AstrometryTransform>(new AstrometryTransformInverse(this, precision, region));
}

AstrometryTransformPolynomial TanSipPixelToRaDec::getPixelToTangentPlane() const {
    if (corr)
        return AstrometryTransformPolynomial(linPixelToTan) * (*corr);
    else
        return linPixelToTan;
}

void TanSipPixelToRaDec::pixToTangentPlane(double xPixel, double yPixel, double &xTangentPlane,
                                           double &yTangentPlane) const {
    // xTangentPlane, yTangentPlane returned in degrees
    if (corr) {
        double xtmp, ytmp;
        corr->apply(xPixel, yPixel, xtmp, ytmp);
        linPixelToTan.apply(xtmp, ytmp, xTangentPlane, yTangentPlane);
    } else
        linPixelToTan.apply(xPixel, yPixel, xTangentPlane, yTangentPlane);
}

std::unique_ptr<AstrometryTransform> TanSipPixelToRaDec::clone() const {
    return std::unique_ptr<AstrometryTransform>(
            new TanSipPixelToRaDec(getLinPart(), getTangentPoint(), corr.get()));
}

void TanSipPixelToRaDec::print(ostream &stream) const {
    stream << " TanSipPixelToRaDec, lin part :" << endl << linPixelToTan;
    Point tp = getTangentPoint();
    stream << " tangent point " << tp.x << ' ' << tp.y << endl;
    Point crpix = getCrPix();
    stream << " crpix " << crpix.x << ' ' << crpix.y << endl;
    if (corr) stream << "PV correction: " << endl << *corr;
}

double TanSipPixelToRaDec::fit(StarMatchList const &) {
    /* OK we could implement this routine, but it is
       probably useless since to do the match, we have to
       project from sky to tangent plane. When a match is
       found, it is easier to carry out the fit in the
       tangent plane, rather than going back to the celestial
       sphere (and reproject to fit...). Anyway if this
       message shows up, we'll think about it.
    */
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                      "TanSipPixelToRaDec::fit is NOT implemented (although it is doable)) ");
    return -1;
}

/***************  reverse transform of TanPixelToRaDec: TanRaDecToPixel ********/

TanRaDecToPixel::TanRaDecToPixel(AstrometryTransformLinear const &tan2Pix, Point const &tangentPoint)
        : linTan2Pix(tan2Pix) {
    setTangentPoint(tangentPoint);
}

void TanRaDecToPixel::setTangentPoint(Point const &tangentPoint) {
    /* the radian to degrees conversion after projection
        is handled in apply */
    ra0 = deg2rad(tangentPoint.x);
    dec0 = deg2rad(tangentPoint.y);
    cos0 = std::cos(dec0);
    sin0 = std::sin(dec0);
}

TanRaDecToPixel::TanRaDecToPixel() : linTan2Pix() {
    ra0 = dec0 = 0;
    cos0 = 1;
    sin0 = 0;
}

Point TanRaDecToPixel::getTangentPoint() const { return Point(rad2deg(ra0), rad2deg(dec0)); }

AstrometryTransformLinear TanRaDecToPixel::getLinPart() const { return linTan2Pix; }

// Use analytic derivatives, computed at the same time as the transform itself
void TanRaDecToPixel::transformPosAndErrors(FatPoint const &in, FatPoint &out) const {
    /* this routine is very similar to apply, but also propagates errors.
       The deg2rad and rad2deg are ignored for errors because they act as
       2 global scalings that cancel each other.
       Derivatives were computed using maple:

       l1 := sin(a - a0)*cos(d);
       m1 := sin(d)*sin(d0)+cos(d)*cos(d0)*cos(a-a0);
       l2 := sin(d)*cos(d0)-cos(d)*sin(d0)*cos(a-a0);
       simplify(diff(l1/m1,a));
       simplify(diff(l1/m1,d));
       simplify(diff(l2/m1,a));
       simplify(diff(l2/m1,d));

       Checked against AstrometryTransform::transformPosAndErrors (dec 09)
    */
    double ra = deg2rad(in.x);
    double dec = deg2rad(in.y);
    if (ra - ra0 > M_PI) ra -= (2. * M_PI);
    if (ra - ra0 < -M_PI) ra += (2. * M_PI);
    // Code inspired from worldpos.c in wcssubs (ancestor of the wcslib)
    // The same code is copied in ::apply()

    double coss = std::cos(dec);
    double sins = std::sin(dec);
    double sinda = std::sin(ra - ra0);
    double cosda = std::cos(ra - ra0);
    double l = sinda * coss;
    double m = sins * sin0 + coss * cos0 * cosda;
    l = l / m;
    m = (sins * cos0 - coss * sin0 * cosda) / m;

    // derivatives
    double deno =
            sq(sin0) - sq(coss) + sq(coss * cos0) * (1 + sq(cosda)) + 2 * sins * sin0 * coss * cos0 * cosda;
    double a11 = coss * (cosda * sins * sin0 + coss * cos0) / deno;
    double a12 = -sinda * sin0 / deno;
    double a21 = coss * sinda * sins / deno;
    double a22 = cosda / deno;

    FatPoint tmp;
    tmp.vx = a11 * (a11 * in.vx + 2 * a12 * in.vxy) + a12 * a12 * in.vy;
    tmp.vy = a21 * a21 * in.vx + a22 * a22 * in.vy + 2. * a21 * a22 * in.vxy;
    tmp.vxy = a21 * a11 * in.vx + a22 * a12 * in.vy + (a21 * a12 + a11 * a22) * in.vxy;

    // l and m are now coordinates in the tangent plane, in radians.
    tmp.x = rad2deg(l);
    tmp.y = rad2deg(m);

    linTan2Pix.transformPosAndErrors(tmp, out);
}

void TanRaDecToPixel::apply(const double xIn, const double yIn, double &xOut, double &yOut) const {
    double ra = deg2rad(xIn);
    double dec = deg2rad(yIn);
    if (ra - ra0 > M_PI) ra -= (2. * M_PI);
    if (ra - ra0 < -M_PI) ra += (2. * M_PI);
    // Code inspired from worldpos.c in wcssubs (ancestor of the wcslib)
    // The same code is copied in ::transformPosAndErrors()
    double coss = std::cos(dec);
    double sins = std::sin(dec);
    double l = std::sin(ra - ra0) * coss;
    double m = sins * sin0 + coss * cos0 * std::cos(ra - ra0);
    l = l / m;
    m = (sins * cos0 - coss * sin0 * std::cos(ra - ra0)) / m;
    // l and m are now coordinates in the tangent plane, in radians.
    l = rad2deg(l);
    m = rad2deg(m);
    linTan2Pix.apply(l, m, xOut, yOut);
}

TanPixelToRaDec TanRaDecToPixel::inverted() const {
    return TanPixelToRaDec(getLinPart().inverted(), getTangentPoint());
}

void TanRaDecToPixel::print(ostream &stream) const {
    Point tp = getTangentPoint();
    stream << "tan2pix " << linTan2Pix << std::endl;
    stream << "tangent point " << tp.x << ' ' << tp.y << endl;
}

std::unique_ptr<AstrometryTransform> TanRaDecToPixel::roughInverse(const Frame &) const {
    return std::unique_ptr<AstrometryTransform>(
            new TanPixelToRaDec(getLinPart().inverted(), getTangentPoint()));
}

std::unique_ptr<AstrometryTransform> TanRaDecToPixel::inverseTransform(const double, const Frame &) const {
    return std::unique_ptr<AstrometryTransform>(
            new TanPixelToRaDec(getLinPart().inverted(), getTangentPoint()));
}

std::unique_ptr<AstrometryTransform> TanRaDecToPixel::clone() const {
    return std::unique_ptr<AstrometryTransform>(new TanRaDecToPixel(*this));
}

double TanRaDecToPixel::fit(StarMatchList const &) {
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                      "TanRaDecToPixel::fit is NOT implemented (although it is doable)) ");
    return -1;
}

/*************  a "run-time" transform, that does not require to modify this file */

UserTransform::UserTransform(AstrometryTransformFun &userFun, const void *userData)
        : _userFun(userFun), _userData(userData) {}

void UserTransform::apply(const double xIn, const double yIn, double &xOut, double &yOut) const {
    _userFun(xIn, yIn, xOut, yOut, _userData);
}

void UserTransform::print(ostream &stream) const {
    stream << "UserTransform with user function @ " << _userFun << "and userData@ " << _userData << endl;
}

double UserTransform::fit(StarMatchList const &) {
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                      "UserTransform::fit is NOT implemented (and will never be)) ");
    return -1;
}

std::unique_ptr<AstrometryTransform> UserTransform::clone() const {
    return std::unique_ptr<AstrometryTransform>(new UserTransform(*this));
}

/*************************************************************/

std::unique_ptr<AstrometryTransform> astrometryTransformRead(const std::string &fileName) {
    ifstream s(fileName.c_str());
    if (!s)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          " astrometryTransformRead : cannot open " + fileName);
    try {
        std::unique_ptr<AstrometryTransform> res(astrometryTransformRead(s));
        s.close();
        return res;
    } catch (pex::exceptions::InvalidParameterError &e) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          std::string(e.what()) + " in file " + fileName);
    }
}

std::unique_ptr<AstrometryTransform> astrometryTransformRead(istream &s) {
    std::string type;
    s >> type;
    if (s.fail())
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "astrometryTransformRead : could not find a AstrometryTransformtype");
    if (type == "AstrometryTransformIdentity") {
        std::unique_ptr<AstrometryTransformIdentity> res(new AstrometryTransformIdentity());
        res->read(s);
        return std::move(res);
    } else if (type == "AstrometryTransformPolynomial") {
        std::unique_ptr<AstrometryTransformPolynomial> res(new AstrometryTransformPolynomial());
        res->read(s);
        return std::move(res);
    } else
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          " astrometryTransformRead : No reader for AstrometryTransform type " + type);
}
}  // namespace jointcal
}  // namespace lsst
