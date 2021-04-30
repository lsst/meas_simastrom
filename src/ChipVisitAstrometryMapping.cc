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

#include "lsst/jointcal/ChipVisitAstrometryMapping.h"
#include "lsst/pex/exceptions.h"

namespace pexExcept = lsst::pex::exceptions;

namespace lsst {
namespace jointcal {

ChipVisitAstrometryMapping::ChipVisitAstrometryMapping(std::shared_ptr<SimpleAstrometryMapping> chipMapping,
                                                       std::shared_ptr<SimpleAstrometryMapping> visitMapping)
        : _m1(chipMapping), _m2(visitMapping) {
    /* Allocate the record of temporary variables, so that they are not
       allocated at every call. This is hidden behind a pointer in order
       to be allowed to alter them in a const routine. */
    //tmp = std::unique_ptr<tmpVars>(new tmpVars);
    setWhatToFit(true, true);
}

std::size_t ChipVisitAstrometryMapping::getNpar() const { return _nPar1 + _nPar2; }

void ChipVisitAstrometryMapping::getMappingIndices(IndexVector &indices) const {
    std::size_t npar = getNpar();
    if (indices.size() < npar) indices.resize(npar);
    // in case we are only fitting one of the two transforms
    if (_nPar1)
        _m1->getMappingIndices(indices);
    else if (_nPar2) {
        _m2->getMappingIndices(indices);
        return;
    }
    // if we get here we are fitting both
    // there is probably a more elegant way to feed a subpart of a std::vector
    IndexVector ind2(_nPar2);
    _m2->getMappingIndices(ind2);
    for (Eigen::Index k = 0; k < _nPar2; ++k) indices.at(k + _nPar1) = ind2.at(k);
}

void ChipVisitAstrometryMapping::computeTransformAndDerivatives(FatPoint const &where, FatPoint &outPoint,
                                                                Eigen::MatrixX2d &H) const {
    // not true in general. Will crash if H is too small.
    //  assert(H.cols()==Npar());

    FatPoint pMid;
    // don't need errors there but no Mapping::Transform() routine.

    if (_nPar1) {
        _m1->computeTransformAndDerivatives(where, pMid, tmp->h1);
        // the last argument is epsilon and is not used for polynomials
        _m2->positionDerivative(pMid, tmp->dt2dx, 1e-4);
        H.block(0, 0, _nPar1, 2) = tmp->h1 * tmp->dt2dx;
    } else
        _m1->transformPosAndErrors(where, pMid);
    if (_nPar2) {
        _m2->computeTransformAndDerivatives(pMid, outPoint, tmp->h2);
        H.block(_nPar1, 0, _nPar2, 2) = tmp->h2;
    } else
        _m2->transformPosAndErrors(pMid, outPoint);
}

/*! Sets the _nPar{1,2} and allocates H matrices accordingly, to
   avoid allocation at every call. If we did not care about dynamic
   allocation, we could just put the information of what moves and
   what doesn't into the SimpleAstrometryMapping. */
void ChipVisitAstrometryMapping::setWhatToFit(const bool fittingT1, const bool fittingT2) {
    if (fittingT1) {
        _nPar1 = _m1->getNpar();
        tmp->h1 = Eigen::MatrixX2d(_nPar1, 2);
    } else
        _nPar1 = 0;
    if (fittingT2) {
        _nPar2 = _m2->getNpar();
        tmp->h2 = Eigen::MatrixX2d(_nPar2, 2);
    } else
        _nPar2 = 0;
}

void ChipVisitAstrometryMapping::transformPosAndErrors(const FatPoint &where, FatPoint &outPoint) const {
    FatPoint pMid;
    _m1->transformPosAndErrors(where, pMid);
    _m2->transformPosAndErrors(pMid, outPoint);
}

void ChipVisitAstrometryMapping::positionDerivative(Point const &where, Eigen::Matrix2d &derivative,
                                                    double epsilon) const {
    Eigen::Matrix2d d1, d2;  // seems that it does not trigger dynamic allocation
    _m1->positionDerivative(where, d1, 1e-4);
    FatPoint pMid;
    _m1->transformPosAndErrors(where, pMid);
    _m2->positionDerivative(pMid, d2, 1e-4);
    /* The following line is not a mistake. It is a consequence
       of chosing derivative(0,1) = d(y_out)/d x_in. */
    derivative = d1 * d2;
}

void ChipVisitAstrometryMapping::freezeErrorTransform() {
    throw LSST_EXCEPT(pexExcept::TypeError,
                      " The routine ChipVisitAstrometryMapping::freezeErrorTransform() was thought to be "
                      "useless and is not implemented (yet)");
}

void ChipVisitAstrometryMapping::print(std::ostream &out) const {
    out << "Composed Astrometry Mapping:" << std::endl;
    out << "Sensor f(x,y)->(u,v): " << *_m1 << std::endl;
    out << "Visit f(u,v)->alpha,delta: " << *_m2 << std::endl;
}

}  // namespace jointcal
}  // namespace lsst
