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

#ifndef LSST_JOINTCAL_CHIP_VISIT_ASTROMETRY_MAPPING_H
#define LSST_JOINTCAL_CHIP_VISIT_ASTROMETRY_MAPPING_H

#include "memory"

#include "lsst/jointcal/AstrometryMapping.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/SimpleAstrometryMapping.h"

namespace lsst {
namespace jointcal {

//! The mapping with two transforms in a row.
class ChipVisitAstrometryMapping : public AstrometryMapping {
public:
    //!
    ChipVisitAstrometryMapping(std::shared_ptr<SimpleAstrometryMapping> chipMapping,
                               std::shared_ptr<SimpleAstrometryMapping> visitMapping);

    /// No copy or move: there is only ever one instance of a given model (i.e.. per ccd+visit)
    ChipVisitAstrometryMapping(ChipVisitAstrometryMapping const &) = delete;
    ChipVisitAstrometryMapping(ChipVisitAstrometryMapping &&) = delete;
    ChipVisitAstrometryMapping &operator=(ChipVisitAstrometryMapping const &) = delete;
    ChipVisitAstrometryMapping &operator=(ChipVisitAstrometryMapping &&) = delete;

    std::size_t getNpar() const override;

    void getMappingIndices(IndexVector &indices) const override;

    //!
    void computeTransformAndDerivatives(FatPoint const &where, FatPoint &outPoint,
                                        Eigen::MatrixX2d &H) const override;
    //!
    void transformPosAndErrors(FatPoint const &where, FatPoint &outPoint) const override;

    /**
     * @copydoc AstrometryMapping::offsetParams
     *
     * @note  this routine is not used when fitting (the Model manages the mappings separately),
     *        but can be useful for debugging
     */
    void offsetParams(Eigen::VectorXd const &delta) override {
        _m1->offsetParams(delta.segment(_m1->getIndex(), _m1->getNpar()));
        _m2->offsetParams(delta.segment(_m2->getIndex() + _m1->getNpar(), _m2->getNpar()));
    }

    //! access to transforms
    AstrometryTransform const &getTransform1() const { return _m1->getTransform(); }

    //! access to transforms
    AstrometryTransform const &getTransform2() const { return _m2->getTransform(); }

    void positionDerivative(Point const &where, Eigen::Matrix2d &derivative, double epsilon) const override;

    //! Currently not implemented
    void freezeErrorTransform();

    void print(std::ostream &out) const override;

private:
    friend class ConstrainedAstrometryModel;
    //!
    void setWhatToFit(const bool fittingT1, const bool fittingT2);

    std::shared_ptr<SimpleAstrometryMapping> _m1, _m2;
    Eigen::Index _nPar1, _nPar2;
    struct tmpVars  // just there to get around constness issues
    {
        Eigen::MatrixX2d h1, h2;
        Eigen::Matrix2d dt2dx;
    } _tmpVars;

    tmpVars *tmp=&_tmpVars;
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_CHIP_VISIT_ASTROMETRY_MAPPING_H
