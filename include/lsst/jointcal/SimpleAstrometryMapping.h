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

#ifndef LSST_JOINTCAL_SIMPLE_ASTROMETRY_MAPPING_H
#define LSST_JOINTCAL_SIMPLE_ASTROMETRY_MAPPING_H

#include <memory>  // for unique_ptr

#include "lsst/jointcal/AstrometryMapping.h"
#include "lsst/jointcal/AstrometryTransform.h"

//! Class for a simple mapping implementing a generic AstrometryTransform
/*! It uses a template rather than a pointer so that the derived
classes can use the specifics of the transform. The class
simplePolyMapping overloads a few routines. */

namespace lsst {
namespace jointcal {

class SimpleAstrometryMapping : public AstrometryMapping {
public:
    SimpleAstrometryMapping(AstrometryTransform const &astrometryTransform, bool toBeFit = true)
            : toBeFit(toBeFit),
              transform(astrometryTransform.clone()),
              errorProp(transform),
              lin(new AstrometryTransformLinear) {}

    /// No copy or move: there is only ever one instance of a given mapping (i.e.. per ccd+visit)
    SimpleAstrometryMapping(SimpleAstrometryMapping const &) = delete;
    SimpleAstrometryMapping(SimpleAstrometryMapping &&) = delete;
    SimpleAstrometryMapping &operator=(SimpleAstrometryMapping const &) = delete;
    SimpleAstrometryMapping &operator=(SimpleAstrometryMapping &&) = delete;

    virtual void freezeErrorTransform() {
        // from there on, updating the transform does not change the errors.
        errorProp = transform->clone();
    }

    /// @copydoc AstrometryMapping::getNpar
    std::size_t getNpar() const override {
        if (toBeFit)
            return transform->getNpar();
        else
            return 0;
    }

    /// @copydoc AstrometryMapping::getMappingIndices
    void getMappingIndices(IndexVector &indices) const override;

    /// @copydoc AstrometryMapping::transformPosAndErrors
    void transformPosAndErrors(FatPoint const &where, FatPoint &outPoint) const override;

    /// @copydoc AstrometryMapping::positionDerivative
    void positionDerivative(Point const &where, Eigen::Matrix2d &derivative, double epsilon) const override;

    /// @copydoc AstrometryMapping::offsetParams
    void offsetParams(Eigen::VectorXd const &delta) override {
        if (toBeFit) transform->offsetParams(delta);
    }

    //! position of the parameters within the grand fitting scheme
    Eigen::Index getIndex() const { return index; }

    /// Set the index of this mapping in the grand fit.
    void setIndex(Eigen::Index i) { index = i; }

    /// @copydoc AstrometryMapping::computeTransformAndDerivatives
    virtual void computeTransformAndDerivatives(FatPoint const &where, FatPoint &outPoint,
                                                Eigen::MatrixX2d &H) const override;

    //! Access to the (fitted) transform
    virtual AstrometryTransform const &getTransform() const { return *transform; }

    /// Get whether this mapping is fit as part of a Model.
    bool getToBeFit() const { return toBeFit; }
    /// Set whether this Mapping is to be fit as part of a Model.
    void setToBeFit(bool value) { toBeFit = value; }

    void print(std::ostream &out) const override;

protected:
    // Whether this Mapping is fit as part of a Model.
    bool toBeFit;
    Eigen::Index index;
    /* inheritance may also work. Perhaps with some trouble because
       some routines in Mapping and AstrometryTransform have the same name */
    std::shared_ptr<AstrometryTransform> transform;

    std::shared_ptr<AstrometryTransform> errorProp;
    /* to avoid allocation at every call of positionDerivative.
       use a pointer for constness */
    std::unique_ptr<AstrometryTransformLinear> lin;
};

//! Mapping implementation for a polynomial transformation.
class SimplePolyMapping : public SimpleAstrometryMapping {
public:
    ~SimplePolyMapping() = default;

    // ! contructor.
    /*! The transformation will be initialized to transform, so that the effective transformation
      reads transform*CenterAndScale */
    SimplePolyMapping(AstrometryTransformLinear const &CenterAndScale,
                      AstrometryTransformPolynomial const &transform);

    /// No copy or move: there is only ever one instance of a given mapping (i.e.. per ccd+visit)
    SimplePolyMapping(SimplePolyMapping const &) = delete;
    SimplePolyMapping(SimplePolyMapping &&) = delete;
    SimplePolyMapping &operator=(SimplePolyMapping const &) = delete;
    SimplePolyMapping &operator=(SimplePolyMapping &&) = delete;

    /* The SimpleAstrometryMapping version does not account for the
       _centerAndScale transform */

    void positionDerivative(Point const &where, Eigen::Matrix2d &derivative, double epsilon) const override;

    //! Calls the transforms and implements the centering and scaling of coordinates
    /* We should put the computation of error propagation and
       parameter derivatives into the same AstrometryTransform routine because
       it could be significantly faster */
    void computeTransformAndDerivatives(FatPoint const &where, FatPoint &outPoint,
                                        Eigen::MatrixX2d &H) const override;

    /// @copydoc AstrometryMapping::transformPosAndErrors
    void transformPosAndErrors(FatPoint const &where, FatPoint &outPoint) const override;

    /// @copydoc SimpleAstrometryMapping::getTransform
    AstrometryTransform const &getTransform() const override;

private:
    /* to better condition the 2nd derivative matrix, the
    transformed coordinates are mapped (roughly) on [-1,1].
    We need both the transform and its derivative. */
    AstrometryTransformLinear _centerAndScale;
    Eigen::Matrix2d preDer;

    /* Where we store the combination. */
    mutable AstrometryTransformPolynomial actualResult;
};

#ifdef STORAGE
/*! "do nothing" mapping. The Ccdimage's that "use" this one impose the
  coordinate system */
class SimpleIdentityMapping : public SimpleAstrometryMapping<AstrometryTransformIdentity> {
public:
    //! nothing to do.
    virtual void computeTransformAndDerivatives(FatPoint const &where, FatPoint &outPoint,
                                                Eigen::MatrixX2d &H) const {
        outPoint = where;
    }
};
#endif
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_SIMPLE_ASTROMETRY_MAPPING_H
