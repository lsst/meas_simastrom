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

#include "pybind11/pybind11.h"
#include "pybind11/eigen.h"
#include "pybind11/stl.h"
#include "ndarray/pybind11.h"
#include "ndarray/eigen.h"
#include "Eigen/Core"

#include "lsst/utils/python.h"

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/AstrometryModel.h"
#include "lsst/jointcal/SimpleAstrometryModel.h"
#include "lsst/jointcal/ConstrainedAstrometryModel.h"
#include "lsst/jointcal/AstrometryTransform.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declareAstrometryModel(py::module &mod) {
    py::class_<AstrometryModel, std::shared_ptr<AstrometryModel>> cls(mod, "AstrometryModel");

    cls.def("getNpar", &AstrometryModel::getNpar);
    cls.def("getMapping", &AstrometryModel::getMapping, py::return_value_policy::reference_internal);
    cls.def("assignIndices", &AstrometryModel::assignIndices);
    cls.def("offsetParams", &AstrometryModel::offsetParams);
    cls.def("getSkyToTangentPlane", &AstrometryModel::getSkyToTangentPlane);
    cls.def("makeSkyWcs", &AstrometryModel::makeSkyWcs);
    cls.def("getTotalParameters", &AstrometryModel::getTotalParameters);
    cls.def("validate", &AstrometryModel::validate);
    utils::python::addOutputOp(cls, "__repr__");
    cls.def("__str__", [](AstrometryModel const &self) { return "AstrometryModel"; });
}

void declareSimpleAstrometryModel(py::module &mod) {
    py::class_<SimpleAstrometryModel, std::shared_ptr<SimpleAstrometryModel>, AstrometryModel> cls(
            mod, "SimpleAstrometryModel");

    cls.def(py::init<CcdImageList const &, const std::shared_ptr<ProjectionHandler const>, bool, unsigned,
                     unsigned>(),
            "ccdImageList"_a, "projectionHandler"_a, "initFromWcs"_a, "nNotFit"_a = 0, "order"_a = 3);

    cls.def("getTransform", &SimpleAstrometryModel::getTransform,
            py::return_value_policy::reference_internal);
    cls.def("__str__", [](SimpleAstrometryModel const &self) { return "SimpleAstrometryModel"; });
}

void declareConstrainedAstrometryModel(py::module &mod) {
    py::class_<ConstrainedAstrometryModel, std::shared_ptr<ConstrainedAstrometryModel>, AstrometryModel> cls(
            mod, "ConstrainedAstrometryModel");

    cls.def(py::init<CcdImageList const &, std::shared_ptr<ProjectionHandler const>, int, int>(),
            "ccdImageList"_a, "projectionHandler"_a, "chipOrder"_a, "visitOrder"_a);

    cls.def("getChipTransform", &ConstrainedAstrometryModel::getChipTransform,
            py::return_value_policy::reference_internal);
    cls.def("getVisitTransform", &ConstrainedAstrometryModel::getVisitTransform,
            py::return_value_policy::reference_internal);
    cls.def("__str__", [](ConstrainedAstrometryModel const &self) { return "ConstrainedAstrometryModel"; });
}

PYBIND11_MODULE(astrometryModels, mod) {
    py::module::import("lsst.jointcal.ccdImage");
    py::module::import("lsst.jointcal.astrometryTransform");
    py::module::import("lsst.jointcal.astrometryMappings");
    declareAstrometryModel(mod);
    declareSimpleAstrometryModel(mod);
    declareConstrainedAstrometryModel(mod);
}
}  // namespace
}  // namespace jointcal
}  // namespace lsst
