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

#ifndef LSST_JOINTCAL_HISTO4D_H
#define LSST_JOINTCAL_HISTO4D_H

namespace lsst {
namespace jointcal {

//! A class to histogram in 4 dimensions. Uses Sparse storage. The number of bin is limited to 256 per
//! dimension. Used in ListMatch.cc
class SparseHisto4d {
public:
    SparseHisto4d() {}
    // obvious meanings. NEntries is used as the size of the primary allocation.
    SparseHisto4d(const int n1, double min1, double max1, const int n2, double min2, double max2,
                  const int n3, double min3, double max3, const int n4, double min4, double max4,
                  const int nEntries);
    //!
    void fill(const double x[4]);
    //!
    void fill(const double x1, const double x2, const double x3, const double x4);
    //!
    int maxBin(double x[4]);

    //!
    void zeroBin(double x[4]);

    //! return the bin limits of dimension idim (0<=idim<4), around point X.
    void binLimits(const double x[4], const int idim, double &xMin, double &xMax) const;

    //!
    int getNEntries() const { return _ndata; }

    ~SparseHisto4d() {}

    // private:
    int code_value(const double x[4]) const;
    void inverse_code(const int code, double x[4]) const;
    void sort();
    void print() const;

private:
    std::vector<int> _data;
    int _ndata;
    int _dataSize;
    int _n[4];
    double _minVal[4], _maxVal[4];
    double _scale[4];
    bool _sorted;
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_HISTO4D_H
