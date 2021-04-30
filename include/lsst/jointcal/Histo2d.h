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

#ifndef LSST_JOINTCAL_HISTO2D_H
#define LSST_JOINTCAL_HISTO2D_H

#include <vector>

namespace lsst {
namespace jointcal {

class Histo2d {
public:
    Histo2d() : data() {}
    Histo2d(int nx, float minx, float maxx, int ny, float miny, float maxy);

    Histo2d(const Histo2d &other);

    void fill(float x, float y, float weight = 1.);

    double maxBin(double &x, double &y) const;

    void binWidth(double &Hdx, double &Hdy) const {
        Hdx = 1. / scalex;
        Hdy = 1. / scaley;
    }

    double binContent(double x, double y) const;

    void zeroBin(double x, double y);

    ~Histo2d() {}

private:
    void operator=(const Histo2d &right);
    bool indices(double x, double y, int &ix, int &iy) const;

    std::vector<float> data;
    int nx, ny;
    float minx, miny;
    float scalex, scaley;
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_HISTO2D_H
