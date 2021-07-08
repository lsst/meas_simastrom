#!/usr/bin/env python
# This file is part of jointcal.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
Use jointcal output to plot distortion as characterized by the radial and
tangential pixel scale as a function of field angle.
"""

import numpy as np

import matplotlib
matplotlib.use('Agg')  # noqa: E402
import matplotlib.pyplot as plt

from lsst.daf.persistence import Butler
from lsst.daf.persistence.butlerExceptions import NoResults
from lsst.geom import Point2D
from lsst.afw.cameraGeom import FOCAL_PLANE


def getPixelScale(wcs, focalRadius):
    """Compute pixel scale in radial and tangential directions as function
    field angle.

    Parameters
    ----------
    wcs : `lsst.afw.geom.SkyWCS`
    focalRadius : `float`
        Maximum radius in FOCAL_PLANE units for which to calculate.

    Returns
    -------
    fieldAngle : `numpy.ndarray`
        Array of field angles in degrees
    radialScale : `numpy.ndarray`
        Array of radial direction pixel scales in arcseconds.
    tangentialScale : `numpy.ndarray`
        Array of tangential direction pixel scales in arcseconds.

    Notes
    -----
    Pixel scales are calculated from finite differences only along the +y
    focal plane direction.
    """
    mapping = wcs.getFrameDict().getMapping('FOCAL', 'PIXELS')
    trans = wcs.getTransform()  # Pixels to Sky as Point2d -> SpherePoint

    boresight = trans.applyForward(Point2D(mapping.applyForward([0, 0])))
    rs = np.linspace(0, focalRadius, 100)
    fieldAngle = []
    radialScale = []
    tangentialScale = []
    for r in rs:
        sp1 = trans.applyForward(Point2D(mapping.applyForward([0, r])))
        sp2 = trans.applyForward(Point2D(mapping.applyForward([0, r+1])))
        sp3 = trans.applyForward(Point2D(mapping.applyForward([1, r])))
        radialScale.append(sp1.separation(sp2).asArcseconds())
        tangentialScale.append(sp1.separation(sp3).asArcseconds())
        fieldAngle.append(boresight.separation(sp1).asDegrees())
    return (np.array(fieldAngle),
            np.array(radialScale),
            np.array(tangentialScale))


def getFocalRadius(camera):
    """Determine maximum field radius of this camera.

    Parameters
    ----------
    camera : `lsst.afw.cameraGeom.Camera`

    Returns
    -------
    focalRadius : `float`
        Maximum focal plane radius in FOCAL_PLANE units.
    """
    radii = []
    for ccd in camera:
        for corner in ccd.getCorners(FOCAL_PLANE):
            radii.append(np.hypot(*corner))
    return np.max(radii)


def main():
    """Commandline entry point for plot_pixelScale.py"""
    import argparse
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("repo", type=str,
                        help="Directory to butler repository "
                             "containing `jointcal_wcs` datasets.")
    parser.add_argument("tract", default=None, type=int,
                        help="Tract to plot.")
    parser.add_argument("--visits", type=int, nargs='+',
                        help="List of visits to process.")
    parser.add_argument("--degree", type=int, default=9,
                        help="Order of polynomial fit.")
    args = parser.parse_args()

    butler = Butler(inputs=args.repo)
    camera = butler.get("camera")
    focalRadius = getFocalRadius(camera)

    results = []
    for visit in args.visits:
        # Results are independent of ccd, so just find the first one that works.
        for ccd in range(250):  # inelegant, but effective and reasonably future-proof?
            try:
                wcs = butler.get("jointcal_wcs", ccd=ccd, visit=visit, tract=args.tract)
            except NoResults:
                continue
            results.append(getPixelScale(wcs, focalRadius))
            break

    results = np.array(results)
    field = results[:, 0]
    tan = results[:, 1]
    rad = results[:, 2]

    meanField = np.mean(field, axis=0)
    meanTan = np.mean(tan, axis=0)
    meanRad = np.mean(rad, axis=0)

    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(5, 5), sharex=True)

    tanLines = axes[0].plot(meanField, tan.T, c='r', lw=0.5, label='tangential')
    radLines = axes[0].plot(meanField, rad.T, c='b', lw=0.5, label='radial')
    axes[1].plot(meanField, 1e6*(tan-meanTan).T, c='r', lw=0.5)
    axes[1].plot(meanField, 1e6*(rad-meanRad).T, c='b', lw=0.5)

    axes[0].set_ylabel("pixel size (arcsec)")
    axes[1].set_ylabel("residual from mean (micro arcsec)")
    axes[1].set_xlabel("field angle (deg)")
    axes[0].legend(handles=[tanLines[0], radLines[0]])
    fig.tight_layout(h_pad=0)
    fig.savefig(f"pixelScale-tract_{args.tract}.png", dpi=300)

    print()
    print("tangential fit coefficients")
    print(np.polyfit(
        meanField,
        meanTan,
        args.degree
    ))

    print("radial fit coefficients")
    print(np.polyfit(
        meanField,
        meanRad,
        args.degree
    ))


if __name__ == "__main__":
    main()
