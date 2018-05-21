"""Tests of astrometryModels (simple, constrained).

Includes tests of producing a Wcs from a model.
"""
import itertools
import os
import numpy as np

import unittest
import lsst.utils.tests

import lsst.afw.cameraGeom
import lsst.afw.geom
import lsst.afw.table
import lsst.afw.image
import lsst.afw.image.utils
import lsst.daf.persistence
import lsst.jointcal
from lsst.jointcal import astrometryModels
from lsst.meas.algorithms import astrometrySourceSelector


def getNParametersPolynomial(order):
    """Number of parameters in an astrometry polynomial model is 2 * (d+1)(d+2)/2."""
    return (order + 1)*(order + 2)


class AstrometryModelTestBase:
    @classmethod
    def setUpClass(cls):
        try:
            cls.dataDir = lsst.utils.getPackageDir('testdata_jointcal')
            # NOTE: the below is to facilitate testing with hsc test data,
            # using chips far apart on the focal plane. See the note in setup()
            # below for details. Using it requires having recently-processed
            # singleFrame output in a rerun directory in validation_data_hsc.
            # cls.dataDir = lsst.utils.getPackageDir('validation_data_hsc')
        except lsst.pex.exceptions.NotFoundError:
            raise unittest.SkipTest("testdata_jointcal not setup")

    def setUp(self):
        np.random.seed(200)

        # Append `msg` arguments to assert failures.
        self.longMessage = True
        # absolute tolerance on positional errors of 10 micro-arcsecond
        self.atol = 10.0 / (60 * 60 * 1e6)

        # Maximum difference (see assertPairsAlmostEqual) for round-trip
        # testing of the inverse for models 1 (simpler) and 2 (more.
        # Replace either one for models that don't have as accurate an inverse.
        self.inverseMaxDiff1 = 1e-5
        self.inverseMaxDiff2 = 1e-5

        self.firstIndex = 0  # for assignIndices
        matchCut = 2.0  # arcseconds
        minMeasurements = 2  # accept all star pairs.

        jointcalControl = lsst.jointcal.JointcalControl("slot_CalibFlux")
        self.associations = lsst.jointcal.Associations()
        # Work around the fact that the testdata_jointcal catalogs were produced
        # before DM-13493, and so have a different definition of the interpolated flag.
        sourceSelectorConfig = astrometrySourceSelector.AstrometrySourceSelectorConfig()
        sourceSelectorConfig.badFlags.append("base_PixelFlags_flag_interpolated")
        sourceSelector = astrometrySourceSelector.AstrometrySourceSelectorTask(config=sourceSelectorConfig)

        # Ensure that the filter list is reset for each test so that we avoid
        # confusion or contamination each time we create a cfht camera below.
        lsst.afw.image.utils.resetFilters()

        # jointcal's cfht test data has 6 ccds and 2 visits.
        inputDir = os.path.join(self.dataDir, 'cfht')
        self.visits = [849375, 850587]
        self.ccds = [12, 13, 14, 21, 22, 23]
        self.badVisit = -12345
        self.badCcd = 888

        # NOTE: the below block is to facilitate testing with validation_data_hsc.
        # NOTE: You need to have recently-processed singleFrame output available in inputDir.
        # ccd 50 is near the center, 3 is in the SW corner, 62 is on the East side.
        # ccds 101 and 103 are rotated by 90deg compared with the above.
        # inputDir = os.path.join(self.dataDir, 'DATA/rerun/20160805')
        # self.visits = [903982, 904828]  # Only need two visits for this test.
        # self.ccds = [50, 3, 62, 101, 103]

        self.butler = lsst.daf.persistence.Butler(inputDir)

        self.catalogs = []
        self.ccdImageList = []
        for (visit, ccd) in itertools.product(self.visits, self.ccds):
            dataRef = self.butler.dataRef('calexp', visit=visit, ccd=ccd)

            src = dataRef.get("src", flags=lsst.afw.table.SOURCE_IO_NO_FOOTPRINTS, immediate=True)
            goodSrc = sourceSelector.run(src)
            # Need memory contiguity to do vector-like things on the sourceCat.
            goodSrc = goodSrc.sourceCat.copy(deep=True)

            visitInfo = dataRef.get('calexp_visitInfo')
            detector = dataRef.get('calexp_detector')
            ccdId = detector.getId()
            wcs = dataRef.get('calexp_wcs')
            bbox = dataRef.get('calexp_bbox')
            filt = dataRef.get('calexp_filter')
            filterName = filt.getName()
            photoCalib = lsst.afw.image.PhotoCalib(100.0, 1.0)

            self.catalogs.append(goodSrc)
            self.associations.createCcdImage(goodSrc,
                                             wcs,
                                             visitInfo,
                                             bbox,
                                             filterName,
                                             photoCalib,
                                             detector,
                                             visit,
                                             ccdId,
                                             jointcalControl)

        # Have to set the common tangent point so projectionHandler can use skyToCTP.
        self.associations.computeCommonTangentPoint()

        self.projectionHandler = lsst.jointcal.OneTPPerVisitHandler(self.associations.getCcdImageList())

        self.associations.associateCatalogs(matchCut)
        self.associations.prepareFittedStars(minMeasurements)
        self.associations.deprojectFittedStars()

    def _prepModels(self):
        """Call this after model1 and model2 are created, to call assignIndices,
        and instantiate the fitters.
        """
        posError = 0.02  # in pixels
        # have to call this once or offsetParams will fail because the transfo indices aren't defined
        self.model1.assignIndices("Distortions", self.firstIndex)
        self.fitter1 = lsst.jointcal.AstrometryFit(self.associations, self.model1, posError)

        # have to call this once or offsetParams will fail because the transfo indices aren't defined
        self.model2.assignIndices("Distortions", self.firstIndex)
        self.fitter2 = lsst.jointcal.AstrometryFit(self.associations, self.model2, posError)

    def testMakeSkyWcsModel1(self):
        self.CheckMakeSkyWcsModel(self.model1, self.fitter1, self.inverseMaxDiff1)

    def testMakeSkyWcsModel2(self):
        self.CheckMakeSkyWcsModel(self.model2, self.fitter2, self.inverseMaxDiff2)

    def CheckMakeSkyWcsModel(self, model, fitter, inverseMaxDiff):
        """Test producing a SkyWcs on a model for every cdImage,
        both post-initialization and after one fitting step.

        Parameters
        ----------
        model : `lsst.jointcal.AstrometryModel`
            The model to test.
        fitter : `lsst.jointcal.FitterBase`
            The fitter to use to step the model to test with new (reasonable) parameters.
        inverseMaxDiff : `float`
            Required accuracy on inverse transform.
            See `lsst.afw.geom.utils.assertPairsAlmostEqual`.

        """
        # first test on as-initialized models
        for ccdImage in self.associations.getCcdImageList():
            self.checkMakeSkyWcsOneCcdImage(model, ccdImage, inverseMaxDiff)

        # now shift the models to non-default, but more reasonable, values by taking one fitting step.
        fitter.minimize("DistortionsVisit")
        fitter.minimize("Distortions")
        for ccdImage in self.associations.getCcdImageList():
            self.checkMakeSkyWcsOneCcdImage(model, ccdImage, inverseMaxDiff)

    def checkMakeSkyWcsOneCcdImage(self, model, ccdImage, inverseMaxDiff):
        """Test converting the model of one ccdImage to a SkyWcs by comparing
        to the original transfo at the tangent plane.

        Parameters
        ----------
        model : `lsst.jointcal.AstrometryModel`
            The model to test.
        ccdImage : `lsst.jointcal.CcdImage`
            The ccdImage to extract from the model and test.
        inverseMaxDiff : `float`
            Required accuracy on inverse transform.
            See `lsst.afw.geom.utils.assertPairsAlmostEqual`.
        """
        skyWcs = model.makeSkyWcs(ccdImage)
        skyToTangentPlane = model.getSky2TP(ccdImage)
        mapping = model.getMapping(ccdImage)
        posErr = 1e-3  # need something for positional error

        bbox = ccdImage.getDetector().getBBox()
        num = 200
        xx = np.linspace(bbox.getMinX(), bbox.getMaxX(), num)
        yy = np.linspace(bbox.getMinY(), bbox.getMaxY(), num)
        points = [lsst.afw.geom.Point2D(*xy) for xy in itertools.product(xx, yy)]

        expects = []
        forwards = []
        inverses = []
        spherePoints = skyWcs.pixelToSky(points)
        inverses = skyWcs.skyToPixel(skyWcs.pixelToSky(points))
        for point, spherePoint in zip(points, spherePoints):
            # TODO: Fix these "Point"s once DM-4044 is done.

            # jointcal's pixel->tangent-plane mapping
            star = lsst.jointcal.star.BaseStar(point.getX(), point.getY(), posErr, posErr)
            tpExpect = mapping.transformPosAndErrors(star)
            expects.append(lsst.afw.geom.Point2D(tpExpect.x, tpExpect.y))

            # skywcs takes pixel->sky, and we then have to go sky->tangent-plane
            onSky = lsst.jointcal.star.BaseStar(spherePoint.getLongitude().asDegrees(),
                                                spherePoint.getLatitude().asDegrees(), posErr, posErr)
            result = skyToTangentPlane.apply(onSky)
            forwards.append(lsst.afw.geom.Point2D(result.x, result.y))

        self.assertPairListsAlmostEqual(forwards, expects)
        self.assertPairListsAlmostEqual(inverses, points, maxDiff=inverseMaxDiff)


class SimpleAstrometryModelTestCase(AstrometryModelTestBase, lsst.utils.tests.TestCase):
    """Test the `SimpleAstrometryModel`, with one mapping per ccd per visit."""
    def setUp(self):
        super().setUp()
        self.order1 = 3
        self.inverseMaxDiff1 = 2e-5
        self.model1 = astrometryModels.SimpleAstrometryModel(self.associations.getCcdImageList(),
                                                             self.projectionHandler,
                                                             True,
                                                             order=self.order1)

        self.order2 = 5
        self.inverseMaxDiff2 = 5e-4
        self.model2 = astrometryModels.SimpleAstrometryModel(self.associations.getCcdImageList(),
                                                             self.projectionHandler,
                                                             False,
                                                             order=self.order2)
        self._prepModels()

    def _testGetNpar(self, model, order):
        for ccdImage in self.associations.getCcdImageList():
            result = model.getNpar(ccdImage)
            self.assertEqual(result, getNParametersPolynomial(order))

    def testGetNpar1(self):
        self._testGetNpar(self.model1, self.order1)

    def testGetNpar2(self):
        self._testGetNpar(self.model2, self.order2)

    def _testGetTotalParameters(self, model, order):
        result = model.getTotalParameters()
        expect = getNParametersPolynomial(order)*len(self.associations.getCcdImageList())
        self.assertEqual(result, expect)

    def testGetTotalParametersModel1(self):
        self._testGetTotalParameters(self.model1, self.order1)

    def testGetTotalParametersModel2(self):
        self._testGetTotalParameters(self.model2, self.order2)


class ConstrainedAstrometryModelTestCase(AstrometryModelTestBase, lsst.utils.tests.TestCase):
    """Test the `ConstrainedAstrometryModel`, with one mapping per ccd and one
    mapping per visit.
    """
    def setUp(self):
        super().setUp()
        self.visitOrder1 = 3
        self.chipOrder1 = 1
        self.inverseMaxDiff1 = 1e-5
        self.model1 = astrometryModels.ConstrainedAstrometryModel(self.associations.getCcdImageList(),
                                                                  self.projectionHandler,
                                                                  chipOrder=self.chipOrder1,
                                                                  visitOrder=self.visitOrder1)

        self.visitOrder2 = 5
        self.chipOrder2 = 2
        self.inverseMaxDiff2 = 5e-5
        self.model2 = astrometryModels.ConstrainedAstrometryModel(self.associations.getCcdImageList(),
                                                                  self.projectionHandler,
                                                                  chipOrder=self.chipOrder2,
                                                                  visitOrder=self.visitOrder2)
        self._prepModels()

        # 22 is closest to the center of the focal plane in this data, so it is not fit.
        self.fixedCcd = 22

    def _polyParams(self, chipOrder, visitOrder):
        """Number of parameters per polynomial is (d+1)(d+2)/2, summed over
        polynomials, times 2 polynomials per dimension.
        The chip transform is fixed for one chip, so only visitOrder matters
        if chipOrder is None.
        """
        params = getNParametersPolynomial(visitOrder)
        if chipOrder is not None:
            params += getNParametersPolynomial(chipOrder)
        return params

    def _testGetNpar(self, model, chipOrder, visitOrder):
        def checkParams(ccdImage, model, chipOrder, visitOrder):
            result = model.getNpar(ccdImage)
            failMsg = "ccdImage: %s, with chipOrder %s and visitOrder %s"%(ccdImage.getName(),
                                                                           chipOrder,
                                                                           visitOrder)
            self.assertEqual(result, self._polyParams(chipOrder, visitOrder), msg=failMsg)

        for ccdImage in self.associations.getCcdImageList():
            realChipOrder = None if ccdImage.getCcdId() == self.fixedCcd else chipOrder
            checkParams(ccdImage, model, realChipOrder, visitOrder)

    def testGetNpar1(self):
        self._testGetNpar(self.model1, self.chipOrder1, self.visitOrder1)

    def testGetNpar2(self):
        self._testGetNpar(self.model2, self.chipOrder2, self.visitOrder2)

    def _testGetTotalParameters(self, model, chipOrder, visitOrder):
        result = model.getTotalParameters()
        # one sensor is held fixed, hence len(ccds)-1
        expect = getNParametersPolynomial(chipOrder)*(len(self.ccds) - 1) + \
            getNParametersPolynomial(visitOrder)*len(self.visits)
        self.assertEqual(result, expect)

    def testGetTotalParametersModel1(self):
        self._testGetTotalParameters(self.model1, self.chipOrder1, self.visitOrder1)

    def testGetTotalParametersModel2(self):
        self._testGetTotalParameters(self.model2, self.chipOrder2, self.visitOrder2)

    def checkGetChipTransfo(self, model):
        # Check valid ccds
        for ccd in self.ccds:
            try:
                model.getChipTransfo(ccd)
            except lsst.pex.exceptions.wrappers.InvalidParameterError:
                self.fail("model: {} raised on ccd: {}, but should not have.".format(model, ccd))

        # Check an invalid ccd
        with self.assertRaises(lsst.pex.exceptions.wrappers.InvalidParameterError) as cm:
            model.getChipTransfo(self.badCcd)
        errMsg = "No such chipId: {} among [{}]".format(self.badCcd, ", ".join(str(ccd) for ccd in self.ccds))
        self.assertIn(errMsg, str(cm.exception))

    def testGetChipTransfo(self):
        """getChipTransfo should get each known transfo, and raise with an appropriate message otherwise."""
        self.checkGetChipTransfo(self.model1)
        self.checkGetChipTransfo(self.model2)

    def checkGetVisitTransfo(self, model):
        # Check valid visits
        for visit in self.visits:
            try:
                model.getVisitTransfo(visit)
            except lsst.pex.exceptions.wrappers.InvalidParameterError:
                self.fail("model: {} raised on visit: {}, but should not have.".format(model, visit))

        # Check an invalid visit
        with self.assertRaises(lsst.pex.exceptions.wrappers.InvalidParameterError) as cm:
            model.getVisitTransfo(self.badVisit)
        errMsg = "No such visitId: {} among [{}]".format(self.badVisit,
                                                         ", ".join(str(v) for v in self.visits))
        self.assertIn(errMsg, str(cm.exception))

    def testGetVisitTransfo(self):
        """getVisitTransfo should get each known transfo, and raise with an appropriate message otherwise."""
        self.checkGetChipTransfo(self.model1)
        self.checkGetChipTransfo(self.model2)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
