from __future__ import division, absolute_import
#
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2011, 2012 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import numpy

import lsst.pex.config as pexConfig
import lsst.coadd.utils as coaddUtils
import lsst.pipe.base as pipeBase
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable

from lsst.afw.fits import FitsError
from lsst.pipe.tasks.selectImages import WcsSelectImagesTask, SelectStruct
from lsst.coadd.utils import CoaddDataIdContainer
from lsst.pipe.tasks.getRepositoryData import DataRefListRunner

from .simastromLib import test, test2, simAstrom

__all__ = ["SimAstromConfig", "SimAstromTask"]

class SimAstromConfig(pexConfig.Config):
    """Config for SimAstromTask
    """

# Keep this config parameter as a place holder
    doWrite = pexConfig.Field(
        doc = "persist SimAstrom output...",
        dtype = bool,
        default = True,
    )

class SimAstromTask(pipeBase.CmdLineTask):
 
    ConfigClass = SimAstromConfig
    RunnerClass = DataRefListRunner
    _DefaultName = "simAstrom"
    
    def __init__(self, *args, **kwargs):
        pipeBase.Task.__init__(self, *args, **kwargs)
#        self.makeSubtask("select")

# We don't need to persist config and metadata at this stage. In this way, we don't need to put a specific entry in the
# camera mapper policy file
    def _getConfigName(self):
        return None
        
    def _getMetadataName(self):
        return None
        
    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser
        """
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)

        parser.add_id_argument("--id", "calexp", help="data ID, e.g. --selectId visit=6789 ccd=0..9")
        return parser

    @pipeBase.timeMethod
    def run(self, ref):
        
#        sourceCat = test2()
        
        srcList = []
        metaList = []
        wcsList = []
        
        for dataRef in ref :
            src = dataRef.get("src", immediate=False)
            calexp = dataRef.get("calexp", immediate=True)
            wcs = calexp.getWcs()
            md = dataRef.get("calexp_md", immediate=False)
            calib = afwImage.Calib(md)
            
            config = StarSelectorConfig()
            ss = StarSelector(config)
            newSrc = ss.select(src, calib)
            print len(newSrc)
            
        # Should call a source selector here in order to send a list
        # of reasonable star to the fitter.
            srcList.append(newSrc)
            
            metaList.append(md)
            wcsList.append(wcs)
            
        simA = simAstrom(srcList, metaList, wcsList)

class StarSelectorConfig(pexConfig.Config):
    
    badFlags = pexConfig.ListField(
        doc = "List of flags which cause a source to be rejected as bad",
        dtype = str,
        default = [ "base_PixelFlags_flag_saturated", 
                    "base_PixelFlags_flag_cr",
                    "base_PixelFlags_flag_interpolated",
                    "base_PsfFlux_flag_edge", 
                    "base_SdssCentroid_flag"],
    )

class StarSelector(object) :
    
    ConfigClass = StarSelectorConfig

    def __init__(self, config):
        """Construct a star selector
        
        @param[in] config: An instance of StarSelectorConfig
        """
        self.config = config
    
    def select(self, srcCat, calib):
# Return a catalog containing only reasonnable stars

        schema = srcCat.getSchema()
        newCat = afwTable.SourceCatalog(schema)
        for src in srcCat :
            # Reject galaxies
            if src.get("base_ClassificationExtendedness_value") > 0.5 :
                continue
            # Do not consider sources with bad flags
            for f in self.config.badFlags :
                rej = 0
                if src.get(f) :
                    rej = 1
                    break
            if rej == 1 :
                continue
            # Reject negative flux
            flux = src.get('base_PsfFlux_flux')
            if flux < 0 :
                continue
            # Reject object with magnitude > 19
            if calib.getMagnitude(flux) > 19 :
                continue
                
            newCat.append(src)
            
        print len(srcCat), len(newCat)
        
        return newCat
