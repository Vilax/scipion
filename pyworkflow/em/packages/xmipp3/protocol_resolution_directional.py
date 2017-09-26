# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Jose Luis Vilas (jlvilas@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
from pyworkflow import VERSION_1_1
from pyworkflow.protocol.params import (PointerParam, StringParam, BooleanParam, FloatParam, LEVEL_ADVANCED)
from pyworkflow.em.protocol.protocol_3d import ProtAnalysis3D
from convert import readSetOfVolumes
from pyworkflow.object import Float
from pyworkflow.em import ImageHandler
from pyworkflow.utils import getExt
import numpy as np
import pyworkflow.em.metadata as md


MONORES_METHOD_URL = 'http://github.com/I2PC/scipion/wiki/XmippProtMonoDir'

OUTPUT_RESOLUTION_FILE = 'mgresolution.vol'
OUTPUT_RESOLUTION_FILE_CHIMERA = 'MG_Chimera_resolution.vol'
OUTPUT_MASK_FILE = 'output_Mask.vol'
FN_MEAN_VOL = 'mean_volume.vol'
METADATA_ANGLES_FILE = 'angles_md.xmd'
OUTPUT_RESOLUTION_MAX_FILE = 'maxResolution.vol'
OUTPUT_RESOLUTION_MIN_FILE = 'minResolution.vol'
OUTPUT_DOA_FILE = 'local_anisotropy.vol'
OUTPUT_VARIANCE_FILE = 'resolution_variance.vol'


class XmippProtMonoDir(ProtAnalysis3D):
    """    
    Given a map the protocol assigns local resolutions to each voxel of the map.
    """
    _label = 'directional ResDir'
    _lastUpdateVersion = VERSION_1_1
    
    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)
        self.min_res_init = Float() 
        self.max_res_init = Float() 
        
    
    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('halfVolumes', BooleanParam, default=False,
                      label="Would you like to use half volumes?",
                      help='The noise estimation for determining the local resolution '
                           'is performed via half volumes.')

        form.addParam('inputVolumes', PointerParam, pointerClass='Volume',
                      label="Input Volume", important=True,
                      condition = 'not halfVolumes',
                      help='Select a volume for determining its local resolution.')

        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Volume Half 1", important=True,
                      condition = 'halfVolumes', 
                      help='Select a volume for determining its local resolution.')

        form.addParam('inputVolume2', PointerParam, pointerClass='Volume',
                      label="Volume Half 2", important=True,
                      condition='halfVolumes',
                      help='Select a second volume for determining a local resolution.')

        form.addParam('Mask', PointerParam, pointerClass='VolumeMask', 
                      condition='(halfVolumes) or (not halfVolumes)',
                      label="Binary Mask", allowsNull=True,
                      help='The mask determines which points are specimen and which ones not')

        group = form.addGroup('Extra parameters')
        group.addParam('symmetry', StringParam, default='c1',
                      label="Symmetry",
                      help='Symmetry group. By default = c1.'
                      'See [[http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry][Symmetry]]'
                      'for a description of the symmetry groups format, If no symmetry is present, give c1.')

        line = group.addLine('Resolution Range (A)',
                            help="If the user knows the range of resolutions or only a"
                                 " range of frequency needs to be analysed")
        
        group.addParam('samplingPoints', FloatParam, default=50, 
                      label="N. Directions",
                      help='Number of directions to be computed')
        
        
        group.addParam('significance', FloatParam, default=0.95, expertLevel=LEVEL_ADVANCED,
                      label="Significance",
                      help='Relution is computed using hipothesis tests, this value determines'
                      'the significance of that test')
        
        group.addParam('isPremasked', BooleanParam, default=False,
                      label="Is the original premasked?",
                      help='Sometimes the original volume is masked inside a spherical mask. In this case'
                      'please select yes')
        
        form.addParam('noiseonlyinhalves', BooleanParam, expertLevel=LEVEL_ADVANCED,
                      default=True,
                      label="Use noise outside the mask?", 
                      condition = 'halfVolumes',
                      help='Select yes if the volume present noise outside the mask.'
                      ' Otherwise, select No.')
        
        group.addParam('volumeRadius', FloatParam, default=-1,
                      label="Spherical mask radius (px)",
                      condition = 'isPremasked and not halfVolumes', 
                      help='When the original volume is originally premasked, the noise estimation ought'
                      'to be performed inside that premask, and out of the provieded mask asked in the previus'
                      'box. The radius value, determines the radius of the spherical premask. By default'
                      'radius = -1 use the half of the volume size as radius')
        
        group.addParam('volumeRadiusHalf', FloatParam, default=-1,
                      label="Spherical mask radius (px)",
                      condition = 'halfVolumes and isPremasked',
                      help='When the origianl volume is originally premasked,'
                      'the noise estimation ought to be performed inside that'
                      'premask, and out of the provieded mask asked in the previus'
                      'box. The radius value, determines the radius in pixels of '
                      'the spherical premask. By default radius = -1 use the half'
                      'of the volume size as radius')

        line.addParam('minRes', FloatParam, default=1, label='High')
        line.addParam('maxRes', FloatParam, default=30, label='Low')
        line.addParam('stepSize', FloatParam, allowsNull=True,
                      expertLevel=LEVEL_ADVANCED, label='Step')


    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        
        self.micsFn = self._getPath()

        if self.halfVolumes:
            self.vol1Fn = self.inputVolume.get().getFileName()
            self.vol2Fn = self.inputVolume2.get().getFileName()
            self.maskFn = self.Mask.get().getFileName()

            self.inputVolumes.set(None)

        else:
            self.vol0Fn = self.inputVolumes.get().getFileName()
            self.maskFn = self.Mask.get().getFileName()
            self.inputVolume.set(None)
            self.inputVolume2.set(None)

            # Convert input into xmipp Metadata format
        convertId = self._insertFunctionStep('convertInputStep', )
        
        MS = self._insertFunctionStep('resolutionMonogenicSignalStep',
                                      prerequisites=[convertId])

        self._insertFunctionStep('createOutputStep', prerequisites=[MS])

        #self._insertFunctionStep("createHistrogram")

    def convertInputStep(self):
        """ Read the input volume.
        """
        if (self.halfVolumes.get() is False):
            extVol0 = getExt(self.vol0Fn)
            if (extVol0 == '.mrc') or (extVol0 == '.map'):
                self.vol0Fn = self.vol0Fn + ':mrc'

        if self.halfVolumes.get() is True:
            extVol1 = getExt(self.vol1Fn)
            extVol2 = getExt(self.vol2Fn)
            if (extVol1 == '.mrc') or (extVol1 == '.map'):
                self.vol1Fn = self.vol1Fn + ':mrc'
            if (extVol2 == '.mrc') or (extVol2 == '.map'):
                self.vol2Fn = self.vol2Fn + ':mrc'

        extMask = getExt(self.maskFn)
        if ((extMask == '.mrc') or (extMask == '.map')):
            self.maskFn = self.maskFn + ':mrc'


    def resolutionMonogenicSignalStep(self):

        # Number of frequencies
        if self.stepSize.hasValue():
            Nfreqs = round((self.maxRes.get() - self.minRes.get())/self.stepSize.get())
        else:
            Nfreqs = 50
  
        if self.halfVolumes:
            if self.isPremasked:
                if self.volumeRadiusHalf == -1:
                    xdim, _ydim, _zdim = self.inputVolume.get().getDim()
                    xdim = xdim*0.5
                else:
                    xdim = self.volumeRadiusHalf.get()
            else:
                xdim, _ydim, _zdim = self.inputVolume.get().getDim()
                xdim = xdim*0.5
        else:
            if self.isPremasked:
                if self.volumeRadius == -1:
                    xdim, _ydim, _zdim = self.inputVolumes.get().getDim()
                    xdim = xdim*0.5
                else:
                    xdim = self.volumeRadius.get()
            else:
                xdim, _ydim, _zdim = self.inputVolumes.get().getDim()
                xdim = xdim*0.5
                
        if self.halfVolumes.get() is False:
            params = ' --vol %s' % self.vol0Fn
            params += ' --mask %s' % self.maskFn
        else:
            params = ' --vol %s' % self.vol1Fn
            params += ' --vol2 %s' % self.vol2Fn
            params += ' --meanVol %s' % self._getExtraPath(FN_MEAN_VOL)
            params += ' --mask %s' % self.maskFn
        if self.Mask.hasValue() is False:
            params += ' --nonmanual_mask '
        params += ' --mask_out %s' % self._getExtraPath(OUTPUT_MASK_FILE)
        params += ' -o %s' % self._getExtraPath(OUTPUT_RESOLUTION_FILE)
        if (self.halfVolumes):
            params += ' --sampling_rate %f' % self.inputVolume.get().getSamplingRate()
            if (self.noiseonlyinhalves is False):
                params += ' --noiseonlyinhalves'
        else:
            params += ' --sampling_rate %f' % self.inputVolumes.get().getSamplingRate()
        params += ' --sampling_points %f' % self.samplingPoints.get()
        params += ' --number_frequencies %f' % Nfreqs
        params += ' --minRes %f' % self.minRes.get()
        params += ' --maxRes %f' % self.maxRes.get()
        params += ' --maxVol %s' % self._getExtraPath(OUTPUT_RESOLUTION_MAX_FILE)
        params += ' --minVol %s' % self._getExtraPath(OUTPUT_RESOLUTION_MIN_FILE)
        params += ' --varVol %s' % self._getExtraPath(OUTPUT_VARIANCE_FILE)
        params += ' --volumeRadius %f' % xdim
        params += ' --chimera_volume %s' % self._getExtraPath(OUTPUT_RESOLUTION_FILE_CHIMERA)
        params += ' --sym %s' % self.symmetry.get()
        params += ' --significance %f' % self.significance.get()
        params += ' --md_resdir %s' % self._getExtraPath(METADATA_ANGLES_FILE)
        params += ' --doa_vol %s' % self._getExtraPath(OUTPUT_DOA_FILE)

        self.runJob('xmipp_resolution_directional', params)


    def createHistrogram(self):

        params = ' -i %s' % self._getExtraPath(OUTPUT_DOA_FILE)
        params += ' --mask binary_file %s' % self._getExtraPath(OUTPUT_MASK_FILE)
        params += ' --steps %f' % 30
        params += ' --range %f %f' % (self.min_res_init, self.max_res_init)#(self.minRes.get(), self.maxRes.get())
        params += ' -o %s' % self._getExtraPath('hist_DoA.xmd')

        self.runJob('xmipp_image_histogram', params)
        
        
    def createOutputStep(self):
        volume_path_max = self._getExtraPath(OUTPUT_RESOLUTION_MAX_FILE)
        volume_path_min = self._getExtraPath(OUTPUT_RESOLUTION_MIN_FILE)
        volume_path_doa = self._getExtraPath(OUTPUT_DOA_FILE)
        volume_path_var = self._getExtraPath(OUTPUT_VARIANCE_FILE)
        
        self.volumesSet_max = self._createSetOfVolumes('resolutionMaxVol')
        self.volumesSet_min = self._createSetOfVolumes('resolutionMinVol')
        self.volumesSet_doa = self._createSetOfVolumes('doaVol')
        self.volumesSet_var = self._createSetOfVolumes('varianceVol')
        
        if (self.halfVolumes):
            self.volumesSet_min.setSamplingRate(self.inputVolume.get().getSamplingRate())
            self.volumesSet_max.setSamplingRate(self.inputVolume.get().getSamplingRate())
            self.volumesSet_doa.setSamplingRate(self.inputVolume.get().getSamplingRate())
            self.volumesSet_var.setSamplingRate(self.inputVolume.get().getSamplingRate())
            
        else:
            self.volumesSet_min.setSamplingRate(self.inputVolumes.get().getSamplingRate())
            self.volumesSet_max.setSamplingRate(self.inputVolumes.get().getSamplingRate())
            self.volumesSet_doa.setSamplingRate(self.inputVolumes.get().getSamplingRate())
            self.volumesSet_var.setSamplingRate(self.inputVolumes.get().getSamplingRate())

        readSetOfVolumes(volume_path_min, self.volumesSet_min)
        readSetOfVolumes(volume_path_max, self.volumesSet_max)
        readSetOfVolumes(volume_path_doa, self.volumesSet_doa)
        readSetOfVolumes(volume_path_var, self.volumesSet_var)
        self._defineOutputs(outputVolume_max=self.volumesSet_max)
        self._defineOutputs(outputVolume_min=self.volumesSet_min)
        self._defineOutputs(outputVolume_doa=self.volumesSet_doa)
        self._defineOutputs(outputVolume_var=self.volumesSet_var)
        
        if (self.halfVolumes):
            self._defineSourceRelation(self.inputVolume, self.volumesSet_max)
            self._defineSourceRelation(self.inputVolume, self.volumesSet_min)
            self._defineSourceRelation(self.inputVolume, self.volumesSet_doa)
            self._defineSourceRelation(self.inputVolume, self.volumesSet_var)
        else:
            self._defineSourceRelation(self.inputVolumes, self.volumesSet_max)
            self._defineSourceRelation(self.inputVolumes, self.volumesSet_min)
            self._defineSourceRelation(self.inputVolumes, self.volumesSet_doa)
            self._defineSourceRelation(self.inputVolumes, self.volumesSet_var)

    # --------------------------- INFO functions ------------------------------

    def _methods(self):
        messages = []
        if hasattr(self, 'outputVolume'):
            messages.append(
                'Information about the method/article in ' + MONORES_METHOD_URL)
        return messages
    
    def _summary(self):
        summary = []

        return summary

    def _citations(self):
        return ['Vilas2017']

