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
OUTPUT_DOA_FILE = 'local_anisotropy.vol'
OUTPUT_VARIANCE_FILE = 'resolution_variance.vol'
OUTPUT_DIRECTIONS_FILE = 'preffered.vol'


class XmippProtMonoDir(ProtAnalysis3D):
    """    
    Given a map the protocol assigns local resolutions to each voxel of the map.
    """
    _label = 'directional ResDir'
    _lastUpdateVersion = VERSION_1_1
    
    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)

    
    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolumes', PointerParam, pointerClass='Volume',
                      label="Input Volume", important=True,
                      help='Select a volume for determining its local resolution.')

        form.addParam('Mask', PointerParam, pointerClass='VolumeMask',
                      label="Binary Mask", allowsNull=True,
                      help='The mask determines which points are specimen and which ones not')

        group = form.addGroup('Extra parameters')
        group.addParam('symmetry', StringParam, default='c1',
                      label="Symmetry",
                      help='Symmetry group. By default = c1.'
                      'See [[http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry][Symmetry]]'
                      'for a description of the symmetry groups format, If no symmetry is present, give c1.')
        
        group.addParam('angularsampling', FloatParam, default=15, expertLevel=LEVEL_ADVANCED,
                      label="Angular Sampling",
                      help='Angular sampling to cover the projection sphere')
        
        group.addParam('significance', FloatParam, default=0.95, expertLevel=LEVEL_ADVANCED,
                      label="Significance",
                      help='Relution is computed using hipothesis tests, this value determines'
                      'the significance of that test')
        
        group.addParam('isPremasked', BooleanParam, default=False,
                      label="Is the original premasked?",
                      help='Sometimes the original volume is masked inside a spherical mask. In this case'
                      'please select yes')
        
        group.addParam('volumeRadius', FloatParam, default=-1,
                      label="Spherical mask radius (px)",
                      condition = 'isPremasked', 
                      help='When the original volume is originally premasked, the noise estimation ought'
                      'to be performed inside that premask, and out of the provieded mask asked in the previus'
                      'box. The radius value, determines the radius of the spherical premask. By default'
                      'radius = -1 use the half of the volume size as radius')

        group.addParam('stepSize', FloatParam, allowsNull=True,
                      expertLevel=LEVEL_ADVANCED, label='Step')


    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        
        self.micsFn = self._getPath()

        self.vol0Fn = self.inputVolumes.get().getFileName()
        self.maskFn = self.Mask.get().getFileName()

        # Convert input into xmipp Metadata format
        convertId = self._insertFunctionStep('convertInputStep', )
        
        MS = self._insertFunctionStep('resolutionMonogenicSignalStep',
                                      prerequisites=[convertId])

        self._insertFunctionStep('createOutputStep', prerequisites=[MS])

        self._insertFunctionStep("createHistrogram")
        
        
    def convertInputStep(self):
        """ Read the input volume.
        """
        extVol0 = getExt(self.vol0Fn)
        if (extVol0 == '.mrc') or (extVol0 == '.map'):
            self.vol0Fn = self.vol0Fn + ':mrc'

        extMask = getExt(self.maskFn)
        if ((extMask == '.mrc') or (extMask == '.map')):
            self.maskFn = self.maskFn + ':mrc'


    def resolutionMonogenicSignalStep(self):

        if self.isPremasked:
            if self.volumeRadius == -1:
                xdim, _ydim, _zdim = self.inputVolumes.get().getDim()
                xdim = xdim*0.5
            else:
                xdim = self.volumeRadius.get()
        else:
            xdim, _ydim, _zdim = self.inputVolumes.get().getDim()
            xdim = xdim*0.5

        # Number of frequencies
        #Nfreqs = xdim
                
        params = ' --vol %s' % self.vol0Fn
        params += ' --mask %s' % self.maskFn
        params += ' -o %s' % self._getExtraPath(OUTPUT_RESOLUTION_FILE)
        params += ' --sampling_rate %f' % self.inputVolumes.get().getSamplingRate()
        params += ' --angular_sampling %f' % self.angularsampling.get()
        #params += ' --number_frequencies %f' % Nfreqs
        params += ' --varVol %s' % self._getExtraPath(OUTPUT_VARIANCE_FILE)
        params += ' --volumeRadius %f' % xdim
        params += ' --sym %s' % self.symmetry.get()
        params += ' --significance %f' % self.significance.get()
        params += ' --md_resdir %s' % self._getExtraPath(METADATA_ANGLES_FILE)
        params += ' --doa_vol %s' % self._getExtraPath(OUTPUT_DOA_FILE)
        params += ' --directions %s' % self._getExtraPath(OUTPUT_DIRECTIONS_FILE)

        self.runJob('xmipp_resolution_directional', params)

    def createHistrogram(self):

        params = ' -i %s' % self._getExtraPath(OUTPUT_DOA_FILE)
        params += ' --mask binary_file %s' % self.maskFn
        params += ' --steps %f' % 30
        params += ' -o %s' % self._getExtraPath('hist_DoA.xmd')
        params += ' --range %f %f' % (0, 1)#(self.minRes.get(), self.maxRes.get())
        
        self.runJob('xmipp_image_histogram', params)
        
        
    def createOutputStep(self):
        volume_path_doa = self._getExtraPath(OUTPUT_DOA_FILE)
        volume_path_var = self._getExtraPath(OUTPUT_VARIANCE_FILE)
        
        self.volumesSet_doa = self._createSetOfVolumes('doaVol')
        self.volumesSet_var = self._createSetOfVolumes('varianceVol')
        
        self.volumesSet_doa.setSamplingRate(self.inputVolumes.get().getSamplingRate())
        self.volumesSet_var.setSamplingRate(self.inputVolumes.get().getSamplingRate())

        readSetOfVolumes(volume_path_doa, self.volumesSet_doa)
        readSetOfVolumes(volume_path_var, self.volumesSet_var)
        self._defineOutputs(outputVolume_doa=self.volumesSet_doa)
        self._defineOutputs(outputVolume_var=self.volumesSet_var)
        
        self._defineSourceRelation(self.inputVolumes, self.volumesSet_doa)
        self._defineSourceRelation(self.inputVolumes, self.volumesSet_var)

    # --------------------------- INFO functions ------------------------------

    def _methods(self):
        messages = []
        if hasattr(self, 'doaVol'):
            messages.append(
                'Information about the method/article in ' + MONORES_METHOD_URL)
        return messages
    
    def _summary(self):
        summary = []

        return summary

    def _citations(self):
        return ['Not yet']

