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

from os.path import join
from pyworkflow.object import Float, String
from pyworkflow.protocol.params import (PointerParam, FloatParam, STEPS_PARALLEL,
                                        StringParam, EnumParam, LEVEL_ADVANCED)
from pyworkflow.em.data import Volume
from pyworkflow.em.protocol import ProtAnalysis3D
from pyworkflow.utils.path import moveFile, makePath
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles, writeSetOfClasses2D
from pyworkflow.em.metadata import MetaData, getBlocksInMetaDataFile



class XmippProtValidateTiltPairs(ProtAnalysis3D):
    """    
    Ranks a set of volumes according to their alignment reliability obtained from a clusterability test.
    """

    _label = 'validate_tiltpairs'
    WEB = 0

    
    def __init__(self, *args, **kwargs):
        ProtAnalysis3D.__init__(self, *args, **kwargs)
        
        if (self.WEB == 1):
            self.stepsExecutionMode = STEPS_PARALLEL
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputVolume', PointerParam, pointerClass='SetOfVolumes, Volume',
                      label="Input volume",  
                      help='Select the input volume.')
        form.addParam('inputTiltPairs', PointerParam, pointerClass='ParticlesTiltPair', 
                      label="Input particles tilt pairs",  
                      help='Select the input projection images .')      
        form.addParam('inputClasses', PointerParam, pointerClass='SetOfClasses2D', 
                      label="Untilted Classes",  
                      help='Select the input projection images .') 
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group", 
                      help='See [[Xmipp Symmetry][http://www2.mrc-lmb.cam.ac.uk/Xmipp/index.php/Conventions_%26_File_formats#Symmetry]] page '
                           'for a description of the symmetry format accepted by Xmipp') 
        form.addParam('angularSampling', FloatParam, default=5, expertLevel=LEVEL_ADVANCED,
                      label="Angular Sampling (degrees)",
                      help='Angular distance (in degrees) between neighboring projection points ')
        form.addParallelSection(threads=0, mpi=4)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):        
        self.classes = self.inputClasses.get()
        self.untilt = self.inputTiltPairs.get().getUntilted()
        self.tilt = self.inputTiltPairs.get().getTilted()
        
        
        convertId = self._insertFunctionStep('convertInputStep', self.classes, self.untilt, self.tilt)
        
        volStepId = self._insertFunctionStep('validationStep')
              
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def convertInputStep(self, classes, untilt, tilt):
        classesFn = self._getTmpPath('input_classes.xmd')
        untiltFn = self._getExtraPath('images_untilted.xmd')
        tiltFn = self._getExtraPath('images_tilted.xmd')
        
        writeSetOfClasses2D(classes, classesFn)
        writeSetOfParticles(untilt, untiltFn)
        writeSetOfParticles(tilt, tiltFn)
        
        self.classes2metadata(classesFn, tiltFn)
        
        

    def classes2metadata(self, classesFn, tiltFn):
        #It is compared always with tilt due to the input is untilted classes
        #The output is a metadata with correspondence untilted particles
        mdClasses = MetaData()
        mdOneBlock = MetaData()
        
        mdClasses.read(classesFn)
        mdBlocks = getBlocksInMetaDataFile(classesFn)

        for mdBlock in mdBlocks:
            if mdBlock =='classes':
                continue
            print mdBlock
            mdOneBlock.read(mdBlock + "@" + classesFn)
            fnblock = self._getExtraPath(mdBlock+"_untilted.xmd")
            mdOneBlock.write(fnblock)
            params =  ' -i %s' % (tiltFn)
            params += ' --set intersection %s '"'itemId'"' '"'itemId'"' ' % (fnblock)
            params += ' -o %s' % (self._getExtraPath(mdBlock+"_tilted.xmd"))
            
            self.runJob('xmipp_metadata_utilities', params)
            print ' '


    def validationStep(self):
        classesFn = self._getTmpPath('input_classes.xmd')
        mdClasses = MetaData()
        mdClasses.read(classesFn)

        mdBlocks = getBlocksInMetaDataFile(classesFn)
        
        for mdBlock in mdBlocks:
            if mdBlock == 'classes':
                continue
            
            params =  ' --untilt %s' %(self._getExtraPath(mdBlock+"_untilted.xmd"))
            params += ' --tilt %s' %(self._getExtraPath(mdBlock+"_tilted.xmd"))
            params += ' --vol %s' %(self.inputVolume.get().getFileName())
            params += ' --angular_sampling %f' %self.angularSampling.get()
#             params += ' --sym %s' %(self.sym.get())
            params += ' -o %s' %(mdBlock+'assigned_out.xmd')
#             params += ' --maxshift %f' % (self.maxshift.get())
            
             
            self.runJob('xmipp_validation_tilt_pairs_new', params)
    
    def createOutputStep(self):
        outputVols = self._createSetOfVolumes()
        
      
        outputVols.setSamplingRate(self.partSet.getSamplingRate())
        self._defineOutputs(outputVolumes=outputVols)
        self._defineTransformRelation(self.inputVolumes, outputVols)
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        validateMsgs = []
        # if there are Volume references, it cannot be empty.
        if not self.inputVolume.get().hasValue():
            validateMsgs.append('Please provide an input reference volume.')
        if not self.inputTiltPairs.get().hasValue():
            validateMsgs.append('Please provide input particles.')            
        return validateMsgs
    
    def _summary(self):
        summary = []

        if  (not hasattr(self,'outputVolumes')):
            summary.append("Output volumes not ready yet.")
        else:
            size = 0
            for i, vol in enumerate(self._iterInputVols()):
                size +=1
            summary.append("Volumes to validate: *%d* " % size)
            summary.append("Angular sampling: %s" % self.angularSampling.get())


        return summary
    
    def _methods(self):
        messages = []
        if (hasattr(self,'outputVolumes')):
            messages.append('The quality parameter(s) has been obtained using the approach [Vargas2014a] with angular sampling of %f and significant value of %f' % (self.angularSampling.get(), self.alpha.get()))
        return messages
    
    def _citations(self):
        return ['Vargas2014a']
    