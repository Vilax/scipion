# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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

from math import floor, ceil
import numpy as np
import os

from pyworkflow import VERSION_1_1
from pyworkflow.protocol.params import PointerParam, StringParam, FloatParam, BooleanParam, IntParam
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.em.constants import ALIGN_PROJ
from pyworkflow.utils.path import cleanPath, cleanPattern, moveFile, makePath
from pyworkflow.em.protocol import ProtAnalysis3D
from pyworkflow.em.data import SetOfClasses2D, Image, SetOfAverages, SetOfParticles, Class2D
from pyworkflow.em.packages.xmipp3.convert import rowToParticle, setXmippAttributes, xmippToLocation
import pyworkflow.em.metadata as md
from pyworkflow.gui.plotter import Plotter
from itertools import izip

import xmipp
from xmipp3 import ProjMatcher
from pyworkflow.em.packages.xmipp3.convert import rowToAlignment
from scipy.io.arff.arffread import MetaData
# from Crypto.Util.number import str2long

     
FN_MD_INPUT_IMAGES = 'inputImages'
MD_ANGLES = 'mdAngles'
FN_GALLERY_STK = 'galleryStk'
FN_GALLERY_DOC = 'galleryDoc'
FN_MD_ANGLESCONT = 'anglesCont'
FN_SSNR_CLUSTER = 'ssnrCluster'
FN_IMAGES_ITER_1 = 'imagesIter1'
FN_ANGLES_ITER_1 = 'anglesIter1'
        
class XmippProtCompareReprojections(ProtAnalysis3D):
    """Compares a set of classes or averages with the corresponding projections of a reference volume.
    The set of images must have a 3D angular assignment and the protocol computes the residues
    (the difference between the experimental images and the reprojections). The zscore of the mean
    and variance of the residues are computed. Large values of these scores may indicate outliers.
    The protocol also analyze the covariance matrix of the residual and computes the logarithm of
    its determinant [Cherian2013]. The extremes of this score (called zScoreResCov), that is
    values particularly low or high, may indicate outliers."""
    _label = 'compare reprojections'
    _lastUpdateVersion = VERSION_1_1
    
    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', PointerParam, label="Input images", important=True, 
                      pointerClass='SetOfClasses2D, SetOfAverages, SetOfParticles')
        form.addParam('inputVolume', PointerParam, label="Volume to compare images to", 
                      important=True,
                      pointerClass='Volume',
                      help='Volume to be used for class comparison')
        form.addParam('useAssignment', BooleanParam, default=True,
                      label='Use input angular assignment (if available)')
        form.addParam('symmetryGroup', StringParam, default="c1",
                      label='Symmetry group', 
                      help='See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry'
                      ' for a description of the symmetry groups format'
                        'If no symmetry is present, give c1')
        form.addParam('angularSampling', FloatParam, default=5, 
                      expertLevel=LEVEL_ADVANCED,
                      label='Angular sampling rate',
                      help='In degrees.'
                      ' This sampling defines how fine the projection gallery'
                      ' from the volume is explored.')
        form.addParam('maxGrayScale', FloatParam, default=0.95, 
                      expertLevel=LEVEL_ADVANCED,
                      label='Maximum gray scale change',
                      help='Typically the gray scale should be 1,'
                      ' but it might change within the range (1-maxGray,1+maxGray).')
        form.addParam('ssnrGroups', IntParam, default=10, 
                      expertLevel=LEVEL_ADVANCED,
                      label='Number of SSNR groups')
        form.addParallelSection(threads=0, mpi=8)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _createFilenameTemplates(self):
        """ Centralize how files are called """
        myDict = {
                 FN_MD_INPUT_IMAGES: self._getExtraPath('input_imgs.xmd'),
                 MD_ANGLES: self._getExtraPath('angles.xmd'),
                 FN_GALLERY_STK: self._getExtraPath('gallery.stk'),
                 FN_GALLERY_DOC: self._getExtraPath('gallery.doc'),
                 FN_MD_ANGLESCONT: self._getExtraPath('anglesCont.xmd'),
                 FN_SSNR_CLUSTER: self._getExtraPath('ssnrClusters.txt'),
                 FN_IMAGES_ITER_1: self._getExtraPath('images_*iter001_00.xmd'),
                 FN_ANGLES_ITER_1: self._getExtraPath('angles_iter001_00.xmd')
                 }
        self._updateFilenamesDict(myDict)  
    
    def _insertAllSteps(self):
        # Convert input images if necessary
        self._createFilenameTemplates()
        self.imgsFn = self._getFileName(FN_MD_INPUT_IMAGES)
        vol = self.inputVolume.get()
        
        self._insertFunctionStep("convertStep", self.imgsFn)
        imgSet = self.inputSet.get()
        if not self.useAssignment or isinstance(imgSet, SetOfClasses2D) or isinstance(imgSet, SetOfAverages) or (isinstance(imgSet, SetOfParticles) and not imgSet.hasAlignmentProj()):
            anglesFn = MD_ANGLES
            self._insertFunctionStep("projMatchStep", 
                                     self.inputVolume.get().getFileName(), 
                                     self.angularSampling.get(), 
                                     self.symmetryGroup.get(), self.imgsFn,
                                     anglesFn, self.inputVolume.get().getDim()[0])
        else:
            anglesFn=self.imgsFn
        self._insertFunctionStep("produceResiduals", vol.getFileName(), 
                                 anglesFn, vol.getSamplingRate())
        self._insertFunctionStep("evaluateSSNR",self.ssnrGroups.get())
        self._insertFunctionStep("evaluateResiduals")
        self._insertFunctionStep("createOutputStep")

    #--------------------------- STEPS functions ---------------------------------------------------
    def projMatchStep(self, volume, angularSampling, symmetryGroup, 
                      images, fnAngles, Xdim):
        # Generate gallery of projections        
        fnGallery = self._getFileName(FN_GALLERY_STK)
        fnGalleryMd = self._getFileName(FN_GALLERY_DOC)
        if volume.endswith('.mrc'):
            volume+=":mrc"
        
        self.runJob("xmipp_angular_project_library", "-i %s -o %s "
                    "--sampling_rate %f --sym %s --method fourier 1 0.25 bspline"
                    " --compute_neighbors --angular_distance -1 --experimental_images %s"\
                   % (volume, fnGallery, angularSampling, symmetryGroup, images))
    
        # Assign angles
        self.runJob("xmipp_reconstruct_significant", "-i %s " 
                    "--odir %s --initgallery %s --maxShift %d "
                    "--dontReconstruct --useForValidation 0"\
                   % (images, self._getExtraPath(), fnGalleryMd, Xdim/2))
        moveFile(self._getFileName(FN_ANGLES_ITER_1),fnAngles)
        
        cleanPattern(self._getExtraPath('gallery*'))
        cleanPattern(self._getFileName(FN_IMAGES_ITER_1))
    
    def convertStep(self, imgsFn):
        from convert import writeSetOfClasses2D, writeSetOfParticles
        imgSet = self.inputSet.get()
        if isinstance(imgSet, SetOfClasses2D):
            writeSetOfClasses2D(imgSet, self.imgsFn, writeParticles=True)
        else:
            writeSetOfParticles(imgSet, self.imgsFn)
        from pyworkflow.em.convert import ImageHandler
        img = ImageHandler()
        fnVol = self._getTmpPath("volume.vol")
        img.convert(self.inputVolume.get(), fnVol)
        xdim=self.inputVolume.get().getDim()[0]
        if xdim!=self._getDimensions():
            self.runJob("xmipp_image_resize","-i %s --dim %d" 
                        %(fnVol,self._getDimensions()))
    
    def produceResiduals(self, fnVol, fnAngles, Ts):
        if fnVol.endswith(".mrc"):
            fnVol+=":mrc"
        anglesOutFn=self._getExtraPath("anglesCont.stk")
        residualsOutFn=self._getExtraPath("residuals.stk")
        projectionsOutFn=self._getExtraPath("projections.stk")
        xdim=self.inputVolume.get().getDim()[0]
        self.runJob("xmipp_angular_continuous_assign2", 
                    "-i %s -o %s --ref %s --optimizeAngles "
                    "--optimizeGray --max_gray_scale %f --optimizeShift "
                    "--max_shift %d --oresiduals %s --oprojections %s "
                    "--sampling %f --ssnr" %\
                    (fnAngles,anglesOutFn,fnVol,self.maxGrayScale.get(),
                     floor(xdim*0.05),residualsOutFn,projectionsOutFn,Ts))
        fnNewParticles=self._getExtraPath("images.stk")
        if os.path.exists(fnNewParticles):
            cleanPath(fnNewParticles)
    
    def evaluateSSNR(self, ssnrGroups):
        fnCont = self._getFileName(FN_MD_ANGLESCONT)

        mdCont = xmipp.MetaData(fnCont)
        ssnrs = mdCont.getColumnValues(xmipp.MDL_SSNR1D)
        from sklearn.cluster import KMeans
        kmeans = KMeans(n_clusters=ssnrGroups, random_state=0).fit(ssnrs)
        i = 0
        for objId in mdCont:
            mdCont.setValue(xmipp.MDL_SSNR1D_GROUP,int(kmeans.labels_[i]+1),objId)
            i+=1
        mdCont.write(fnCont)
        np.savetxt(self._getFileName(FN_SSNR_CLUSTER), kmeans.cluster_centers_)
#         print(kmeans.cluster_centers_)
    
    def evaluateResiduals(self):
        # Evaluate each image
        fnAutoCorrelations = self._getExtraPath("autocorrelations.xmd")
        stkAutoCorrelations = self._getExtraPath("autocorrelations.stk")
        stkResiduals = self._getExtraPath("residuals.stk")
        anglesOutFn=FN_MD_ANGLESCONT
        self.runJob("xmipp_image_residuals", 
                    " -i %s -o %s --save_metadata_stack %s" 
                    % (stkResiduals, stkAutoCorrelations, fnAutoCorrelations),
                     numberOfMpi=1)
        self.runJob("xmipp_metadata_utilities", 
                    '-i %s --operate rename_column "image imageResidual"' 
                    % fnAutoCorrelations, numberOfMpi=1)
        self.runJob("xmipp_metadata_utilities", 
                    '-i %s --set join %s imageResidual' 
                    % (anglesOutFn, fnAutoCorrelations), numberOfMpi=1)
        cleanPath(fnAutoCorrelations)
    
    def createOutputStep(self):
        self.plotsDir = self._getExtraPath('plots')
        makePath(self.plotsDir)
        filessnr = self._getFileName(FN_SSNR_CLUSTER)
        
        fnImgs = self._getExtraPath('images.stk')
        if os.path.exists(fnImgs):
            cleanPath(fnImgs)

        imgSet = self.inputSet.get()
        imgSet.setSamplingRate(self.inputSet.get().getSamplingRate())
        imgFn = self._getFileName(FN_MD_ANGLESCONT)
        self.newAssignmentPerformed = os.path.exists(MD_ANGLES)
        self.samplingRate = self.inputSet.get().getSamplingRate()

        outputSet = self._createSetOfClasses2D(imgSet)
        self.createPlotsRepresentatives(filessnr, outputSet)

        self.fillingSetOfClasses(imgFn, outputSet)

        self._defineOutputs(outputClasses=outputSet)
        self._defineSourceRelation(self.inputSet, outputSet)

        
    def fillingSetOfClasses(self, mdfn, myclasses):
        mtda = md.MetaData(mdfn)
        for imgRow in md.iterRows(mtda, sortByLabel=xmipp.MDL_SSNR1D_GROUP):
            classid = imgRow.getValue(xmipp.MDL_SSNR1D_GROUP)
            myclass = myclasses[classid]
            part = rowToParticle(imgRow)
            self._processRow(part,imgRow)
            myclass.enableAppend()
            myclass.append(part)
            myclasses.update(myclass)
            
        

    def createPlotsRepresentatives(self, path, myclasses):
        fileSSNR = open(path, 'r')
        count = 1
        for line in fileSSNR: 
            Y = np.fromstring(line, dtype=float, sep=' ')
            X = range(len(Y))
            plotter = self.createSSNRPlot(X, Y, count)
            filename = 'plots/class_%i' % count
            plotter.savefig(self._getExtraPath(filename))
            img = Image(self._getExtraPath(filename+'.png'))
            newclass = Class2D(objId = count)
            newclass.setRepresentative(img)
            myclasses.append(newclass)
            count = count + 1
            
        fileSSNR.close()    

    def _processRow(self, particle, row):
        setXmippAttributes(particle, row,
                           xmipp.MDL_ZSCORE_RESVAR, xmipp.MDL_ZSCORE_RESMEAN, 
                           xmipp.MDL_ZSCORE_RESCOV, xmipp.MDL_IMAGE_ORIGINAL,
                           xmipp.MDL_COST, xmipp.MDL_CONTINUOUS_GRAY_A, 
                           xmipp.MDL_CONTINUOUS_GRAY_B, xmipp.MDL_CONTINUOUS_X,
                            xmipp.MDL_CONTINUOUS_Y, xmipp.MDL_SSNR1D_GROUP)
        
        def __setXmippImage(label):
            attr = '_xmipp_' + xmipp.label2Str(label)
            if not hasattr(particle, attr):
                img = Image()
                setattr(particle, attr, img)
                img.setSamplingRate(particle.getSamplingRate())
            else:
                img = getattr(particle, attr)
            img.setLocation(xmippToLocation(row.getValue(label)))

        __setXmippImage(xmipp.MDL_IMAGE)
        __setXmippImage(xmipp.MDL_IMAGE_REF)
        __setXmippImage(xmipp.MDL_IMAGE_RESIDUAL)
        __setXmippImage(xmipp.MDL_IMAGE_COVARIANCE)



    def createSSNRPlot(self, X, Y, classnumber):
        """ Create a plotter with the cumulative shift per frame. """
        figureSize = (1, 1)
        plotter = Plotter(*figureSize)
        figure = plotter.getFigure()
    
        ax = figure.add_subplot(111)
        ax.grid()
        titlestr = 'SSNR plot %i' %classnumber
        ax.set_title(titlestr)
        ax.set_xlabel(' ')
        ax.set_ylabel('SSNR a.u.')
        ax.plot(0, 0, 'r')
    
        ax.plot(X, Y, color='b')
    
        plotter.tightLayout()
    
        return plotter


    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        summary.append("Images evaluated: %i" % self.inputSet.get().getSize())
        summary.append("Volume: %s" % self.inputVolume.getNameId())
        return summary
    
    def _methods(self):
        methods = []
        if hasattr(self, 'outputClasses'):
            methods.append("We evaluated %i input images %s regarding to volume %s."\
                           %(self.inputSet.get().getSize(), 
                             self.getObjectTag('inputSet'), 
                             self.getObjectTag('inputVolume')) )
        methods.append("The residuals were evaluated according to their mean,"
                       " variance and covariance structure [Cherian2013].")
        return methods
    
    #--------------------------- UTILS functions --------------------------------------------
    def _getDimensions(self):
        imgSet = self.inputSet.get()
        if isinstance(imgSet, SetOfClasses2D):
            xDim = imgSet.getImages().getDim()[0]
        else:
            xDim = imgSet.getDim()[0]
        return xDim
