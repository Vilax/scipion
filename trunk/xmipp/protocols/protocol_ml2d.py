#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using maximum-likelihood principles
#
#   Author:  Sjors Scheres, January 2008
#  Updated:  J. M. de la Rosa Trevin July 2011
#

from xmipp import MetaData, MDL_ITER, MDL_LL, MDL_REF, MDValueEQ
from protlib_base import XmippProtocol, protocolMain
from config_protocols import protDict
import os
from protlib_utils import runImageJPlugin

class ProtML2D(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.ml2d.name, scriptname, project)
        self.Import = 'from protocol_ml2d import *'
        self.fnIterLogs = self.workingDirPath('ml2d_iter_logs.xmd')
        self.fnIterRefs = self.workingDirPath('ml2d_iter_refs.xmd')
        
    def validate(self):
        errors = []
        if self.DoCorrectAmplitudes and not os.path.exists(self.InCtfDatFile):
            errors.append("Missing '%s' file for correcting amplitudes" % self.InCtfDatFile)
        return errors
    
    def summary(self):
        return []
        lines = [('Input images            ', "%s (%u)" % (self.ImgMd, MetaData(input).size())),
                 ('Reference image', self.RefMd)]
            
        if os.path.exists(self.fnIterLogs):
            md = MetaData(self.fnIterLogs)
            id = md.lastObject()
            iteration = md.getValue(MDL_ITER, id)
            lines.append(('Iteration                   ', str(iteration)))
            LL = md.getValue(MDL_LL, id)
            lines.append(('LogLikelihood          ', str(LL)))
        
        output = ["%s : %s" % (k.ljust(20),  v) for k, v in lines]
        return output
    
    def defineSteps(self):
        progId = "ml"
        if (self.DoMlf):
            progId += "f"  
        
        program = "xmipp_%s_align2d" % progId

        restart = False
        if (restart):
            pass 
            #Not yet implemented
            #params= ' --restart ' + utils_xmipp.composeFileName('ml2d_it', RestartIter,'log')
        else: 
            # Dictionary with boolean options and the cmd options
            booleanDict = {'DoMirror': '--mirror', 'DoNorm': '--norm', 'ZeroOffsets': '--zero_offsets',
                           'FixSigmaNoise': '--fix_sigma_noise', 'FixSigmaOffset': '--fix_sigma_offset',
                           'FixFractions': '--fix_fractions'}
            
            prefix = '%s2d' % progId
            oroot = self.workingDirPath(prefix)
            params = ' -i %s --oroot %s' % (self.ImgMd, oroot)
            # Number of references will be ignored if -ref is passed as expert option
            if self.DoGenerateReferences:
                params += ' --nref %d' % self.NumberOfReferences
            
            if (self.DoFast and not self.DoMlf):
                params += ' --fast'
            if (self.NumberOfThreads > 1  and not self.DoMlf):
                params += ' --thr %i' % self.NumberOfThreads
            if (self.DoMlf):
                if (self.DoCorrectAmplitudes):
                    params += ' --ctfdat %s' % self.InCtfDatFile
                else:
                    params += ' --no_ctf --pixel_size %f' % self.PixelSize
                if (not self.ImagesArePhaseFlipped):
                    params += ' --not_phase_flipped'
                if (self.HighResLimit > 0):
                    params += ' --high %f' % self.HighResLimit
            if self.MaxIters != 100:
                params += " --iter %d" % self.MaxIters
            #Add all boolean options if true
            for k, v in booleanDict.iteritems():
                if self.__dict__[k]:
                    params += " " + v
            
            #Add extra options
            #params += ' ' + self.ExtraParams
        self.insertRunJobStep(program, params)

        self.Db.insertStep('collectResults', WorkingDir=self.WorkingDir, Prefix=prefix,
                           verifyfiles=[self.workingDirPath(f) for f in ['result_images.xmd', 'result_classes.xmd']])

    def visualize(self):
        if os.path.exists(self.fnIterRefs):
            from xmipp import getBlocksInMetaDataFile
            blocks = getBlocksInMetaDataFile(self.fnIterRefs)
            lastBlock = blocks[-1]
            try:
                runImageJPlugin("512m", "XmippBrowser.txt", "-i %s@%s --mode metadata --render" 
                            % (lastBlock, self.fnIterRefs), batchMode=True)
            except Exception, e:
                from protlib_gui_ext import showError
                showError("Error launching java app", str(e))
            launchML2DPlots(self)
        
def collectResults(log, WorkingDir, Prefix):
    oroot = os.path.join(WorkingDir, Prefix)
    mdImgs = MetaData(oroot + '_final_images.xmd')
    outImages = os.path.join(WorkingDir, 'result_images.xmd')
    mdImgs.write('images@' + outImages)
    mdRefs = MetaData(oroot + '_final_refs.xmd')
    outRefs = os.path.join(WorkingDir, 'result_classes.xmd')
    mdRefs.write('classes@' + outRefs)
    mdGroup = MetaData()
    for idx in mdRefs:
        ref = mdRefs.getValue(MDL_REF, idx)
        mdGroup.importObjects( mdImgs, MDValueEQ(MDL_REF, ref))
        mdGroup.writeBlock(outRefs, 'class%06d_images' % ref)

''' Launch some plot for an ML2D protocol run '''
def launchML2DPlots(prot):
    import matplotlib
    matplotlib.use('TkAgg') # do this before importing pylab
    import matplotlib.ticker as ticker
    import matplotlib.gridspec as gridspec
    import matplotlib.pyplot as plt
#    import Tkinter as tk
    
    gs = gridspec.GridSpec(2, 2, height_ratios=[7,4])
    gs.update(left=0.15, right=0.95, hspace=0.25, wspace=0.4)#, top=0.8, bottom=0.2)
    
    #f = plt.figure(figsize=(6, 6), dpi=100)
    #f = plt.figure()
    
    def addIterationsPlot(a, title, xlabel, ylabel, mdlabel, yformat=None, color='blue'):
        #a.get_label().set_fontsize(12)
        a.set_title(title)
        a.set_xlabel(xlabel)
        a.set_ylabel(ylabel)
        
        if yformat:
            formatter = ticker.FormatStrFormatter('%1.2e')
            a.yaxis.set_major_formatter(formatter)
        
        a.xaxis.get_label().set_fontsize(10)
        a.yaxis.get_label().set_fontsize(10)
        labels = a.xaxis.get_ticklabels() + a.yaxis.get_ticklabels()
        for label in labels:
            label.set_fontsize(8)
    
        # Get data from iteration logs metadata
        md = MetaData(prot.fnIterLogs)
        iters = []
        lls = []
        for id in md:
            iters.append(md.getValue(MDL_ITER, id))
            lls.append(md.getValue(mdlabel, id))
        
        a.plot(iters[1:], lls[1:], color=color)
    
    from xmipp import MDL_PMAX, MDL_SIGNALCHANGE
    a = plt.subplot(gs[:-1, :])
    addIterationsPlot(a, 'Log-likelihood (should increase)', 'iterations', 
                      'LL', MDL_LL, '%1.2e')
    a = plt.subplot(gs[1, 0])
    addIterationsPlot(a, 'Probabilities distribution', 'iterations', 
                      'Pmax/Psum',MDL_PMAX, color='green')
    a = plt.subplot(gs[1, 1])
    addIterationsPlot(a, 'Maximum signal change', 'iterations', 
                      'signal change',MDL_SIGNALCHANGE, color='green')    # a tk.DrawingArea
    plt.show()
    
    
