#!/usr/bin/env python
# -*- coding: utf-8 -*-

# copyright (C) 2011 Jean-Louis Durrieu
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

#####################################################
# LICENSE TO USE CODE FROM TUTORIAL EXAMPLE
#####################################################
#
#This file is part of the KDE project
#Copyright (C) 2007 Matthias Kretz <kretz@kde.org>
#adapted 2008 by Thorsten Staerk
#ported to PyQt4 by Christoph Burgmer
#
#Permission to use, copy, modify, and distribute this software
#and its documentation for any purpose and without fee is hereby
#granted, provided that the above copyright notice appear in all
#copies and that both that the copyright notice and this
#permission notice and warranty disclaimer appear in supporting
#documentation, and that the name of the author not be used in
#advertising or publicity pertaining to distribution of the
#software without specific, written prior permission.
#
#The author disclaim all warranties with regard to this
#software, including all implied warranties of merchantability
#and fitness.  In no event shall the author be liable for any
#special, indirect or consequential damages or any damages
#whatsoever resulting from loss of use, data or profits, whether
#in an action of contract, negligence or other tortious action,
#arising out of or in connection with the use or performance of
#this software.


import numpy as np
import sys
import matplotlib
if matplotlib.__version__ < '1.0':
    print Warning("Your Matplotlib version is "+matplotlib.__version__+".\n"+\
                  "This program was tested for Matplotlib version 1.0.1.\n"+\
                  "For < 1.0.1, there might be problems with zooming and\n"+\
                  "selecting the leading source.")
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
from matplotlib.pyplot import setp

from PySide import QtCore, QtGui
from PySide.phonon import Phonon

# that's the GUI from designer file
from separateLeadQtGUI2 import Ui_separateLeadMainWindow

# the class implementing the decomposition/separation algorithm
import SeparateLeadStereo

def db(positiveValue):
    """
    db(positiveValue)
    
    Returns the decibel value of the input positiveValue
    """
    return 10 * np.log10(np.abs(positiveValue))

class SeparateLeadMainWindow(QtGui.QMainWindow, Ui_separateLeadMainWindow):
    """GUI to separate Lead/Accompaniment"""
    
    fileLoaded = False
    previous_point = []
    windowLength = 0.04644 # s
    minF0 = 100 # Hz
    maxF0 = 800 # Hz
    minValHF0 = -100 # dB
    maxValHF0 = 100 # dB
    nbIter = 30
    nbFilters = 10
    nbAccElts = 40
    
    def __init__(self, parent=None):
        """SeparateLeadMainWindow()
        
        Should be instanciated after a QtGui.QApplication(sys.argv)
        call.
        
        """
        super(SeparateLeadMainWindow, self).__init__(parent)
        self.setupUi(self)
        # the sound object:
        self.m_media = None
        self.progressBar=None
        
        # connecting signals to slots:
        QtCore.QObject.connect(self.filenameWidget,
                               QtCore.SIGNAL("clicked()"), self.selectFile)
        QtCore.QObject.connect(self.openFileButton,
                               QtCore.SIGNAL("clicked()"), self.selectFile)
        QtCore.QObject.connect(self.loadFileButton,
                               QtCore.SIGNAL("clicked()"), self.loadFile)
        QtCore.QObject.connect(self.actionOpen,
                               QtCore.SIGNAL('triggered()'),
                               self.selectFile)
        QtCore.QObject.connect(self.actionQuit,
                               QtCore.SIGNAL('triggered()'),
                               QtGui.qApp, QtCore.SLOT("quit()"))
        QtCore.QObject.connect(self.separateButton,
                               QtCore.SIGNAL("clicked()"),
                               self.separate)
        QtCore.QObject.connect(self.separateAutoButton,
                               QtCore.SIGNAL("clicked()"),
                               self.separateAuto)
        
        # for the radio button to normalise display
        QtCore.QObject.connect(self.normaliseNoneButton,
                               QtCore.SIGNAL("clicked()"),
                               self.toggleNormalisation)
        QtCore.QObject.connect(self.normaliseMaxButton,
                               QtCore.SIGNAL("clicked()"),
                               self.toggleNormalisation)
        QtCore.QObject.connect(self.normaliseSumButton,
                               QtCore.SIGNAL("clicked()"),
                               self.toggleNormalisation)
        self.normaliseNoneButton.setEnabled(False)
        self.normaliseMaxButton.setEnabled(False)
        self.normaliseSumButton.setEnabled(False)
        
        # for the radio buttons for the source selection:
        if self.selectLeadButton.isDown():
            self.selectLeadButton.toggle()
        self.selectedPart = 'lead'
        #QtCore.QObject.connect(self.selectAccButton,
        #                       QtCore.SIGNAL("clicked()"),
        #                       self.changeSelectedPart)
        QtCore.QObject.connect(self.selectLeadButton,
                               QtCore.SIGNAL("clicked()"),
                               self.changeSelectedPart)
        QtCore.QObject.connect(self.selectDelButton,
                               QtCore.SIGNAL("clicked()"),
                               self.changeSelectedPart)
        #self.selectAccButton.setEnabled(False)
        self.selectDelButton.setEnabled(False)
        self.selectLeadButton.setEnabled(False)
        
        # connecting for mouse events:
        self.mplHF0Widget.canvas.mpl_connect('motion_notify_event', \
                                             self.motion_notify_callback)
        self.mplHF0Widget.canvas.mpl_connect('button_release_event', \
                                             self.mouseReleased)
        self.mplHF0Widget.canvas.mpl_connect('button_press_event', \
                                             self.mousePressed)
        
        self.mplHF0Widget.canvas.fig.subplots_adjust(right=1.,
                                                     top=0.99,
                                                     left=0.05,
                                                     bottom=0.05)
        
        self.setParamChoices()
    
    def setParamChoices(self):
        # window length choices:
        self.windowLenMS.setRange(0.02322, 0.2322)
        self.windowLenMS.setSingleStep(0.02322)
        self.windowLenMS.setSuffix(" s")
        self.windowLenMS.setDecimals(5)
        self.windowLenMS.setValue(self.windowLength)
        QtCore.QObject.connect(self.windowLenMS,
                               QtCore.SIGNAL("valueChanged(double)"),
                               self.chgValWinLen)
        
        # number of filter components
        self.filtNbSpin.setRange(1, 100)
        self.filtNbSpin.setValue(self.nbFilters)
        QtCore.QObject.connect(self.filtNbSpin,
                               QtCore.SIGNAL("valueChanged(int)"),
                               self.chgValFiltNb)
        
        # number of accompaniment components
        self.accompNbSpin.setRange(1, 400)
        self.accompNbSpin.setValue(self.nbAccElts)
        QtCore.QObject.connect(self.accompNbSpin,
                               QtCore.SIGNAL("valueChanged(int)"),
                               self.chgValAccNb)
        
        # number of iterations
        self.iterSpin.setRange(1, 10000)
        self.iterSpin.setValue(self.nbIter)
        QtCore.QObject.connect(self.iterSpin,
                               QtCore.SIGNAL("valueChanged(int)"),
                               self.chgValNbIter)
        
        # possible values for F0s:
        self.maxF0Spin.setRange(200, 3000)
        self.minF0Spin.setRange(20,500)
        self.maxF0Spin.setSuffix(" Hz")
        self.minF0Spin.setSuffix(" Hz")
        self.maxF0Spin.setValue(self.maxF0)
        self.minF0Spin.setValue(self.minF0)
        QtCore.QObject.connect(self.maxF0Spin,
                               QtCore.SIGNAL("valueChanged(int)"),
                               self.chgValMaxF0)
        QtCore.QObject.connect(self.minF0Spin,
                               QtCore.SIGNAL("valueChanged(int)"),
                               self.chgValMinF0)
        
        # parameters for display control:
        self.maxValEdit.setRange(-1000,1000)
        self.minValEdit.setRange(-1000,1000)
        self.maxValEdit.setSuffix(" dB")
        self.minValEdit.setSuffix(" dB")
        self.maxValEdit.setValue(self.maxValHF0)
        self.minValEdit.setValue(self.minValHF0)
        self.maxValEdit.setEnabled(False)
        self.minValEdit.setEnabled(False)
        QtCore.QObject.connect(self.maxValEdit,
                               QtCore.SIGNAL("valueChanged(int)"),
                               self.chgValMaxVal)
        QtCore.QObject.connect(self.minValEdit,
                               QtCore.SIGNAL("valueChanged(int)"),
                               self.chgValMinVal)
        
        # for the height selection of the "pencil"
        self.heightSelectSpin.setRange(0.,20.)
        # self.heightSelectSpin.setSuffix(" semitone")
        self.heightSelectSpin.setValue(1.)
        self.heightSelectSpin.setEnabled(False)
        QtCore.QObject.connect(self.heightSelectSpin,
                               QtCore.SIGNAL("valueChanged(double)"),
                               self.chgValWidthF0)
    
    def chgValWinLen(self, value):
        self.windowLength = value
        # print self.windowLength
    
    def chgValFiltNb(self, value):
        self.nbFilters = value
    
    def chgValAccNb(self, value):
        self.nbAccElts = value
    
    def chgValNbIter(self, value):
        self.nbIter = value
    
    def chgValMinVal(self, value):
        self.minValHF0 = value
        self.updateClimmplHF0()
    
    def chgValMaxVal(self, value):
        self.maxValHF0 = value
        self.updateClimmplHF0()
    
    def chgValMinF0(self, value):
        if value < self.maxF0:
            self.minF0 = value
        else:
            self.minF0Spin.setValue(self.minF0)
    
    def chgValMaxF0(self, value):
        if value > self.minF0:
            self.maxF0 = value
        else:
            self.maxF0Spin.setValue(self.maxF0)
    
    def chgValWidthF0(self, value):
        if self.fileLoaded:
            self.widthF0 = np.int(value*\
                           self.separateInstance.SIMMParams['stepNotes'])
    
    def selectFile(self):
        """opens a file select dialog"""
        # open the dialog and get the selected file
        file = QtGui.QFileDialog.getOpenFileName()
        # if a file is selected
        if file:
            # update the lineEdit widget text with the selected filename
            self.filenameWidget.setText(file[0])
    
    def loadFile(self):
        """loads the file whose name is in the filenameWidget"""
        self.statusbar.showMessage("Loading file...")
        # the sound object:
        # (reinitialize sound object, so that if changing file,
        #  also changing sound)
        self.m_media = None
        #self.windowLenMS.setEnabled(False)
        #self.maxF0Spin.setEnabled(False)
        #self.minF0Spin.setEnabled(False)
        outputDirSuffix = str(self.outputDirSuffixEdit.text())
        self.separateInstance = SeparateLeadStereo.SeparateLeadProcess(\
                                    self.filenameWidget.text(),
                                    imageCanvas=self.mplHF0Widget.canvas,
                                    wavCanvas=self.mplHF0Widget.canvas,
                                    progressBar=self.progressBar,
                                    nbIter=self.nbIter,
                                    minF0=self.minF0, 
                                    maxF0=self.maxF0,
                                    windowSize=self.windowLength,
                                    verbose=True,
                                    outputDirSuffix=outputDirSuffix,
                                    K_numFilters=self.nbFilters,
                                    numCompAccomp=self.nbAccElts)
        self.statusbar.showMessage("Estimating SIMM parameters...")
        self.separateInstance.estimSIMMParams()
        self.mplHF0Widget.canvas.draw()
        self.statusbar.showMessage("Estimating melody line...")
        self.separateInstance.runViterbi()
        self.mplHF0Widget.canvas.draw()
        self.separateInstance.initiateHF0WithIndexBestPath()
        self.statusbar.showMessage("Ready! Please select area to separate, "+\
                                   "or run Separate (Auto) for automatic "+\
                                   "melody.")
        self.fileLoaded = True
        imgHF0clim = self.mplHF0Widget.canvas.ax.get_images()[0].get_clim()
        self.maxValHF0 = imgHF0clim[1]
        self.minValHF0 = imgHF0clim[0]
        self.maxValEdit.setValue(self.maxValHF0)
        self.minValEdit.setValue(self.minValHF0)
        self.maxValEdit.setEnabled(True)
        self.minValEdit.setEnabled(True)
        self.heightSelectSpin.setValue(1.)
        self.heightSelectSpin.setEnabled(True)
        
        # creating masks for lead and accompaniment
        F = self.separateInstance.SIMMParams['NF0']
        N = self.separateInstance.N
        # using the estimated mask instead:
        ##self.leadMask = np.zeros([F,N], dtype=np.bool)
        self.leadMask = (self.separateInstance.SIMMParams['HF00']>0)
        self.mplHF0Widget.canvas.ax.contour(self.leadMask,
                   levels=[False,True],
                   colors='r',
                   animated=False)
        self.mplHF0Widget.canvas.draw()
        # dropping the accompaniment mask
        ##self.accMask = np.zeros([F,N], dtype=np.bool) # TODO: not used?
        
        # for the display, set the width of annotation
        self.widthF0 = self.separateInstance.SIMMParams['stepNotes']
        
        self.normaliseNoneButton.setEnabled(True)
        self.normaliseMaxButton.setEnabled(True)
        self.normaliseSumButton.setEnabled(True)
        self.selectDelButton.setEnabled(True)
        self.selectLeadButton.setEnabled(True)
        
        self.drawMusicalLines()
        setp(self.mplHF0Widget.canvas.ax.get_yticklabels() + \
             self.mplHF0Widget.canvas.ax2.get_xticklabels() + \
             self.mplHF0Widget.canvas.ax1.get_xticklabels(),
             visible=False)
    
    def separate(self):
        """Runs the separation script"""
        outputDirSuffix = str(self.outputDirSuffixEdit.text())
        self.separateInstance.setOutputFileNames(outputDirSuffix=\
                                                 outputDirSuffix)
        self.separateInstance.SIMMParams['HF0'] *= self.leadMask
        self.separateInstance.automaticMelodyAndSeparation()
    
    def separateAuto(self):
        """Runs the separation script"""
        outputDirSuffix = str(self.outputDirSuffixEdit.text())
        self.separateInstance.setOutputFileNames(outputDirSuffix=\
                                                 outputDirSuffix)
        self.separateInstance.automaticMelodyAndSeparation()
    
    def changeSelectedPart(self):
        # the following code should "unset" the navigation toolbar options:
        if self.mplHF0Widget.navbar._active != 'PAN':
            self.mplHF0Widget.navbar.pan()
        self.mplHF0Widget.navbar.pan()
        
        checkedButton = self.buttonGroup.checkedButton()
        if checkedButton is self.selectDelButton:
            self.selectedPart = None
        elif checkedButton is self.selectLeadButton:
            self.selectedPart = 'lead'
        ## print self.selectedPart
    
    def toggleNormalisation(self):
        checkedButton = self.buttonGroup_2.checkedButton()
        HF0 = self.separateInstance.SIMMParams['HF0']
        NF0 = self.separateInstance.SIMMParams['NF0']
        if checkedButton is self.normaliseNoneButton:
            self.mplHF0Widget.canvas.ax.get_images(\
                )[0].set_data(db(HF0))
            self.mplHF0Widget.canvas.ax.get_images(\
                )[0].set_clim(db(HF0).max()-50, db(HF0).max())
        elif checkedButton is self.normaliseMaxButton:
            self.mplHF0Widget.canvas.ax.get_images(\
                )[0].set_data(db(HF0) -
                              np.outer(np.ones(NF0),
                                       np.max(db(HF0),axis=0)))
            self.mplHF0Widget.canvas.ax.get_images(\
                )[0].set_clim(-50, 0)
        elif checkedButton is self.normaliseSumButton:
            self.mplHF0Widget.canvas.ax.get_images(\
                )[0].set_data(db(HF0) - 
                              db(np.outer(np.ones(NF0),
                                          np.sum(HF0,axis=0))))
            self.mplHF0Widget.canvas.ax.get_images(\
                )[0].set_clim(-50, 0)
            
        self.mplHF0Widget.canvas.draw()
        imgHF0clim = self.mplHF0Widget.canvas.ax.get_images()[0].get_clim()
        self.maxValHF0 = imgHF0clim[1]
        self.minValHF0 = imgHF0clim[0]
        self.maxValEdit.setValue(self.maxValHF0)
        self.minValEdit.setValue(self.minValHF0)
    
    def mouseReleased(self, event):
        x, y = event.xdata, event.ydata
        if event.inaxes:
            self.previous_point = [x, y]
    
    def mousePressed(self, event):            
        if event.inaxes in (self.mplHF0Widget.canvas.ax,
                            self.mplHF0Widget.canvas.ax2) and \
               event.button == 3 and self.fileLoaded:
            self.playStop()
    
    def motion_notify_callback(self, event):
        """motion_notify_callback
        
        To catch what the mouse is doing on the canvas
        
        Inspired from:
        http://www.mail-archive.com/...
            ...matplotlib-users@lists.sourceforge.net/msg00661.html
        """
        if event.inaxes==self.mplHF0Widget.canvas.ax and \
               self.mplHF0Widget.navbar._active is None and \
               self.fileLoaded: 
            ax = event.inaxes
            # maybe here is the reason why it does not draw everything back
            # after: ax is only the HF0 axes object. 
            self.mplHF0bg = \
                self.mplHF0Widget.canvas.copy_from_bbox(self.mplHF0Widget.canvas.fig.bbox)
            x, y = event.xdata, event.ydata
            NF0 = self.separateInstance.SIMMParams['NF0']
            x = max(min(x,self.separateInstance.N-1),0)
            y = max(min(y, NF0 - 1), 0)
            if event.button is None:
                self.previous_point = [x, y]
            elif event.button == 1 and len(self.previous_point) == 0:
                self.previous_point = [x, y]
            elif event.button == 1 and \
                   hasattr(self, 'leadMask'): # Free Hand Drawing
                if self.previous_point[0] < x:
                    x0 = self.previous_point[0]
                    x1 = x
                    xstep = 1
                else:
                    x1 = self.previous_point[0]
                    x0 = x
                    xstep = -1
                if self.previous_point[1] < y:
                    y0 = self.previous_point[1]
                    y1 = y
                    ystep = 1
                else:
                    y1 = self.previous_point[1]
                    y0 = y
                    ystep = -1
                
                ## line = Line2D([self.previous_point[0], x],
                ##               [self.previous_point[1], y])
                indexX = np.arange(x0, x1+1, dtype=np.int32)[::xstep]
                indexY = np.arange(y0, y1+1, dtype=np.int32)[::ystep]
                if indexX.size > indexY.size:
                    indexY = np.linspace(y0, y1+1, num=indexX.size)[::ystep]
                    indexY = np.int32(indexY)
                else:
                    indexX = np.linspace(x0, x1+1, num=indexY.size)[::xstep]
                    indexX = np.int32(indexX)
                
                indexX0 = indexX.tolist()
                indexX = indexX.tolist()
                indexY0 = np.copy(indexY)
                indexY = indexY.tolist()
                for n in range(2*self.widthF0+1):
                    indexX.extend(indexX0)
                    nuIndexY = indexY0+n-self.widthF0
                    nuIndexY = np.maximum(nuIndexY, 0)
                    nuIndexY = np.minimum(nuIndexY, NF0 - 1)
                    indexY.extend(nuIndexY.tolist())
                
                ## indexXGraph = []
                ## indexYGraph = []
                ## indexXGraph.extend(indexX0)
                ## indexYGraph.extend((indexY0-self.widthF0).tolist())
                ## indexXGraph.extend(indexX0)
                ## indexYGraph.extend((indexY0+self.widthF0).tolist())
                ## line = Line2D(indexXGraph, indexYGraph)
                if self.selectedPart == 'lead':
                    self.leadMask[indexY,indexX] = True
                    ## line.set_color('k')
                else:
                    self.leadMask[indexY,indexX] = False
                    ## self.accMask[indexY,indexX] = False
                    ## line.set_color('w')
                
                ## line.set_linewidth(1)
                ## ax.add_line(line)
                self.mplHF0Widget.canvas.restore_region(\
                    self.mplHF0bg)
                nuContour = ax.contour(self.leadMask,
                                       levels=[False,True],
                                       colors='r',
                                       animated=True)
                ax.collections = []
                for c in nuContour.collections:
                    ax.add_collection(c)
                
                self.previous_point = [x, y]
            
            self.mplHF0Widget.canvas.blit(self.mplHF0Widget.canvas.fig.bbox)
            self.mplHF0Widget.canvas.draw()
    
    def updateClim(self, canvas, minVal, maxVal):
        canvas.ax.get_images()[0].set_clim(minVal, maxVal)
        canvas.draw()
    
    def updateClimmplHF0(self):
        self.updateClim(self.mplHF0Widget.canvas,
                        self.minValHF0,
                        self.maxValHF0)
    
    def getX0X1inSecondFromView(self):
        #x0, x1, y0, y1 = self.mplHF0Widget.canvas.ax.get_images()[0].get_clip_box()
        x0, x1 = self.mplHF0Widget.canvas.ax.get_xlim()
        
        x0InSec = x0 * self.separateInstance.stftParams['hopsize'] * \
                  1. / self.separateInstance.fs
        x1InSec = x1 * self.separateInstance.stftParams['hopsize'] * \
                  1. / self.separateInstance.fs
        
        return x0InSec, x1InSec
    
    def play(self):
        """
        from http://techbase.kde.org/Development/Tutorials/Phonon/Introduction/Python
        """
        self.delayedInit()
        #self.m_media.setCurrentSource(self.m_model.filePath(index))
        self.m_media.setCurrentSource(
            Phonon.MediaSource(self.separateInstance.files['inputAudioFilename']))
        self.m_media.play()
    
    def playStop(self, msecToEnd=0):
        """
        if stopped, plays the sound, if played, stops the music.
        """
        
        self.delayedInit()
        ## self.m_media.setCurrentSource(
        ##     Phonon.MediaSource(\
        ##         self.separateInstance.files['inputAudioFilename']))
        
        ## print self.m_media.isSeekable() #DEBUG
        ## print self.m_media.currentTime() #DEBUG
        
        if self.m_media.state() in (Phonon.StoppedState,
                                    Phonon.PausedState):
            x0InSec, x1InSec = self.getX0X1inSecondFromView()
            ## print x0InSec, x1InSec #DEBUG!
            self.m_media.seek(np.maximum(0,np.int64(x0InSec*1e3)))
            # cheating: the prefinish mark is used to stop the audio playback
            self.m_media.setPrefinishMark(np.maximum(\
                np.int64(self.m_media.totalTime() - x1InSec * 1e3),0))
            ## print np.maximum(0,np.int64(x0InSec*1e3)),\
            ##       np.maximum(np.int64(self.m_media.totalTime()-\
            ##                           x1InSec*1e3),0) #DEBUG!
            ## print self.m_media.currentTime() #DEBUG
            self.m_media.play()
            ## print self.m_media.currentTime() #DEBUG
            
        elif self.m_media.state() in (Phonon.PlayingState,):
            self.m_media.stop()
    
    def delayedInit(self):
        """
        from http://techbase.kde.org/Development/Tutorials/Phonon/Introduction/Python
        """
        if not self.m_media:
            self.m_media = Phonon.MediaObject(self)
            audioOutput = Phonon.AudioOutput(Phonon.MusicCategory, self)
            Phonon.createPath(self.m_media, audioOutput)
            self.m_media.setCurrentSource(
                Phonon.MediaSource(\
                    self.separateInstance.files['inputAudioFilename']))
            self.m_media.prefinishMarkReached.connect(self.playStop)
            ##QtCore.QObject.connect(self.m_media,
            ##                      QtCore.SIGNAL("prefinishMarkReached()"),
            ##                      self.playStop)
    
    def drawMusicalLines(self):
        """
        draws musical score staff on the left side axes.
        
        TODO: improve the look of this musical staff...
        """
        minF0 = self.separateInstance.SIMMParams['minF0']
        stepNotes = self.separateInstance.SIMMParams['stepNotes']
        A4 = 440 # Hz
        A3 = 220 # Hz
        laNotesCenters = stepNotes * np.array([-12, 0, 12, 24]) + \
                         12 * stepNotes * np.log2(A3/np.double(minF0))
        laNotesLabels = ['A2', 'A3', 'A4', 'A5']
        longLinesCters = stepNotes * np.array([-14, -10, -7, -4, 0,
                                   7, 10, 14, 17, 20]) + \
                         12 * stepNotes * np.log2(A3/np.double(minF0))
        smallLinesCtrs = stepNotes * np.array([3, 24]) + \
                         12 * stepNotes * np.log2(A3/np.double(minF0))
        # drawing long lines:
        for y in longLinesCters:
            self.mplHF0Widget.canvas.ax1.plot([0, 1], [y, y], color='k')
        for y in smallLinesCtrs:
            self.mplHF0Widget.canvas.ax1.plot([0.55, 0.95], [y, y], color='k')
        for i, y in enumerate(laNotesCenters):
            circ = Ellipse([0.75, y],
                           height=stepNotes, width=.15, angle=-.5,
                           color='b')
            self.mplHF0Widget.canvas.ax1.add_patch(circ)
            self.mplHF0Widget.canvas.ax1.text(0.85, y, laNotesLabels[i])
        
        self.mplHF0Widget.canvas.ax1.set_xlim(-0.2, 1.2)
        self.mplHF0Widget.canvas.draw()

    def reset(self):
        """TODO: reset function, to reinitialize everything (especially HF0
        pane)
        """
        pass

def main():
    # create the GUI application
    app = QtGui.QApplication(sys.argv)
    app.setApplicationName("GUI for Audio Separation")
    slmw = SeparateLeadMainWindow()
    slmw.show()
    sys.exit(app.exec_())

if __name__=='__main__':
    main()
