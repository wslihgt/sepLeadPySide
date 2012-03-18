#!/usr/bin/env python
# This is part of the example python scripts from:
#     S. Tosi, ``Matplotlib for Python Developers''
#     Ed. Packt Publishing
#     http://www.packtpub.com/matplotlib-python-development/book
#

# Python Qt4 bindings for GUI objects
# import PySide as PyQt4
#from PyQt4 import QtGui
from PySide import QtGui

import matplotlib

import numpy as np

# import the Qt4Agg FigureCanvas object, that binds Figure to
# Qt4Agg backend. It also inherits from QWidget
from matplotlib.backends.backend_qt4agg \
     import FigureCanvasQTAgg as FigureCanvas
#from matplotlib.backends.backend_macosx \
#     import FigureCanvasMac as FigureCanvas

from matplotlib.backends.backend_qt4agg \
     import NavigationToolbar2QTAgg as NavigationToolbar

# Matplotlib Figure object
from matplotlib.figure import Figure

class MplCanvas(FigureCanvas):
    """Class to represent the FigureCanvas widget"""
    def __init__(self):
        # setup Matplotlib Figure and Axis
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        
        # initialization of the canvas
        FigureCanvas.__init__(self, self.fig)
        # we define the widget as expandable
        #FigureCanvas.setSizePolicy(self,
        #                           QtGui.QSizePolicy.Expanding,
        #                           QtGui.QSizePolicy.Expanding)
        # notify the system of updated policy
        FigureCanvas.updateGeometry(self)

class MplWidget(QtGui.QWidget):
    """Widget defined in Qt Designer"""
    def __init__(self, parent = None, withNavBar=False):
        # initialization of Qt MainWindow widget
        QtGui.QWidget.__init__(self, parent)
        # set the canvas to the Matplotlib widget
        self.canvas = MplCanvas()
        # create a vertical box layout
        self.vbl = QtGui.QVBoxLayout()
        # add mpl widget to the vertical box
        if withNavBar:
            self.navbar = NavigationToolbar(self.canvas, self)
            self.vbl.addWidget(self.navbar)
        
        self.vbl.addWidget(self.canvas)
        # set the layout to the vertical box
        self.setLayout(self.vbl)

class MplWidgetNVB(QtGui.QWidget):
    """Widget defined in Qt Designer"""
    previous_point = []
    N = 20
    NF0 = 20
    width = 2
    def __init__(self, parent = None, withNavBar=True):
        # initialization of Qt MainWindow widget
        QtGui.QWidget.__init__(self, parent)
        # set the canvas to the Matplotlib widget
        self.canvas = MplCanvas()
        # create a vertical box layout
        self.vbl = QtGui.QVBoxLayout()
        # add mpl widget to the vertical box
        if withNavBar:
            self.navbar = NavigationToolbar(self.canvas, self)
            self.vbl.addWidget(self.navbar)
        
        self.vbl.addWidget(self.canvas)
        # set the layout to the vertical box
        self.setLayout(self.vbl)
        
        self.mask = np.zeros([self.NF0, self.N], dtype=np.bool)
        self.canvas.ax.imshow(self.mask)
        self.canvas.draw()
        
        self.canvas.mpl_connect('motion_notify_event', \
                                self.motion_notify_callback)
        self.canvas.mpl_connect('button_release_event', \
                                self.mouseReleased)
    
    def mouseReleased(self, event):
        x, y = event.xdata, event.ydata
        if event.inaxes:
            self.previous_point = [x, y]
    
    def motion_notify_callback(self, event):
        """motion_notify_callback
        
        To catch what the mouse is doing on the canvas
        
        Inspired from:
        http://www.mail-archive.com/...
            ...matplotlib-users@lists.sourceforge.net/msg00661.html
        """
        if event.inaxes and self.navbar._active is None: 
            ax = event.inaxes
            self.bg = self.canvas.copy_from_bbox(self.canvas.ax.bbox)
            x, y = event.xdata, event.ydata
            NF0 = self.NF0
            x = max(min(x,self.N-1),0)
            y = max(min(y, NF0 - 1), 0)
            if event.button is None:
                self.previous_point = [x, y]
            elif (event.button == 1 or event.button == 3) and \
                     len(self.previous_point) == 0:
                self.previous_point = [x, y]
            elif event.button == 1 or event.button == 3: # Free Hand Drawing
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
                for n in range(2*self.width+1):
                    indexX.extend(indexX0)
                    nuIndexY = indexY0+n-self.width
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
                if event.button == 1:
                    self.mask[indexY,indexX] = True
                elif event.button == 3:
                    self.mask[indexY,indexX] = False
                
                ## line.set_linewidth(1)
                ## ax.add_line(line)
                nuContour = self.canvas.ax.contour(self.mask,
                                       levels=[False,True],
                                       colors='r',
                                       animated=True,
                                       hold=False)
                ## self.canvas.restore_region(self.bg)
                self.canvas.ax.collections = []
                for c in nuContour.collections:
                    self.canvas.ax.add_collection(c)
                
                self.previous_point = [x, y]
                
            # self.canvas.blit(ax.bbox)
            self.canvas.draw()

class SLFilename(QtGui.QLineEdit):
    def __init__(self, parent):
        super(SLFilename, self).__init__(parent)
        self.setAcceptDrops(True)
    
    def dragEnterEvent(self, event):
        if event.mimeData().hasFormat('text/uri-list'):
            event.accept()
        else:
            event.ignore()
    
    def dropEvent(self, event):
        self.setText(str(event.mimeData().urls()[0].toLocalFile()))

if __name__=='__main__':
    import sys
    app = QtGui.QApplication(sys.argv)
    slmw = MplWidgetNVB()
    slmw.show()
    sys.exit(app.exec_())
