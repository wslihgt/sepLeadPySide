#!/usr/bin/env python
"""
This is part of the example python scripts from:
    S. Tosi, ``Matplotlib for Python Developers''
    Ed. Packt Publishing
    http://www.packtpub.com/matplotlib-python-development/book

"""
# Python Qt4 bindings for GUI objects
from PySide import QtGui
import matplotlib
# to force the use of PySide, instead of PyQt4 (default with mpl?)
matplotlib.rcParams['backend.qt4'] = 'PySide'
import matplotlib.pyplot as plt
# import the Qt4Agg FigureCanvas object, that binds Figure to
# Qt4Agg backend. It also inherits from QWidget
from matplotlib.backends.backend_qt4agg \
     import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg \
     import NavigationToolbar2QTAgg as NavigationToolbar
# Matplotlib Figure object
from matplotlib.figure import Figure
# grid to make different sizes of axes in a subplot:
from mpl_toolkits.axes_grid1 import make_axes_locatable

class MplCanvas(FigureCanvas):
    """Class to represent the FigureCanvas widget
    
    Attributes
    ----------
    fig : the handler for the matplotlib figure
    ax : the main axe object for the figure
    
    Methods
    -------
    
    Notes
    -----
    taken from S. Tosi, ``Matplotlib for Python Developers''
    Ed. Packt Publishing
    http://www.packtpub.com/matplotlib-python-development/book
    
    """
    def __init__(self):
        """MplCanvas()
        
        creates an object that contains a figure with an axe object in,
        both accessible as attributes.
        
        """
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

class MplCanvas3Axes(FigureCanvas):
    """Class to represent the FigureCanvas widget
    
    Attributes
    ----------
    fig : the figure handler
    ax : the main axe object, in the figure
    ax1 : the axes which are displayed on the left of ax
    ax2 : the axes which are displayed on top of ax
    
    """
    def __init__(self):
        # setup Matplotlib Figure and Axis
        self.fig = Figure()
        # better try with mpl_toolkits.axes_grid1.AxesGrid or ImageGrid?
        self.ax = self.fig.add_subplot(111) # axes for the image
        # axes for the left pane notes
        ## self.ax1 = self.fig.add_subplot(221, sharey=self.ax)
        ## pos = self.ax1.get_position().get_points()
        ## #print pos
        ## pos = [pos[0][0], pos[0][1], 0.1, pos[1][1]-pos[0][1]]
        ## # pos[2] = 40 # changing width
        ## self.ax1.set_position(pos)
        ## # axes for the waveform
        ## self.ax2 = self.fig.add_subplot(224, sharex=self.ax)
        ## pos = self.ax2.get_position().get_points()
        ## #print pos
        ## pos = [pos[0][0], pos[0][1], pos[1][0]-pos[0][0], 0.1]
        ## # pos[3] = 40 # changing height
        ## self.ax2.set_position(pos)
        
        divider = make_axes_locatable(self.ax)
        # help(divider)
        self.ax2 = divider.append_axes("top", size=0.5, pad=0.2,
                                       sharex=self.ax)
        self.ax1 = divider.append_axes("left", size=1., pad=0.3,
                                       sharey=self.ax)
        plt.setp(self.ax.get_yticklabels() + \
                 self.ax2.get_xticklabels() + \
                 self.ax1.get_xticklabels(),
                 visible=False)
        
        # initialization of the canvas
        FigureCanvas.__init__(self, self.fig)
        # we define the widget as expandable
        #FigureCanvas.setSizePolicy(self,
        #                           QtGui.QSizePolicy.Expanding,
        #                           QtGui.QSizePolicy.Expanding)
        # notify the system of updated policy
        FigureCanvas.updateGeometry(self)

class MplWidget(QtGui.QWidget):
    """Widget defined in Qt Designer
    
    Attributes
    ----------
    canvas : the MplCanvas object, which contains the figure, axes, etc.
    vbl : a vertical box layout from QtGui
    navbar : if instanciated with a withNavBar=True, then this attribute
        is the navigation toolbar object (from matplotlib), to allow
        exploration within the axes. By default, withNavBar=False.
    """
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
    """Widget defined in Qt Designer
    
    Attributes
    ----------
    canvas : the MplCanvas object, which contains the figure, axes, etc.
    vbl : a vertical box layout from QtGui
    navbar : if instanciated with a withNavBar=True, then this attribute
        is the navigation toolbar object (from matplotlib), to allow
        exploration within the axes. By default, withNavBar=True.
    """
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

class MplWidget3Axes(QtGui.QWidget):
    """Widget defined in Qt Designer
    
    Attributes
    ----------
    canvas : the MplCanvas3Axes object, which contains the figure, axes, etc.
    vbl : a vertical box layout from QtGui
    navbar : if instanciated with a withNavBar=True, then this attribute
        is the navigation toolbar object (from matplotlib), to allow
        exploration within the axes. By default, withNavBar=True.
        
    """
    def __init__(self, parent = None, withNavBar=True):
        # initialization of Qt MainWindow widget
        QtGui.QWidget.__init__(self, parent)
        # set the canvas to the Matplotlib widget
        self.canvas = MplCanvas3Axes()
        # create a vertical box layout
        self.vbl = QtGui.QVBoxLayout()
        # add mpl widget to the vertical box
        if withNavBar:
            self.navbar = NavigationToolbar(self.canvas, self)
            self.vbl.addWidget(self.navbar)
        
        self.vbl.addWidget(self.canvas)
        # set the layout to the vertical box
        self.setLayout(self.vbl)

class MplMotion3AxesTest(QtGui.QWidget, ):
    
    previous_point = []
    
    def __init__(self, parent=None, withNavBar=True):
        # initialization of Qt MainWindow widget
        QtGui.QWidget.__init__(self, parent)
        # set the canvas to the Matplotlib widget
        self.canvas = MplCanvas3Axes()
        # create a vertical box layout
        self.vbl = QtGui.QVBoxLayout()
        # add mpl widget to the vertical box
        if withNavBar:
            self.navbar = NavigationToolbar(self.canvas, self)
            self.vbl.addWidget(self.navbar)
        
        self.vbl.addWidget(self.canvas)
        # set the layout to the vertical box
        self.setLayout(self.vbl)
        self.canvas.mpl_connect('motion_notify_event', \
                                             self.motion_notify_callback)
        self.canvas.mpl_connect('button_release_event', \
                                             self.mouseReleased)
        self.canvas.mpl_connect('button_press_event', \
                                             self.mousePressed)
        
        self.canvas.ax.imshow([[1,2],[3,4]])
        
    def mouseReleased(self, event):
        x, y = event.xdata, event.ydata
        if event.inaxes:
            self.previous_point = [x, y]
    
    def mousePressed(self, event):
        pass
    
    def motion_notify_callback(self, event):
        """motion_notify_callback
        
        To catch what the mouse is doing on the canvas
        
        Inspired from:
        http://www.mail-archive.com/...
            ...matplotlib-users@lists.sourceforge.net/msg00661.html
        """
        if event.inaxes==self.canvas.ax and \
               self.navbar._active is None: 
            ax = event.inaxes
            # maybe here is the reason why it does not draw everything back
            # after: ax is only the HF0 axes object. 
            self.bg = \
                self.canvas.copy_from_bbox(self.canvas.fig.bbox)
            x, y = event.xdata, event.ydata
            if event.button is None:
                self.previous_point = [x, y]
            elif event.button == 1 and len(self.previous_point) == 0:
                self.previous_point = [x, y]
            elif event.button == 1: # Free Hand Drawing
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
                
                line = plt.Line2D([self.previous_point[0], x],
                                  [self.previous_point[1], y])
                
                line.set_linewidth(1)
                ax.add_line(line)
                self.canvas.restore_region(\
                    self.bg)
                
                self.previous_point = [x, y]
                #print x,y
            
            self.canvas.blit(self.canvas.fig.bbox)
            self.canvas.draw()
    
class SLFilename(QtGui.QLineEdit):
    """ a line edit widget that accepts drag n drops
    
    20110903 : is that not default behaviour anyway?...
    
    Attributes
    ----------
    same as QtGui.QLineEdit

    Methods
    -------
    dragEnterEvent : accepts the dragged object if it is a 'text/uri-list'
    dropEvent : sets the text in the field as the full path to the file which is
        drag and dropped in the field.
    """
    def __init__(self, parent):
        super(SLFilename, self).__init__(parent)
        self.setAcceptDrops(True)
    
    def dragEnterEvent(self, event):
        """dragEnterEvent(self, event)
        
        accepts the dragged object if it is a 'text/uri-list'
        """
        if event.mimeData().hasFormat('text/uri-list'):
            event.accept()
        else:
            event.ignore()
    
    def dropEvent(self, event):
        """dropEvent(self, event)
        
        action when the object is dropped in the widget:
            if it is a file in the system, then the text in the widget is
            set to the full path towards the file.
        """
        self.setText(str(event.mimeData().urls()[0].toLocalFile()))

# An example of use of the above classes:
if __name__=='__main__':
    import sys
    app = QtGui.QApplication(sys.argv)
    ## slmw = MplWidgetNVB()
    #slmw = MplWidget3Axes(withNavBar=True)
    slmw = MplMotion3AxesTest()
    slmw.show()
    sys.exit(app.exec_())
