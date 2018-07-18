#!/usr/bin/python

"""
Visualising bayesian localisation clustering data.

Michael Hirsch Feb 2017
"""
import sys, os, argparse
from PyQt4.QtCore import SIGNAL, Qt, SLOT, QObject, QSettings, QString, QSize, QEvent
from PyQt4.QtCore import QVariant
from PyQt4.QtGui import QIcon, QSizePolicy, QPixmap, QAction, QMessageBox
from PyQt4 import QtGui


import matplotlib
matplotlib.use('QT4AGG')

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from PIL import Image
from io import BytesIO
from matplotlib import image
from functools import partial

from shapely.geometry.polygon import Polygon
from shapely.geometry import Point, MultiPoint
from shapely.ops import unary_union
from shapely import affinity

import logging

import numpy as np
from scipy.spatial import KDTree
import csv, json

import qrc_resources

COLOURS = ['red', 'green']
NANOMICRON = 1000000.0
DILATE = 750
THRESHOLD = 0.025 # density threshold for segmentation

logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
logFormatter = logging.Formatter('%(levelname)s: %(message)s')
log = logging.getLogger()

def rotate(list_of_pairs):
    ''' rotates the pts vector so that the lower left corner comes first'''
    idx, pt0 = min(enumerate(list_of_pairs), key = lambda p : p[1][1])
    return list_of_pairs[idx:] + list_of_pairs[0:idx]

def angle(list_of_pairs, unit='rad'):
    ''' anti-clockwise rotation around the lower left corner '''
    pt0 = list_of_pairs[0]
    pt1 = list_of_pairs[1]
    q = PointVector(1, 0)
    base = PointVector(*pt1) - PointVector(*pt0)
    if unit == 'rad':
        return np.arccos((q*base)/base.norm())
    elif unit == 'deg':
        return np.arccos((q*base)/base.norm())*180.0/np.pi
    else:
        raise Exception("Unknown unit for angle (known: 'rad' or 'deg')")

class Shape(Polygon):
    ''' a shapely polygon + a matplotlib polygon '''
    POLYGON, RECTANGLE = range(2)
    def __init__(self, coords, shape_type):
        super(Shape, self).__init__(coords)
        self._patch = None
        self._shape_type = shape_type
        self._highlighted = False

    def __del__(self):
        try:
            if not self._patch is None and self._patch.is_figure_set():
                self._patch.remove()
        except:
            print('no patch to remove or something went wrong')
        super(Shape, self).__del__()

    def isHighlighted(self):
        return self._highlighted

    def update(self, **kwargs):
        #assert(not self._patch is None)
        ops = {
            'fc'    : self._patch.set_fc,
            'alpha' : self._patch.set_alpha,
        }
        for key, val in kwargs.items():
            assert(key) in ops
            ops[key](val)

    def highlight(self):
        self.update(fc='b', alpha=0.2)
        self._highlighted = True

    def shade(self):
        self.update(fc='k', alpha=0.2)
        self._highlighted = False

    def draw(self, artist):
        if self._patch is None:
            self._patch = matplotlib.patches.Polygon(self.exterior)
        artist.add_patch(self._patch)

    def shape_type(self):
        return self._shape_type

    def remove(self):
        ''' removes the patch from the figure '''
        self._patch.remove()

    def angle_origin(self, unit='rad'):
        ''' (angle, (x0, y0), width, heigth) for the lowest, most left corner of the shape'''
        exterior = self.exterior # LineRing
        coords = list(exterior.coords) # list of pairs
        if not exterior.is_ccw:
            coords.reverse()
        p0 = coords.pop() # last pair
        assert(p0 == coords[0]) # is doubled
        coords = rotate(coords) # puts the lowest, most left point first
        return (
            angle(coords, unit),
            coords[0],
            (PointVector(*coords[1]) - PointVector(*coords[0])).norm(),
            (PointVector(*coords[-1]) - PointVector(*coords[0])).norm()
            )

class ShapeManger(QObject):
    def __init__(self, parent=None):
        super(ShapeManger, self).__init__(parent)
        self._shapes = []

    def __len__(self):
        return len(self._shapes)

    def __iter__(self):
        return self._shapes.__iter__()

    def pop(self, idx=-1):
        ''' remove a shape '''
        shape = self._shapes.pop(idx)
        self.emit(SIGNAL("storedShapesChanged"))
        return shape
    def push(self, shape):
        ''' add a new shape '''
        self._shapes.append(shape)
        self.emit(SIGNAL("storedShapesChanged"))
    def clear(self):
        ''' remove all shapes '''
        something_to_remove = len(self) > 0
        self._shapes = []
        if something_to_remove:
            self.emit(SIGNAL("storedShapesChanged"))
    def unselectAll(self):
        ''' un-highlight all shapes '''
        for shape in self._shapes:
            shape.shade()
    def select(self, idx):
        ''' highlight the shape with the given index '''
        self.unselectAll()
        self._shapes[idx].highlight()

    def highlighted(self):
        ''' list of highlighted shapes '''
        return [ shape for shape in  self._shapes if shape.isHighlighted() ]

    def popHighlighted(self):
        ''' removes all highlighted shapes '''
        remove_this = self.highlighted()
        something_to_remove = len(remove_this) > 0
        for shape in remove_this:
            self._shapes.remove(shape)
        if something_to_remove:
            self.emit(SIGNAL("storedShapesChanged"))
    def remove(self, shape):
        if shape in self._shapes:
            self._shapes.remove(shape)
            self.emit(SIGNAL("storedShapesChanged"))

class MyKdTree(object):
    def __init__(self, dataSet):
        self._data = np.array(zip(dataSet.x(), dataSet.y()))
        self._tree = KDTree(self._data)
        self._area = self._defaultArea(1.0)

    def _defaultArea(self, default = 1.0):
        if len(self._data) == 0:
            return default
        xmin, ymin = self._data.min(axis=0)
        xmax, ymax = self._data.max(axis=0)
        return (xmax-xmin) * (ymax-ymin)

    def __len__(self):
        return len(self._data)

    def area(self):
        return self._area

    def ripleyK(self, r):
        return float(sum([ len(self.ballPoint(x, r)) for x in self._data ]))/(len(self)**2)*self.area()

    def ripleyL(self, r):
        return np.sqrt(self.ripleyK(r)/np.pi)

    def ballPoint(self, x, r):
        ''' indices of all points with radius r of x '''
        return self._tree.query_ball_point(x, r)

    def point(self, i):
        if isinstance(i, int):
            return self._data[i]
        if isinstance(i, (list, tuple, set)):
            return [ self._data[j] for j in i ]

class GaussianTypeKernel(object):
    def __init__(self, sigma):
        self._fact = 1.0/(sigma*np.sqrt(2*np.pi))
        self._two_sigma_square = 2.0*sigma*sigma

    def _square_dist(self, x, y):
        return (x[0]-y[0])**2 + (x[1]-y[1])**2

    def __call__(self, x, y):
        return self._fact*np.exp(-self._square_dist(x, y)/self._two_sigma_square)

class PointVector(object):
    ''' point as vector '''
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def __add__(self, q):
        return PointVector(self.x+q.x, self.y+q.y)
    def __sub__(self, q):
        return PointVector(self.x-q.x, self.y-q.y)
    def __mul__(self, q):
        ''' scalar product or scaling, depending on the type of q '''
        if isinstance(q, PointVector):
            return float(self.x*q.x + self.y*q.y)
        elif isinstance(q, (int, float)):
            return PointVector(self.x*q, self.y*q)
        else:
            raise Exception("invalid operand")

    def __truediv__(self, n):
        return self.__div__(n)
    def __div__(self, n):
        if isinstance(n, (int, float)):
            return PointVector(1.0*self.x/n, 1.0*self.y/n)
        else:
            raise Exception("invalid operand")
    def norm(self):
        ''' Euclidean  norm '''
        return np.sqrt(self*self)

    def __str__(self):
        return "[{0}, {1}]".format(self.x, self.y)

def writeShapeConfig(folder, w, h):
    config = BlcConfig()
    config['xlim'] = (0, int(np.ceil(w)))
    config['ylim'] = (0, int(np.ceil(h)))
    config.write(os.path.join(folder, 'config.txt'))


def get_data_fname(folder0, colour):
    assert(colour in ['red', 'green'])
    folder = folder0 + colour
    fname = os.path.join(folder, '__data.txt')
    if os.path.exists(fname):
        return fname
    fname = os.path.join(folder, "1", 'data.txt')
    if os.path.exists(fname):
        return fname
    return None

def sortPoints(p1, p2):
    ''' return the most left and then most lower point first '''
    if p1.x < p2.x or (p1.x == p2.x and p1.y < p2.y):
        return p1, p2
    else:
        return p2, p1

class Rectangler(object):
    ''' A rectangle created by three points '''
    def __init__(self, p1, p2, q):
        self._p0, self._p3 = sortPoints(p1, p2)
        self._q = q
        self._pts = []
        self._calc()

    def _calc(self):
        ''' calculate  a rectangle based on 3 points, p1, p2 and q:
        Side s1 of the rectangle is p1->p2, its parallel side s3 contains q,
        the sides s2 and s4 contain p1 and p2 respective and are perpenticular to s1.
        '''
        p0, p3, q = self._p0, self._p3, self._q
        v30 = p3-p0
        n30 = v30/v30.norm() # normalised P0->P3
        v0q = p0-q # vector P0->Q
        lpq = v0q*n30 # length of Q projected on P0->P3
        npq = (q-p0) - n30*(-lpq) # vector perpenticular to P0->P3; later P0->P1
        p1 = p0 + npq
        p2 = p1 + n30*v30.norm()
        assert((p3-p2)*v30 < 0.0001)  # (P2->P3) * (P0->P3) should be 0
        self._pts= [(p0.x, p0.y), (p1.x, p1.y), (p2.x, p2.y), (p3.x, p3.y)]
        self._rotate()
        self._orient()

    def _orient(self):
        ''' orients the points anti-clockwise assumes pts have been rotated '''
        x0 = self.points()[0][0]
        xp = self.points()[-1][0]
        xs = self.points()[1][0]
        if xp < x0 or x0 < xs:
            return
        self._pts.reverse()
        p0 = self._pts.pop()
        self._pts.insert(0, p0)

    def _rotate(self):
        ''' rotates the pts vector so that the lower left corner comes first'''
        idx, pt0 = min(enumerate(self._pts), key = lambda p : p[1][1])
        self._pts = self._pts[idx:] + self._pts[0:idx]

    def angle(self, unit='rad'):
        ''' anti-clockwise rotation around the lower left corner '''
        pt0 = self.points()[0]
        pt1 = self.points()[1]
        q = PointVector(1, 0)
        base = PointVector(*pt1) - PointVector(*pt0)
        if unit == 'rad':
            return np.arccos((q*base)/base.norm())
        elif unit == 'deg':
            return np.arccos((q*base)/base.norm())*180.0/np.pi
        else:
            raise Exception("Unknown unit for angle (known: 'rad' or 'deg')")

    def points(self):
        ''' return the four points of the rectangle as a list of pairs '''
        return self._pts

    def width(self):
        pt0 = PointVector(*self.points()[0])
        pt1 = PointVector(*self.points()[1])
        return (pt1-pt0).norm()

    def height(self):
        pt0 = PointVector(*self.points()[0])
        pt1 = PointVector(*self.points()[3])
        return (pt1-pt0).norm()

class BlcOps(QObject):
    ''' container for operations on data '''
    def __init__(self, parent=None):
        super(BlcOps, self).__init__(parent)
        self._operations = {
            'cut' : self._cut,
            'crop' : self._crop,
            'shift' : self._shift,
            'rotate' : self._rotate,
            'count' : self._count,
        }

    def _cut(self, data, polygon):
        if data is None:
            return []
        return [ item for item in data if not polygon.contains(Point(*item[0:2])) ]

    def _crop(self, data, polygon):
        if data is None:
            return []
        return [ item for item in data if polygon.contains(Point(*item[0:2])) ]

    def _count(self, data, shape):
        if data is None:
            return 0
        return len([ 1 for item in data if shape.contains(Point(*item[0:2])) ])

    def _shift(self, data, point):
        if data is None:
            return []
        x0, y0 = point
        return [ [x-x0, y-y0, sd ] for x, y, sd in data ]

    def _rotate(self, data, angle):
        ''' angle is in "rad" '''
        asin = np.sin(angle)
        acos = np.cos(angle)
        rotx = lambda x, y : x*acos - y*asin
        roty = lambda x, y : x*asin + y*acos
        return [ [rotx(x, y), roty(x, y), sd ] for x, y, sd in data ]

    def __getitem__(self, key):
        return self._operations[key]

class BlcConfig(QObject):
    ''' BLC config file handler '''
    def __init__(self, parent=None):
        super(BlcConfig, self).__init__(parent)
        self._data = {
            'model' : 'Gaussian(prec)',
            'alpha' : 20,
            'pbackground' : .5,
            'rseq' : (5, 150, 5),
            'thseq' : (5, 350, 5),
            'histbins' : (1, 100),
            'histvalues' : (1, 1),
            'xlim' : (0, 25000),
            'ylim' : (0, 25000),
        }

    def _writeRow(self, key, val):
        if isinstance(val, (str, int, float)):
            return "{0}={1}\n".format(key, val)
        else:
            tmp = [ str(x) for x in val ]
            return "{0}={1}\n".format(key, ','.join(tmp))

    def write(self, fname):
        with open(fname, 'w') as fh:
            for key, val in self._data.items():
                fh.write(self._writeRow(key, val))

    def __getitem__(self, key):
        return self._data[key]

    def __setitem__(self, key, val):
        assert(key in self._data)
        assert(isinstance(val, (str, int, float, tuple, list)))
        self._data[key] = val

class BlcData(QObject):
    ''' a single [x,y,sd] data set '''
    def __init__(self, parent=None):
        super(BlcData, self).__init__(parent)
        self._data = []
        self._orig_data = []
        self._header = None
        self._have_data = False
        self._fname = None
        self._targetFile = None

    @staticmethod
    def init(header, data):
        blcData = BlcData()
        blcData._data = data
        blcData._have_data = True
        blcData._header = header
        return blcData

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        return self._data.__iter__()

    def __getitem__(self, index):
        return self._data[index]

    def points(self, index_list):
        return [ [x, y] for x, y, sd in  [ self._data[i] for i in index_list ] ]


    def header(self):
        return self._header

    def loadData(self, fname):
        with open(fname) as fh:
            reader = csv.reader(fh)
            self._header = reader.next()
            if not self.testHeader():
                return False
            self._data = [ [float(x), float(y), float(sd)] for x, y, sd in reader ]
            self._orig_data = list(self._data)
            self._have_data = True
            self._fname = fname
            return True

    def saveData(self, fname):
        with open(fname, 'w') as fh:
            writer = csv.writer(fh)
            writer.writerow(self._header)
            writer.writerows(self._data)
            self._targetFile = fname

    def testHeader(self):
        return self._header == "x,y,sd".split(',') or \
            self._header == "Position X [nm],Position Y [nm],Precision [nm]".split(',') or \
            self._header == "X ,Y ,sd".split(',') or self._header == "X,Y,sd".split(',')

    def hasData(self):
        return self._have_data

    def x(self):
        if len(self._data) > 0:
            return zip(*self._data)[0]
        else:
            return []

    def y(self):
        if len(self._data) > 0:
            return zip(*self._data)[1]
        else:
            return []

    def sd(self):
        if len(self._data) > 0:
            return zip(*self._data)[2]
        else:
            return []

    def modifyData(self, operator, shape):
        self._data = operator(self._data, shape)
        return self

    def applyTo(self, operator, shape):
        return operator(self._data, shape)

class BlcComputer(QObject):
    ''' class that manages the data sets and the operations on them '''
    def __init__(self, args, parent=None):
        super(BlcComputer, self).__init__(parent)
        self._lineX = []
        self._lineY = []
        self._ops = BlcOps()
        self._dataSets = {}
        self._properties = {} # some options
        self._trees = {} # KdTrees by colour
        self._shapes = ShapeManger()

    def __getitem__(self, key):
        return self._dataSets[key]

    def shapes(self):
        return self._shapes

    def getTree(self, colour):
        ''' returns the k-d tree '''
        if not colour in self._trees:
            self._trees[colour] = MyKdTree(self._dataSets[colour])
        return self._trees[colour]

    def createDataSet(self, colour):
        self._dataSets[colour] = BlcData()
        if colour in self._trees:
            del self._trees[colour]

    def setProp(self, **kwargs):
        self._properties.update(kwargs)

    def prop(self, key):
        return self._properties[key]

    def addPoint(self, x, y):
        self._lineX.append(x)
        self._lineY.append(y)

    def remPoint(self):
        ''' removes the last point and returns the number of remaining points '''
        if len(self._lineX) == 0:
            return -1
        self._lineX.pop()
        self._lineY.pop()
        len(self._lineX)

    def findBoundingBoxGeom(self, geom):
        bbox = lambda x : Polygon([(x[0], x[1]), (x[2], x[1]), (x[2], x[3]), (x[0], x[3]) ])
        best_box = bbox(geom.bounds)
        best_angle = 0
        for angle in range(1, 180):
            poly=affinity.rotate(geom, angle, origin='centroid')
            bb = bbox(poly.bounds)
            if bb.area < best_box.area:
                best_box = bb
                best_angle = angle
        return affinity.rotate(best_box, -best_angle, origin=geom.centroid)

    def getObsCount(self, col, shape):
        ''' counts the observations of data set "col" in the shape given by "mode" '''
        count = self._ops['count']
        dataSet = self[col]
        return dataSet.applyTo(count, shape)

    def getCroppedDataFromShape(self, colour, shape):
        dataSet = self[colour]
        crop = self._ops['crop']
        header = dataSet.header()
        return BlcData.init(header, dataSet.applyTo(crop, shape) )

    def getParallelData(self, data, shape):
        angle, origin, _, _ = shape.angle_origin()
        data.modifyData(self._ops['shift'], origin)
        data.modifyData(self._ops['rotate'], -angle)
        return data

    def getGeneralRectangleAsPoly(self):
        ''' Returns the current shape as general rectangle '''
        return Shape(self.getGeneralRectangle().points(), Shape.RECTANGLE)

    def getGeneralRectangle(self):
        ''' calc a rectangle from the first 3 points '''
        num = self.numPoints()
        if num < 3:
            raise Exception('''"general rectangle" - Don't have enough points.''')
        return Rectangler(PointVector(self.lineX()[0], self.lineY()[0]),
                          PointVector(self.lineX()[1], self.lineY()[1]),
                          PointVector(self.lineX()[2], self.lineY()[2]))

    def getPolygon(self):
        ''' returns a polygon made from the points as vertices '''
        num = self.numPoints()
        if num < 3:
            raise Exception('''"polygon" - Don't have enough points.''')
        return Shape(zip(self.lineX(), self.lineY()), Shape.POLYGON)

    def clearPoints(self):
        self._lineX = []
        self._lineY = []

    def numPoints(self):
        return len(self.lineX())

    def lineX(self):
        return self._lineX

    def lineY(self):
        return self._lineY

class Modes(object):
    NOMODE, POLYDRAW = range(2)

MODEICON = ['applications-accessories', 'applications-development',
            'applications-office', 'applications-other',
            'applications-utilities']

class MessageWidget(QtGui.QWidget):
    ''' Basic Widget to display and save messages '''
    def __init__(self, parent=None):
        super(MessageWidget, self).__init__(parent)
        layout = QtGui.QVBoxLayout()
        self._textWidget = QtGui.QTextEdit()
        buttons = QtGui.QHBoxLayout()
        self._buttonSave = QtGui.QPushButton("Save")
        self._buttonSave.setIcon(QIcon(":/filesave.png"))
        self._buttonSave.setToolTip('Save messages to file.')
        self._buttonClear = QtGui.QPushButton("Clear")
        self._buttonClear.setIcon(QIcon(':/editclear.png'))
        self._buttonClear.setToolTip('Clear all messages.')
        buttons.addWidget(self._buttonSave)
        buttons.addWidget(self._buttonClear)
        buttons.addStretch()
        layout.addLayout(buttons)
        layout.addWidget(self._textWidget)
        self.setLayout(layout)
        self._connect()

    def _connect(self):
        self.connect(self._buttonClear, SIGNAL('clicked()'), self._textWidget.clear)
        self.connect(self._buttonSave, SIGNAL('clicked()'), self.onSave)

    def onSave(self):
        ''' save messages to file '''
        filename = QtGui.QFileDialog.getSaveFileName(self, "Save Messages")
        if filename:
            with open(filename, 'w') as fh:
                fh.write(self._textWidget.toPlainText())

    def append(self, *args):
        for item in args:
            self._textWidget.append(item)
            #self._textWidget.append('\n')

APPNAME = 'ROI Selection'

class AppForm(QtGui.QMainWindow):
    ''' user interface and main control '''
    MODES = ['Polygon', 'Rectangle']
    POLYMODE, GENRECTMODE = MODES
    BASEOUTPUTDIR = "/mnt/rclsfserv005/labinacell/bayes_cluster_analysis"
    def __init__(self, parent=None, args=None):
        QtGui.QMainWindow.__init__(self, parent)
        self.setWindowTitle(APPNAME)
        self._args = args
        self.artist = GmmPlotArtist(alpha=self._args.alpha)
        self._select_mode = Modes.NOMODE # this needs to be changed
        self._modeAction = None # Polygon, Rectangle etc.
        self.comp = BlcComputer(args)
        self.comp.setProp(squareLength=2000, dilate=DILATE, threshold=THRESHOLD)
        self.artist.set_computer(self.comp)
        self._logFile = None
        self._folderBase = None
        self._scatter = False # scatter plot OR imshow
        self._currentShape = None # what is currently shown
        self._shapes_saved = True
        self._path = os.getcwd()
        self._defIcon = QIcon(':/default.png')
        self._outputFolder = None

        self._densOpts = { 'colour' : COLOURS[0] }

        self.makeActions()
        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()
        self.createToolbar()
        self.createOptionsToolbar()
        self.create_dock1()
        #self.create_shape_dock()
        self.make_icons()

        self.connect(self.comp.shapes(), SIGNAL("storedShapesChanged"), self.onSavedState)
        #self.connect(self.shapeWidget, SIGNAL("currentRowChanged(int)"), self.selectShape)

        settings = QSettings("blc", "vis")
        self.restoreState(settings.value("MainWindow/State").toByteArray())

        if self._args.file and os.path.exists(self._args.file):
            self.loadDir(self._args.file)

    def shortCutDialog(self):
        dialog = QtGui.QDialog()
        dialog.setWindowTitle("Shortcuts")
        #dialog.setMinimumWidth(80)
        layout = QtGui.QVBoxLayout()
        textEdit = QtGui.QTextEdit()
        textEdit.setMinimumWidth(800)
        textEdit.setMinimumHeight(20)
        textEdit.setReadOnly(True)
        textEdit.setSizePolicy(QSizePolicy(QSizePolicy.MinimumExpanding,
                                           QSizePolicy.MinimumExpanding))
        buttons = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok)
        layout.addWidget(textEdit)
        for action in self.actions() + self._toolbar.actions():
            text = action.text().remove("&")
            shortcut = action.shortcut().toString()
            if len(text) and len(shortcut):
                string = QString("%1 (%2)").arg(text).arg(shortcut)
                tip = action.toolTip()
                if len(tip):
                    string.append(QString(": %1").arg(tip))
                textEdit.append(string)
        layout.addWidget(buttons)
        self.connect(buttons, SIGNAL("accepted()"), dialog, SLOT("accept()"))
        dialog.setLayout(layout)
        dialog.exec_()

    def make_icons(self):
        size = QSize(20*10, 20)
        self.pix_red = QPixmap(size)
        self.pix_red.fill(Qt.red)
        self.pix_grey = QPixmap(size)
        self.pix_grey.fill(Qt.gray)
        self.pix_green = QPixmap(size)
        self.pix_green.fill(Qt.green)

    def create_dock1(self):
        self.messageDock = QtGui.QDockWidget("Messages") # creates a dock widget with title 'Items'
        self.messageDock.setObjectName("MessagesDock") # this is required by QSettings to uniquely identify the dock widget
        self.messageWidg = MessageWidget() # class attribute to make referal easier
        self.messageDock.setWidget(self.messageWidg) # assign the list widget to the dock widget
        self.addDockWidget(Qt.LeftDockWidgetArea, self.messageDock) # add the dock widget to the main window (allows more options)

    def create_shape_dock(self):
        self.shapeDock =  QtGui.QDockWidget("Shapes")
        self.shapeDock.setObjectName("ShapeDock")
        self.shapeWidget = QtGui.QListWidget()
        self.shapeDock.setWidget(self.shapeWidget)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.shapeDock)

    def setModeAction(self, activate = True):
        self._modeAction = self.modeGroup.checkedAction().text()
        if activate:
            self.onSelection()

    def load_cmd(self):
        if self._folderBase is None:
            raise Exception("Why ...?")
        for colour in COLOURS:
            path = get_data_fname(self._folderBase, colour)
            self.comp.createDataSet(colour)
            if not path is None:
                if not self.comp[colour].loadData(path):
                    msg = 'File does not look like a "data.txt" file:\n{0}'.format(path)
                    QMessageBox.information(self, "Warning", msg)
        self.on_draw()
        self.setWindowTitle(APPNAME + " - " + os.path.basename(self._folderBase))
        self.status_text.setText(self._folderBase)
        self.message("Data set {0}".format(self._folderBase))

    def event(self, event):
        if event.type() == QEvent.KeyPress:
            if event.key() == Qt.Key_Control:
                self._holdControl = True
        if event.type() == QEvent.KeyRelease:
            if event.key() == Qt.Key_Control:
                self._holdControl = False
        return super(AppForm, self).event(event)

    def loadDir(self, path):
        ok = False
        path = path.rstrip('/')
        for colour in COLOURS:
            if path.endswith(colour):
                self._folderBase = path[0:-len(colour)]
                ok = True
                break
        if not ok:
            QMessageBox.warning(self, 'warning',
                '''folder doesn't seem to be a valid data folder "{0}"'''.format(path))
            return
        self.load_cmd()

    def load_dir_dlg(self):
        if not self._shapes_saved:
            response = QtGui.QMessageBox.question(self, "Data not saved",
                                "There are unsaved shapes. Save extracted data now?",
                                QMessageBox.Save, QMessageBox.Discard, QMessageBox.Cancel)
            if response == QMessageBox.Save:
                self.onSaveData()
                if not self._shapes_saved:
                    return
            elif response == QMessageBox.Discard:
                ''' nothing to do, just continue '''
            elif response == QMessageBox.Cancel:
                return
            else:
                raise Exception('Unknown response')
        path = unicode(QtGui.QFileDialog.getExistingDirectory(self, 'Load directory'))
        if path == "":
            return
        self.loadDir(path)

    def message(self, *args):
        self.messageWidg.append(*args)

    def onAnything(self, *args):
        self.message("Anything!!")
        for arg in args:
            self.message(arg)

    def on_about(self):
        msg = """
A GUI to select Region of interests.
Michael Hirsch. CLF, STFC. 2017
"""
        QMessageBox.about(self, "About", msg.strip())

    def highlightShape(self, x, y):
        #print("{0}, {1}".format(x,y))
        if not isinstance(x, (float, int)) or not isinstance(y, (float, int)):
            return
        point = Point(x, y)
        for idx, shape in enumerate(self.comp.shapes()):
            if shape.contains(point):
                self.comp.shapes().select(idx)
                self.canvas.draw()
                break

    def onSavedState(self):
        self._shapes_saved = len(self.comp.shapes()) == 0

    def onClick(self, event):
        ''' what happens if one clicks on the image '''
        if self._select_mode != Modes.POLYDRAW:
            ''' not in drawing mode '''
            self.highlightShape(event.xdata, event.ydata)
            return
        if self._modeAction == AppForm.POLYMODE:
            ''' any number of points are possible '''
            self.comp.addPoint(event.xdata, event.ydata)
            self.artist.drawPoint()
            self.canvas.draw()
            return
        if self._modeAction == AppForm.GENRECTMODE:
            ''' exactly three points shoulb be set '''
            self.comp.addPoint(event.xdata, event.ydata)
            if self.comp.numPoints() < 3:
                self.artist.drawPoint()
                self.canvas.draw()
                return
            else:
                self.onRectangle2()

    def onCount(self):
        ''' shows a message with number of observation per channel
        for the whole data set.'''
        stringList = ["observations:"]
        for col in COLOURS:
            size = len(self.comp[col])
            stringList.append("{0} -> {1}".format(col, size))
        self.message("Observation count", "\n".join(stringList))

    def onPolygonProper(self):
        self.artist.clearLines()
        if self.comp.numPoints() < 3:
            return
        self.addShape(self.comp.getPolygon())
        self._currentShape = AppForm.POLYMODE
        self.countObs()
        self.canvas.draw()

    def _density2multipoly(self, points, densities, threshold, buffer_size):
        ''' creates a MultiPolygon from a list of points given a threshold for a  density '''
        segments = []

        for pt, dens in zip(points, densities):
            if dens <= threshold:
                continue
            segments.append(Point(*pt).buffer(buffer_size))
        return unary_union(segments)

    def _density2convexhull(self, points, densities, threshold, buffer_size):
        segments = MultiPoint([ pt for pt, dens in zip(points, densities) if dens > threshold ])
        return segments.convex_hull.buffer(buffer_size)

    def onGetKdtree(self):
        for colour in COLOURS:
            tree = self.comp.getTree(colour)
            self.message("Kdtree for {0}, size = {1}".format(colour, len(tree)))

    def onDensitySegmentation(self):
        shapes = self.comp.shapes().highlighted()
        if not shapes:
            QtGui.QMessageBox.warning(self, "Warning", "No shape!")
            return
        shape = shapes[0]
        MAXSIZE = 100000
        sigma = (float(DILATE)/2.0)/2.0
        kernel = GaussianTypeKernel(sigma)
        all_pts = []
        all_dens = []
        for colour in COLOURS:
            data = self.comp.getCroppedDataFromShape(colour, shape)
            size = len(data)
            if size > MAXSIZE:
                QtGui.QMessageBox.warning(self, "Warning", "ROI to large: {0}. Maximum = {1}".format(
                    size, MAXSIZE ))
                return
            if size == 0:
                continue
            points = np.array([ [x, y] for x, y, _ in data ])
            kdtree = self.comp.getTree(colour)
            #kdtree = MyKdTree(data) # dodgy, goes wrong if order of points is not preserved
            densities = []

            for pnt in points:
                idx = kdtree.ballPoint(pnt, r=DILATE)
                #densities.append( sum([ kernel(pnt, x) for x in points[idx] ]) )
                densities.append( sum([ kernel(pnt, x) for x in self.comp[colour].points(idx) ]) )
            all_pts.extend(points)
            all_dens.extend(densities)
        if len(all_pts) == 0:
            return

        ''' turn points into shapes '''
        buffer_size = sigma/5.0

        #combined = self._density2multipoly(all_pts, all_dens, THRESHOLD, buffer_size)
        threshold=self.comp.prop('threshold')
        self.message("Threshold={0}".format(threshold))
        combined = self._density2convexhull(all_pts, all_dens, threshold, buffer_size)


        new_shape = Shape(self.comp.findBoundingBoxGeom(combined).exterior, Shape.RECTANGLE)
        ''' adopt changes '''
        self.comp.shapes().remove(shape)
        self.addShape(new_shape)
        shape.remove()
        self.canvas.draw()

    def countObs(self):
        ''' Counts the observations in the current shape '''
        if len(self.comp.shapes().highlighted()) == 0:
            QMessageBox.warning(self, 'Note:', 'No highlighted shapes.')
            return
        stringList = ["observations in highlighted:"]
        for shape in self.comp.shapes().highlighted():
            if shape.shape_type() == Shape.RECTANGLE:
                angle, origin, width, height = shape.angle_origin(unit='deg')
                stringList.append(u"origin=({0:.1f}, {1:.1f}), angle={2:.1f}\u00b0".format(
                    origin[0], origin[1], angle))
                stringList.append("{0:.1f}x{1:.1f}={2:.1f}".format(width, height, shape.area))
            else:
                stringList.append("bounds = {0:.1f}, {1:.1f}, {2:.1f}, {3:.1f}".format(
                    *shape.bounds))
                stringList.append("area = {0:.1f}".format(shape.area))
            for col in COLOURS:
                size = self.comp.getObsCount(col, shape)
                stringList.append(u"{0} = {1} obs; density = {2:.2f}/\u03bcm\u00b2".format(
                    col, size, NANOMICRON*size/shape.area))
        self.message(*stringList)

    def selectShape(self, idx):
        self.comp.shapes().select(idx)
        self.canvas.draw()

    def addShape(self, shape):
        self.comp.shapes().push(shape)
        shape.draw(self.artist)
        self.comp.shapes().select(-1)

    def onRemoveShape(self):
        self.comp.shapes().popHighlighted()
        self.canvas.draw()

    def onRectangle2(self):
        self.setMode(Modes.NOMODE)
        self.artist.clearLines()
        if self.comp.numPoints() < 3:
            return
        self.addShape(self.comp.getGeneralRectangleAsPoly())
        self.canvas.draw()
        self._currentShape = AppForm.GENRECTMODE
        self.countObs()

    def onGrid(self):
        self.axes.grid(self.grid_cb.isChecked())
        self.canvas.draw()

    def onAlpha(self):
        self.artist.set_alpha(self.alphaSpinBox.value())
        self.on_draw()

    def on_draw(self, **kwargs):
        """ Redraws the figure
        """
        if not self._shapes_saved:
            response = QMessageBox.question(self, "Warning",
                        "Data not saved. Continue?", QMessageBox.Yes, QMessageBox.Abort)
            if response == QMessageBox.Yes:
                ''' nothing to do; continue '''
            elif response == QMessageBox.Abort:
                return
            else:
                raise Exception("Unknown response")
        self.comp.clearPoints()
        self.artist.clearLines()
        self.comp.shapes().clear()
        self.axes.clear()
        self.axes.grid(self.grid_cb.isChecked())

        if self._scatter:
            self.artist.plot_sampleplot()
        else:
            self.artist.plot_scatter_as_image()
        self.artist.set_title(self._folderBase)
        self.mpl_toolbar.push_current()
        self.mpl_toolbar.draw()

    def setMode(self, mode):
        ''' sets and unsets the drawing mode '''
        self._select_mode = mode
        if mode == Modes.NOMODE:
            self.selectAction.setIcon(QIcon(':/start.png'))
            self.selectAction.setText("Start")
            self._indicator.setPixmap(QIcon(':/face-smile.png').pixmap(32, 32))
            self._indicator.setToolTip("not in drawing mode")
            self._indicator.setStatusTip("not in drawing mode")
        if mode == Modes.POLYDRAW:
            self.selectAction.setIcon(QIcon(':/stop.png'))
            self.selectAction.setText("Finish drawing")
            self._indicator.setPixmap(QIcon(':/face-cool.png').pixmap(32, 32))
            self._indicator.setToolTip("Drawing mode")
            self._indicator.setStatusTip("Drawing mode")

    def _setOutputFolder(self):
        username = os.environ['USER']
        base = AppForm.BASEOUTPUTDIR
        folder = os.path.join(base, username)
        if not os.path.isdir(folder):
            os.mkdir(folder)
        self._outputFolder = unicode(QtGui.QFileDialog.getExistingDirectory(
            caption="Select base directory for output data", directory=folder))
        if not os.path.isdir(self._outputFolder):
            self._outputFolder = None
        else:
            self.message('Output folder = {0}'.format(self._outputFolder))
        return

    def onSaveData(self):
        ''' the general "Save" proxy. Calls specific save functions '''
        if self._outputFolder is None:
            self._setOutputFolder()
        if self._outputFolder is None:
            QMessageBox.information(self, "Warning", "No valid output directory selected.")
            return
        for shape in self.comp.shapes():
            self.onShapeSave(shape)

    def onSelection(self):
        if self._select_mode == Modes.POLYDRAW:
            self.setMode(Modes.NOMODE)
            if self._modeAction == AppForm.POLYMODE:
                self.onPolygonProper()
                return
        else:
            self.onClear() # This needs to come before the mode setting
            self.setMode(Modes.POLYDRAW)

    def onShapeSave(self, shape):
        if shape.shape_type() != Shape.RECTANGLE:
            return
        if self._outputFolder is None:
            QMessageBox.information(self, "warning", "No output directory")
            return
        for colour in COLOURS:
            data = self.comp.getCroppedDataFromShape(colour, shape)
            if len(data) == 0:
                self.message('Warning', 'No data to save for {0}.'.format(colour))
                continue
            data = self.comp.getParallelData(data, shape)
            angle, origin, width, height = shape.angle_origin()
            sign = "r{0}-{1}_{2}x{3}_{4:.0f}".format(*map(
                int, [origin[0], origin[1], width, height, angle*180/np.pi]))
            #sign = "x{0}y{1}d{2}a{3}".format(*map(
            #    int, [origin[0], origin[1], angle*180/np.pi, shape.area]))
            basefolder = os.path.split(self._folderBase)[1] + sign + colour
            basefolder = os.path.join(self._outputFolder, basefolder)
            folder = os.path.join(basefolder, str(1))
            if not os.path.isdir(folder):
                os.makedirs(folder)
            assert(os.path.isdir(folder))
            #
            shape_log = "shape_log.json"
            log_data = {
                "type" : shape.shape_type(),
                "exterior" : list(shape.exterior.coords)[0: -1]
            }
            with open(os.path.join(basefolder, shape_log), 'w') as fh:
                json.dump(log_data, fh)
            #
            fname = os.path.join(folder, "data.txt")
            data.saveData(fname)
            writeShapeConfig(basefolder, width, height)
            self.message("{0} points saved to {1}".format(len(data), fname))
        self._shapes_saved = True

    def _remPoint(self):
        remaining = self.comp.remPoint()
        if remaining == -1:
            self.onSelection()
            return
        self.artist.drawPoint()
        self.canvas.draw()

    def onClear(self):
        ''' removes the polygon without cutting '''
        if self._select_mode == Modes.POLYDRAW:
            self._remPoint()
            return
        if self._select_mode == Modes.NOMODE:
            self.comp.clearPoints()
            self.artist.clearLines()
            self.canvas.draw()
            self._currentShape = None

    def onSaveLog(self, state):
        if not state == Qt.Checked:
            return
        path = QtGui.QFileDialog.getSaveFileName(self, 'Log file name',
                                           options = QtGui.QFileDialog.DontConfirmOverwrite)
        if len(path) > 0:
            self._logFile = unicode(path)

    def onScatter(self, state):
        if state == Qt.Checked:
            self._scatter = True
        elif state == Qt.Unchecked:
            self._scatter = False
        else:
            return
        self.on_draw()

    def makeActions(self):
        self.modeGroup = QtGui.QActionGroup(self)
        for name, icon in zip(AppForm.MODES, MODEICON):
            action = self.create_action(name, self.setModeAction, checkable=True, icon=icon)
            self.modeGroup.addAction(action)
        defaultIndex = AppForm.MODES.index(AppForm.GENRECTMODE)
        self.modeGroup.actions()[defaultIndex].setChecked(True)
        self.setModeAction(False)

        self.obsCountAction = self.create_action("Observation count", slot=self.onCount,
                            icon='info', tip=" ".join(self.onCount.__doc__.split()))

        self.densitySegmentationAction = self.create_action("Improve ROI",
                            slot=self.onDensitySegmentation, shortcut = 's',
                            tip='Find the optimal rectangle for the highlighted shape',
                            icon='preferences-desktop-theme')

        self.countObsAction = self.create_action("Shape counts",
                            slot=self.countObs,
                            tip='Message size and observation count of highlighted shapes',
                            icon='face-glasses')

        self.saveRectAction = self.create_action("Save data", self.onSaveData,
            shortcut="Ctrl+S",
            icon='filesaveas', tip='Save the data in the rectangular ROIs')

        self.reloadAction = self.create_action("Redraw", self.on_draw, icon='reload',
                    tip='Clears and redraws the plot.')

        self.load_file_action = self.create_action("&Load data file",
            shortcut="Ctrl+L", slot=self.load_dir_dlg,
            tip="Load a file that contains super-resolution data",
            icon = 'fileopen')

        self.quit_action = self.create_action("&Quit", slot=self.close,
            shortcut="Ctrl+Q", tip="Quit the application", icon='exit')

        self.selectAction = self.create_action("Start", slot = self.onSelection,
                            tip = "Start/Stop drawing ROIs", shortcut = "d")
        self.selectAction.setIcon(QIcon(':/start.png'))
        self.remShapeAction = self.create_action("Delete Shape", slot=self.onRemoveShape,
                            tip='Delete highlighted shape(s)', shortcut='Del',
                            icon='edit-delete')

        self.clearAction = self.create_action("Clear", slot = self.onClear,
                            tip = 'Clear the last point.',
                            shortcut = "Esc", icon='edit-clear')
        self.selectOutputFolderAction = self.create_action('Set base output directory',
                                                    self._setOutputFolder,
                tip="Data will be stored in sub-directories of this directory.")

    def createToolbar(self):
        self._toolbar = QtGui.QToolBar("Main Toolbar")
        self._toolbar.setObjectName("mainToolBar")

        self._toolbar.addAction(self.load_file_action)
        self._toolbar.addAction(self.saveRectAction)
        self._toolbar.addSeparator()
        self._toolbar.addAction(self.selectAction)
        self._toolbar.addSeparator()
        self._toolbar.addActions(self.modeGroup.actions())
        self._toolbar.addSeparator()
        self._toolbar.addAction(self.clearAction)
        self._toolbar.addAction(self.reloadAction)
        self._toolbar.addSeparator()
        self._toolbar.addAction(self.densitySegmentationAction)
        self._toolbar.addAction(self.countObsAction)
        self._toolbar.addAction(self.remShapeAction)
        self._toolbar.addSeparator()
        self._toolbar.addAction(self.quit_action)
        self.addToolBar(self._toolbar)

    def createOptionsToolbar(self):
        self._optionToolbar = QtGui.QToolBar("Options Toolbar")
        self._optionToolbar.setObjectName("optionsToolBar")
        self._optionToolbar.setAllowedAreas(Qt.TopToolBarArea|Qt.BottomToolBarArea)

        self.grid_cb = QtGui.QCheckBox("Show &Grid")
        self.grid_cb.setChecked(False)
        self.grid_cb.setToolTip("De-/Activate the plot grid")
        self.grid_cb.setStatusTip("De-/Activate the plot grid")
        self.connect(self.grid_cb, SIGNAL('stateChanged(int)'), self.onGrid)

        self.alphaSpinBox = QtGui.QDoubleSpinBox()
        self.alphaLabel = QtGui.QLabel("Opacity: ")
        self.alphaSpinBox.setRange(0.0, 1.0)
        self.alphaSpinBox.setDecimals(2)
        self.alphaSpinBox.setSingleStep(0.01)
        self.alphaSpinBox.setValue(self.artist.alpha())
        self.alphaSpinBox.setToolTip("Opacity of the plot")
        self.alphaSpinBox.setStatusTip("Opacity of the plot")
        self.connect(self.alphaSpinBox, SIGNAL("editingFinished()"), self.onAlpha)

        self.threshSpinBox = QtGui.QDoubleSpinBox()
        self.threshLabel = QtGui.QLabel("Density threshold: ")
        self.threshLabel.setBuddy(self.threshSpinBox)
        self.threshSpinBox.setRange(0.0, 1.0)
        self.threshSpinBox.setDecimals(3)
        self.threshSpinBox.setSingleStep(0.001)
        self.threshSpinBox.setValue(THRESHOLD)
        self.threshSpinBox.setToolTip('Treshold for ROI optimisation.\n'
            'Larger values result in smaller ROIs')
        self.threshSpinBox.setStatusTip('Treshold for ROI optimisation. '
            'Larger values result in smaller ROIs')
        self.connect(self.threshSpinBox, SIGNAL("valueChanged(double)"),
                         lambda x : self.comp.setProp(threshold=x))
        self.threshButton = QtGui.QPushButton("Default")
        self.threshButton.setToolTip("Set the density treshold to the default value")
        self.threshButton.setStatusTip("Set the density treshold to the default value")

        self.connect(self.threshButton, SIGNAL("clicked()"),
                         partial(self.threshSpinBox.setValue, THRESHOLD))



        self.scatterCheckBox = QtGui.QCheckBox('Scatter')
        self.scatterCheckBox.setToolTip('Use scatter plot (slower but higher resolution)')
        self.scatterCheckBox.setStatusTip('Use the slower scatter plot')
        self.connect(self.scatterCheckBox, SIGNAL("stateChanged(int)"), self.onScatter)

        self._indicator = QtGui.QLabel()
        self._indicator.setPixmap(QIcon(':/face-smile.png').pixmap(32, 32))
        self._indicator.setToolTip("not in drawing mode")
        self._indicator.setStatusTip("not in drawing mode")

        for w in [ self.grid_cb, self.alphaLabel, self.alphaSpinBox, self.scatterCheckBox,
                  self.threshLabel, self.threshSpinBox, self.threshButton]:
            self._optionToolbar.addWidget(w)
            #hbox.setAlignment(w, Qt.AlignLeft)
        self._optionToolbar.addSeparator()
        self._optionToolbar.addWidget(self._indicator)
        #hbox.addStretch()

        self.addToolBar(Qt.BottomToolBarArea, self._optionToolbar)


    def create_main_frame(self):
        self.main_frame = QtGui.QWidget()
        self.cb = None

        # Create the mpl Figure and FigCanvas objects.
        # 5x4 inches, 100 dots-per-inch
        #
        self.dpi = 100
        self.fig = Figure((8.0, 6.4), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)

        self.axes = self.fig.add_subplot(111)
        self.artist.set_axes(self.axes)

        self.canvas.mpl_connect('button_release_event', self.onClick)
        # Create the navigation toolbar, tied to the canvas
        #
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)

        action = QAction("Canvas Pan", self)
        action.setShortcut(Qt.Key_F2)
        self.connect(action, SIGNAL("triggered()"), self.mpl_toolbar.pan)
        self.addAction(action)

        action = QAction("Canvas Zoom", self)
        action.setShortcut(Qt.Key_F3)
        self.connect(action, SIGNAL("triggered()"), self.mpl_toolbar.zoom)
        self.addAction(action)

        action = QAction("Canvas Home", self)
        action.setShortcut(Qt.Key_F1)
        self.connect(action, SIGNAL("triggered()"), self.mpl_toolbar.home)
        self.addAction(action)

        action = QAction("Canvas Back", self)
        action.setShortcut(Qt.Key_Back)
        self.connect(action, SIGNAL("triggered()"), self.mpl_toolbar.back)
        self.addAction(action)

        action = QAction("Canvas Forward", self)
        action.setShortcut(Qt.Key_Forward)
        self.connect(action, SIGNAL("triggered()"), self.mpl_toolbar.forward)
        self.addAction(action)

        self._navToolBar = QtGui.QToolBar("Naviagtion Toolbar")
        self._navToolBar.setObjectName("navToolBar")
        self._navToolBar.addWidget(self.mpl_toolbar)
        self._navToolBar.setAllowedAreas(Qt.TopToolBarArea|Qt.BottomToolBarArea)
        self.addToolBar(self._navToolBar)

        # Other GUI controls
        #


        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.canvas, stretch = 1)
        #vbox.addLayout(hbox)

        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)

    def create_paradock_widget(self):
        self.paradock = QtGui.QDockWidget("Files")
        self.paralist = QtGui.QListWidget()
        self.paradock.setWidget( self.paralist )
        self.addDockWidget( Qt.RightDockWidgetArea, self.paradock )
        self.connect(self.paralist, SIGNAL('itemClicked (QListWidgetItem *)'), self.on_filesel)


    def create_status_bar(self):
        self.status_text = QtGui.QLabel("Baysian localisation clustering tool.")
        self.statusBar().addWidget(self.status_text, 1)

    def create_menu(self):
        # file menu
        self.file_menu = self.menuBar().addMenu("&File")


        self.add_actions(self.file_menu,
            (self.load_file_action, self.saveRectAction, self.selectOutputFolderAction,
             None, self.quit_action))

        # selection menu
        self.selection_menu = self.menuBar().addMenu("&ROI")
        self.add_actions(self.selection_menu,
                         (self.selectAction, self.clearAction, None,
                          self.densitySegmentationAction,
                          self.obsCountAction))

        # help menu
        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About",
            shortcut=None, slot=self.on_about,
            tip='About the program')

        action = QAction("Keyboard short cuts", self)
        action.setStatusTip('Shows a list of the keybord short-cuts')
        self.connect(action, SIGNAL('triggered()'), self.shortCutDialog)
        self.add_actions(self.help_menu, (action, about_action,))

    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_action(self, text, slot=None, shortcut=None,
                        icon=None, tip=None, checkable=False,
                        signal="triggered()"):
        action = QAction(text, self)
        if icon is not None:
            if isinstance(icon, QIcon):
                action.setIcon(icon)
            else:
                action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action

    def closeEvent(self, event):
        ''' the overridden close event '''
        settings = QSettings("blc", "vis") # open the config file (~/.config/qt_demo/demo.conf). The concrete path depends on the system and options given
        settings.setValue("MainWindow/State", QVariant(self.saveState()))  # store the current states of windows

class GmmPlotArtist(object):
    ''' class that creates the plots (but doesn't render it) '''
    def __init__(self, **kwargs):
        self.axes  = kwargs.get("axes",  None)
        self.comp  = kwargs.get("comp",  None)
        self.set_alpha(kwargs.get("alpha", None))
        self._lines = []
        #self._mainPlots = []
        self._dpi = 80
        self._width = 8*2
        self._height = 6*2
        self._storedPatches = []

    def set_axes(self, axes):
        self.axes = axes
    def set_computer(self, comp):
        self.comp = comp

    def alpha(self):
        return self._alpha
    def set_alpha(self, alpha):
        if alpha is None or 1.0 < alpha or alpha <= 0.0:
            self._alpha = 1.0
        else:
            self._alpha = alpha

    def clearLines(self):
        while len(self._lines):
            self._lines.pop().remove()

    def drawPoint(self):
        self.clearLines()
        self._lines = self.axes.plot(self.comp.lineX(), self.comp.lineY(), marker = 'x', color='k')

    def add_patch(self, patch):
        self.axes.add_patch(patch)

    def plot_scatter_as_image(self):
        ''' Displays the scatter plot as image
        '''
        plt.clf()
        for colour in COLOURS:
            data = self.comp[colour]
            if not data.hasData():
                continue
            plt.scatter( data.x(), data.y(), s=1, alpha=self.alpha(), color=colour)
        plt.axis('off')
        plt.subplots_adjust(0, 0, 1, 1, 0, 0)
        plt.gcf().set_figwidth(self._width)
        plt.gcf().set_figheight(self._height)
        stream = BytesIO()
        xmin, xmax = plt.xlim()
        ymin, ymax = plt.ylim()
        plt.savefig(stream, dpi=self._dpi, format='raw')
        pilImage = Image.frombytes('RGBA',size=(self._width * self._dpi, self._height * self._dpi),
                                    data = stream.getvalue())
        self.axes.imshow(image.pil_to_array(pilImage), extent = (xmin, xmax, ymin, ymax))
        self.axes.set_xlabel("x (nm)")
        self.axes.set_ylabel("y (nm)")

    def set_title(self, title):
        self.axes.set_title(title)

    def plot_sampleplot(self):
        ''' Creates a scatter plot
        '''
        for colour in COLOURS:
            data = self.comp[colour]
            if not data.hasData():
                continue
            self.axes.scatter( data.x(), data.y(), s=2, alpha=self.alpha(), color=colour)
        self.axes.set_xlabel("x (nm)")
        self.axes.set_ylabel("y (nm)")


def arg_parser():
    parser = argparse.ArgumentParser(
        description=globals()['__doc__'],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter    )
    parser.add_argument("file", nargs='?', help='data folder (e.g. "s20a03red")')
    parser.add_argument("--alpha", "-a", type=float, default=0.3,
        help = 'Transparency factor that is used in some plots. Range (0, 1]')
    return parser.parse_args()

def main():
    args = arg_parser()
    app = QtGui.QApplication(sys.argv)
    app.setWindowIcon(QIcon(":/applications-accessories.png"))
    form = AppForm(args=args)
    form.show()
    app.exec_()

if __name__ == "__main__":
    main()
