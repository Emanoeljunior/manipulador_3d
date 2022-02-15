from PyQt5 import QtWidgets, QtCore, QtGui
import sys
from math import pi
import simpy.rt
import model
import control
import numpy as np
import math
import pyqtgraph.opengl as gl
from skspatial.objects import Plane
from skspatial.objects import Point

import ui_main
import pyqtgraph
from threading import Thread

class ExampleApp(QtWidgets.QMainWindow, ui_main.Ui_MainWindow):
    signal = QtCore.pyqtSignal()
    def __init__(self, parent=None):
        pyqtgraph.setConfigOption('background', 'w') #before loading widget
        pyqtgraph.setConfigOption('antialias', True)
        super(ExampleApp, self).__init__(parent)
        self.setupUi(self)

        self.mainTh = Thread(target=self.simTh)

        self.grPlot.setBackgroundColor('w')

        xygrid = gl.GLGridItem(color=(200, 200, 200, 200))
        xygrid.scale(20, 20, 20)
        xygrid.setSpacing(0.01, 0.01, 0.01)
        self.grPlot.addItem(xygrid)

        # controller object
        self.ctrl = control.control()

        # model object
        iniPos = [0, 45, -45]
        self.mdl = model.model([iniPos[0] * pi / 180, iniPos[1] * pi / 180, iniPos[2] * pi / 180])
        self.gripRef = None
        self.jointRef = None

        self.th1Slider.valueChanged.connect(self.slidersChange)
        self.th1Slider.setValue(iniPos[0])
        self.th2Slider.valueChanged.connect(self.slidersChange)
        self.th2Slider.setValue(iniPos[1])
        self.th3Slider.valueChanged.connect(self.slidersChange)
        self.th3Slider.setValue(iniPos[2])

        self.setRefBtn.clicked.connect(self.setRef)
        self.dynBtn.clicked.connect(self.dynClick)
        self.ctrlBtn.clicked.connect(self.ctrlClick)

        self.signal.connect(self.updateInterface)

        self.mainTh.start()


    def simulation(self, env):
        tm = 0
        while True:
            # verify if dynamics is enable
            if self.dynBtn.isChecked():
                if self.ctrlBtn.isChecked():
                    tau = self.ctrl.control(self.mdl, self.jointRef)
                    
                else:
                    tau = np.array([[0, 0, 0]]).transpose()

                thdd = self.mdl.dynamics(tau)
               
                thdd = tau
            else:
                self.mdl.thd = np.array([[0, 0, 0]]).transpose()
                thdd = np.array([[0, 0, 0]]).transpose()
            
            print("thdd")
            print(thdd)
            self.mdl.integrate(thdd)

            if tm > 0.1:
                self.signal.emit()
                tm = 0

            tm = tm + self.mdl.dt

            yield env.timeout(self.mdl.dt)

    def simTh(self):
        # simulation environment
        self.env = simpy.rt.RealtimeEnvironment(factor=1, strict=0)
        # simulation process definition
        proc = self.env.process(self.simulation(self.env))
        # simulation start
        self.env.run()

    def draw(self, pointList, ref=None):
        armPoints = []
        projPoints = []

        for idx in range(0, len(pointList)-1):
            p0 = np.array([pointList[idx]])
            p1 = np.array([pointList[idx+1]])
            res = self.intersect(p0, p1)
            if res is not None:
                if (p0[0][2] < p1[0][2]):
                    armPoints.append([np.concatenate((p0, res), axis=0), -1])
                    armPoints.append([np.concatenate((res, p1), axis=0), 2])
                    projPoints.append(self.z_proj(res.squeeze(), p1.squeeze()))
                else:
                    armPoints.append([np.concatenate((p0, res), axis=0), 2])
                    armPoints.append([np.concatenate((res, p1), axis=0), -1])
                    projPoints.append(self.z_proj(p0.squeeze(), res.squeeze()))
            else:
                if (p0[0][2] < 0):
                    armPoints.append([np.concatenate((p0, p1), axis=0), -1])
                else:
                    armPoints.append([np.concatenate((p0, p1), axis=0), 2])
                    projPoints.append(self.z_proj(p0.squeeze(), p1.squeeze()))

        while(len(self.grPlot.items)) > 1:
            self.grPlot.removeItem(self.grPlot.items[-1])

        for idx, elem in enumerate(armPoints, start=1):
            item = gl.GLLinePlotItem(pos=elem[0], width=10, color=QtGui.QColor(0,50,idx*int(255/(len(armPoints))),255), antialias=True)
            item.setDepthValue(elem[1])
            self.grPlot.addItem(item)

        for elem in projPoints:
            item = gl.GLLinePlotItem(pos=elem, width=5, color=QtGui.QColor(150,150,150,100), antialias=True)
            item.setDepthValue(-1)
            self.grPlot.addItem(item)

        if ref is not None:
            md = gl.MeshData.sphere(rows=10, cols=20, radius=0.05)

            item = gl.GLMeshItem(meshdata=md, smooth=True, color=QtGui.QColor(200,0,0,255))
            item.translate(*ref)
            item.setDepthValue(2)
            self.grPlot.addItem(item)

    def z_proj(self, p0, p1):

        plane = Plane(point=[0, 0, 0], normal=[0, 0, 1])

        p1_p = plane.project_point(Point(p0)).to_array()
        p2_p = plane.project_point(Point(p1)).to_array()

        return np.concatenate(([p1_p], [p2_p]), axis=0)

    def intersect(self, p0, p1, epsilon=1e-6):
        """
        p0, p1: Define the line.
        p_co, p_no: define the plane:
            p_co Is a point on the plane (plane coordinate).
            p_no Is a normal vector defining the plane direction;
                 (does not need to be normalized).

        Return a Vector or None (when the intersection can't be found).
        """
        p_co = np.array([[10, 10, 0]])
        p_no = np.array([[0, 0, 1]])

        u = np.subtract(p1, p0)
        dot = np.dot(np.squeeze(p_no), np.squeeze(u))

        if abs(dot) > epsilon:
            # The factor of the point between p0 -> p1 (0 - 1)
            # if 'fac' is between (0 - 1) the point intersects with the segment.
            # Otherwise:
            #  < 0.0: behind p0.
            #  > 1.0: infront of p1.
            w = np.subtract(p0, p_co)
            fac = -np.dot(np.squeeze(p_no), np.squeeze(w)) / dot
            p2 = np.add(p0, fac * u)

            # verifica se o p2 pertence ao elo
            AB = math.dist(p0.squeeze(), p1.squeeze())
            AC = math.dist(p0.squeeze(), p2.squeeze())
            CB = math.dist(p2.squeeze(), p1.squeeze())

            if abs(AB - (AC+CB)) < epsilon:
                return p2

        return None

    def updateInterface(self):

        points = np.vstack(self.mdl.kinematics()).transpose()

        self.draw(points, self.gripRef)

        self.th1Label.setText("{:.2f}".format(self.mdl.th[0][0]))
        self.th2Label.setText("{:.2f}".format(self.mdl.th[1][0]))
        self.th3Label.setText("{:.2f}".format(self.mdl.th[2][0]))

        self.xLabel.setText("{:.2f}".format(points[-1][0]))
        self.yLabel.setText("{:.2f}".format(points[-1][1]))
        self.zLabel.setText("{:.2f}".format(points[-1][2]))

    def slidersChange(self):
        self.mdl.th[0][0] = self.th1Slider.value() * pi / 180
        self.th1SliderLabel.setText(str(self.th1Slider.value()) + 'º')
        self.mdl.th[1][0] = self.th2Slider.value() * pi / 180
        self.th2SliderLabel.setText(str(self.th2Slider.value()) + 'º')
        self.mdl.th[2][0] = self.th3Slider.value() * pi / 180
        self.th3SliderLabel.setText(str(self.th3Slider.value()) + 'º')

    def setRef(self):
        def is_number(string):
            try:
                float(string)
                return True
            except ValueError:
                return False

        print("set ref")
        if is_number(self.xRefEdit.text()) and is_number(self.yRefEdit.text()) and is_number(self.zRefEdit.text()):
            self.xRef.setText(self.xRefEdit.text())
            self.yRef.setText(self.yRefEdit.text())
            self.zRef.setText(self.zRefEdit.text())
            xRef = float(self.xRef.text())
            yRef = float(self.yRef.text())
            zRef = float(self.zRef.text())
            self.gripRef = np.array([[xRef, yRef, zRef]]).transpose()
            kin = np.vstack(self.mdl.kinematics()).transpose()
            self.jointRef = self.mdl.invKinematics(self.gripRef)
            print("join")
            print(self.gripRef)
            print("join")
            self.xRefEdit.clear()
            self.yRefEdit.clear()
            self.zRefEdit.clear()



    def dynClick(self):
        state = self.dynBtn.isChecked()

        self.th1Slider.setEnabled(not state)
        self.th2Slider.setEnabled(not state)
        self.th3Slider.setEnabled(not state)
        self.ctrlBtn.setEnabled(state)

        if not state:
            self.xRefEdit.setEnabled(False)
            self.yRefEdit.setEnabled(False)
            self.zRefEdit.setEnabled(False)
            self.setRefBtn.setEnabled(False)
            self.setRefBtn.setChecked(False)
            self.xRef.setText('-')
            self.yRef.setText('-')
            self.zRef.setText('-')
            self.xRefEdit.clear()
            self.yRefEdit.clear()
            self.zRefEdit.clear()
            self.gripRef = None

    def ctrlClick(self):
        state = self.ctrlBtn.isChecked()

        self.xRefEdit.setEnabled(state)
        self.yRefEdit.setEnabled(state)
        self.zRefEdit.setEnabled(state)
        self.setRefBtn.setEnabled(state)

        if not state:
            self.xRef.setText('-')
            self.yRef.setText('-')
            self.zRef.setText('-')
            self.xRefEdit.clear()
            self.yRefEdit.clear()
            self.zRefEdit.clear()
            self.gripRef = None
        else:
            ref = np.vstack(self.mdl.kinematics()).transpose()[-1]
            self.gripRef = np.array([ref]).transpose()
            self.jointRef = self.mdl.th
            self.xRef.setText("{:.2f}".format(ref[0]))
            self.yRef.setText("{:.2f}".format(ref[1]))
            self.zRef.setText("{:.2f}".format(ref[2]))


if __name__=="__main__":
    app = QtWidgets.QApplication(sys.argv)
    form = ExampleApp()
    form.show()
    form.update() #start with something
    app.exec_()
    print("DONE")
