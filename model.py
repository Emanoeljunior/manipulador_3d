from math import sin, cos, pi
from re import A
import numpy as np
import math


class model(object):
    def __init__(self, th=[0, 0, 0], dt=0.002):

        self.th = np.array([th]).transpose()
        self.thd = np.array([0, 0, 0]).transpose()

        self.l1 = 2
        self.l2 = 2
        self.l3 = 1

        self.m1 = 1
        self.m2 = 1
        self.m3 = 1

        self.g = 9.8

        self.dt = dt

    def integrate(self, thdd=np.array([[0, 0, 0]]).transpose()):

        # thdd=np.array([[0, 0, 0]]).transpose()
        self.thd = np.add(self.thd, thdd*self.dt)
        self.th = np.add(self.th, self.thd*self.dt)

    def kinematics(self, th=None):

        if th is None:
            th = self.th

        x0 = 0
        y0 = 0
        z0 = 0

        x1 = 0
        y1 = 0
        z1 = self.l1

        x2 = self.l2 * (cos(th[0][0]) * cos(th[1][0]))
        y2 = self.l2 * (sin(th[0][0]) * cos(th[1][0]))
        z2 = z1 + self.l2 * sin(th[1][0])

        x3 = x2 + self.l3 * cos(th[0][0]) * cos(th[1][0] + th[2][0])
        y3 = y2 + self.l3 * sin(th[0][0]) * cos(th[1][0] + th[2][0])
        z3 = z2 + self.l3 * sin(th[1][0] + th[2][0])

        return [[x0, x1, x2, x3], [y0, y1, y2, y3], [z0, z1, z2, z3]]

    def invKinematics(self, p = np.array([[0], [0], [0]])):
        # implementar a cinemática
        xi = np.array([np.array(self.kinematics())[:,-1]]).transpose()
        erro = math.sqrt((p[0] - xi[0][0])**2 + (p[1]-xi[1][0])**2 + (p[2]-xi[2][0])**2)
        th = self.th
        beta = 0.1
        inter = 0

        while(erro > 1e-10 and inter < 10000):
       
            inter = inter +1
            
            J = self.jacobian(th)
            
            Jinv = J.transpose()
            # Jinv = np.linalg.inv(J)
        #   Jinv = J.transpose()
            dxi = beta*np.add(np.array([p]).transpose(), -xi)
            dth = np.matmul(Jinv, dxi)[0] 
            th = th + dth
            xi = np.array([np.array(self.kinematics(th))[:,-1]]).transpose()
            erro = math.sqrt((p[0] - xi[0][0])**2 + (p[1]-xi[1][0])**2 + (p[2]-xi[2][0])**2)
        # print(inter)

        return th

        # return self.th  # substituir esta linha

    def dynamics(self, tau = np.array([[0], [0], [0]])):
        # implementar a dinâmica
        return np.array([[0, 0, 0]]).transpose()  # substituir esta linha

    def jacobian(self, th):
        j11 = - self.l2*sin(th[0][0])*cos(th[1][0]) - self.l3 *sin(th[0][0])* cos(th[1][0] + th[2][0])
        j12 = - self.l2*cos(th[0][0])*sin(th[1][0]) - self.l3 *cos(th[0][0])* sin(th[1][0] + th[2][0])
        j13 = - self.l3*cos(th[0][0])*sin(th[1][0] + th[2][0])
        j21 =   self.l2*cos(th[0][0])*cos(th[1][0]) + self.l3 *cos(th[0][0])* cos(th[1][0] + th[2][0])
        j22 = - self.l2*sin(th[0][0])*sin(th[1][0]) - self.l3 *sin(th[0][0])* sin(th[1][0] + th[2][0])
        j23 = - self.l3*sin(th[0][0])*sin(th[1][0] + th[2][0])
        j31 =   0
        j32 =   self.l2*cos(th[1][0]) + self.l3 * cos(th[1][0] + th[2][0])
        j33 =   self.l3 * cos(th[1][0] + th[2][0])


        

        return np.array([[j11,j12, j13],[j21,j22,j23], [j31,j32,j33]])