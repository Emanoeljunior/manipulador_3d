from math import sin, cos, pi
import math
from re import A
import numpy as np
import math


class model(object):
    def __init__(self, th=[0, 0, 0], dt=0.0005):

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

        self.M = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
        self.V = np.array([[0], [0], [0]])
        self.G = np.array([[0], [0], [0]])

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
        # Cinemática
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
            dxi = beta*np.add(np.array([p]).transpose(), -xi)
            dth = np.matmul(Jinv, dxi)[0] 
            th = th + dth
            xi = np.array([np.array(self.kinematics(th))[:,-1]]).transpose()
            erro = math.sqrt((p[0] - xi[0][0])**2 + (p[1]-xi[1][0])**2 + (p[2]-xi[2][0])**2)

        return th


    def dynamics(self, tau = np.array([[0], [0], [0]])):
        print('th ',self.th)
        th1 = self.th[0][0]
        th2 = self.th[1][0]
        th3 = self.th[2][0]

        g = self.g

        d = self.l1
        e = self.l2 
        f = self.l3 
        
        m_1 = self.m1
        m_2 = self.m2
        m_3 = self.m3


        m00 = e**2*m_2*math.cos(th2)**2 + m_3*(e*math.cos(th2) + f*math.cos(th2 + th3))**2
        m01 = 0
        m02 = 0
        m10 = 0
        m11 = e**2*m_2 + m_3*(e**2 + 2*e*f*math.cos(th3) + f**2)
        m12 = 1.0*f*m_3*(e*math.cos(th3) + f)
        m20 = 0
        m21 = 1.0*f*m_3*(e*math.cos(th3) + f)
        m22 = 1.0*f**2*m_3



        M = np.array([[m00, m01, m02], [m10, m11, m12], [m20, m21, m22]])
        # M = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
        self.M = M


        v00 = -2*f*m_3*th3**2*(e*math.cos(th2) + f*math.cos(th2 + th3))*math.sin(th2 + th3) - th2**2*(e**2*m_2*math.sin(2*th2) + e**2*m_3*math.sin(2*th2) + 2*e*f*m_3*math.sin(2*th2 + th3) + f**2*m_3*math.sin(2*th2 + 2*th3))
        v01 = 0
        v02 = 0
        v10 = 0
        v11 = -2*e*f*m_3*th3**2*math.sin(th3)
        v12 = -1.0*e*f*m_3*th3**2*math.sin(th3)
        v20 = 0
        v21 = e*f*m_3*(1.0*th2**2 - 0.5*th3**2)*math.sin(th3)
        v22 = 0.5*e*f*m_3*th2**2*math.sin(th3)

        V = np.array([[v00 + v01 + v02], [v10 + v11 + v12], [v20 + v21 + v22]])
        self.V = V
        # V = np.array([[0], [0], [0]])

        
        g0 = 0
        g1 = g*(e*m_2*math.cos(th2) + e*m_3*math.cos(th2) + f*m_3*math.cos(th2 + th3))
        g2 = f*g*m_3*math.cos(th2 + th3)

        G = np.array([[g0], [g1], [g2]])
        self.G = G
        # G = np.array([[0], [0], [0]])

        return np.array(np.matmul(np.linalg.inv(M), np.add(tau, -np.add(V, G))))
        # Dinâmica
        # return np.array([[0, 0, 0]]).transpose()  # substituir esta linha

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