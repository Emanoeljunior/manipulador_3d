import numpy as np
from math import cos, sin

class control (object):
    def __init__(self):

        self.e = np.array([[0.0, 0.0 , 0.0]]).transpose()
        self.ed = np.array([[0.0, 0.0, 0.0]]).transpose()
        self.ei = np.array([[0.0, 0.0, 0.0]]).transpose()

        self.thRef = np.array([[0.0, 0.0, 0.0]]).transpose()

        return

    def control(self, mdl, ref=np.array([[0.0, 0.0, 0.0]]).transpose()):
        if ref is not None: self.thRef = ref
        # return self.thRef
        return self.control_1(mdl)
        # implementar controlador
        # return np.array([[0], [0], [0]])  # substituir esta linha

    def control_1(self,mdl):
        kp = 1000
        kd = 200
        ki = 1000

        
        eOld = self.e
        self.e = self.thRef - mdl.th
        self.ed = (self.e -eOld)/ (2*mdl.dt)
        self.ei += self.e*mdl.dt

        u = kp*self.e + ki*self.ei + kd*self.ed
        print("u")
        print(u)
        
        return np.array([[u[0][0]],[u[1][0]],[u[2][0]]])