import numpy as np
from math import cos, sin
import math

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
        return self.control_3(mdl)
        # implementar controlador
        # return np.array([[0], [0], [0]])  # substituir esta linha

    def control_1(self,mdl):
        #Controlador 100
        kp = 200
        kd = 100
        ki = 100

        

        eOld = self.e
        self.e = self.thRef - mdl.th
        self.ed = (self.e -eOld)/ (2*mdl.dt)
        self.ei += self.e*mdl.dt

        u = kp*self.e + ki*self.ei + kd*self.ed
        
        return np.array([[u[0][0]],[u[1][0]],[u[2][0]]])

    def control_2(self,mdl):
        #Controlador PID com compensação de  gravidade
            kp = 300
            kd = 100
            ki = 100

          
            eOld = self.e
            self.e = self.thRef - mdl.th
            self.ed = (self.e -eOld)/ (2*mdl.dt)
            self.ei += self.e*mdl.dt

            u = kp*self.e + ki*self.ei + kd*self.ed + mdl.G
            
            return np.array([[u[0][0]],[u[1][0]],[u[2][0]]])
    
    def control_3(self,mdl):
        #Controlador PD com compensação de torque e gravidade
            kp = 200
            kd = 100
            ki = 0

            
            eOld = self.e
            self.e = self.thRef - mdl.th
            self.ed = (self.e -eOld)/ (2*mdl.dt)
            self.ei += self.e*mdl.dt

            

            u = np.matmul(mdl.M,kp*self.e + ki*self.ei + kd*self.ed) + mdl.V + mdl.G
            
            return np.array([[u[0][0]],[u[1][0]],[u[2][0]]])
    
    