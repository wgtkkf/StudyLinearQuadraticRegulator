# Tanh-Sinh quadrature class
# Coded by Takuro Tokunaga
# Last updated:
# March 21st, 2026

import numpy as np

class TanhSinh:
    # integral parameters
    SN = 1                   # n
    SH = np.power(10.0,-5)   # h, -6
    CONV = np.power(10.0,-4) # convergence criteria, -8
    DIFF = 1
    PI = np.pi

    def __init__(self):
        pass

    def PHIP(self, sn):
        y = np.tanh(0.5*self.PI*np.sinh(sn*self.SH))
        return y

    def DPHIP(self, sn):
        y = 0.5*self.PI*np.cosh(sn*self.SH)/np.power(np.cosh(0.5*self.PI*np.sinh(sn*self.SH)),2)
        return y

    def PHIE(self, sn):
        y = np.exp(0.5*self.PI*np.sinh(sn*self.SH))
        return y

    def DPHIE(self, sn):
        y = 0.5*self.PI*np.cosh(sn*self.SH)*np.exp(0.5*self.PI*np.sinh(sn*self.SH))
        return y        