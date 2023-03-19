import numpy as np

class Rotor2:
    '''Rotor2 object with default test rotor for hw3'''
    def __init__(self):
        self._Nb = 1  # number of blades
        self._R = 30  # ft
        self._vtip = 780  # ft/sec

        # Coefficient properties related to linear radial variations of airfoil cl_slope, theta twist, and chord
        self._AV = 4*np.pi
        self._cl_slope_75 = 5*np.pi

        self._theta_tw = -np.deg2rad(8)
        self._theta_75 = np.deg2rad(5)

        self._TR = 2
        self._chord_75 = 2  # ft

    @property
    def Nb(self):
        return self._Nb

    @property
    def R(self):
        return self._R

    @property
    def vtip(self):
        return self._vtip
    
    @property
    def AV(self):
        return self._AV
    
    @property
    def cl_slope_75(self):
        return self._cl_slope_75
    
    @property
    def theta_tw(self):
        return self._theta_tw
    
    @property
    def theta_75(self):
        return self._theta_75
    
    @property
    def TR(self):
        return self._TR
    
    @property
    def chord_75(self):
        return self._chord_75

    @Nb.setter
    def Nb(self, val):
        self._Nb = val

    @R.setter
    def R(self, val):
        self._R = val

    @vtip.setter
    def vtip(self, val):
        self._vtip = val

    @AV.setter
    def AV(self, val):
        self._AV = val

    @cl_slope_75.setter
    def cl_slope_75(self, val):
        self._cl_slope_75 = val

    @theta_tw.setter
    def theta_tw(self, val):
        self._theta_tw = val

    @theta_75.setter
    def theta_75(self, val):
        self._theta_75 = val

    @TR.setter
    def TR(self, val):
        self._TR = val

    @chord_75.setter
    def chord_75(self, val):
        self._chord_75 = val
