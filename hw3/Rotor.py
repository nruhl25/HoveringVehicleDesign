import numpy as np

class Rotor:
    '''Rotor object with default test rotor for hw3'''
    def __init__(self):
        self._Nb = 1  # number of blades
        self._c = 2.0  # Chord length, ft
        self._R = 30  # ft
        self._vtip = 650  # ft/sec
        self._cl_slope = 2*np.pi
        self._sigma = self._Nb*self._c/(np.pi*self._R)

    @property
    def Nb(self):
        return self._Nb
    
    @property
    def c(self):
        return self._c

    @property
    def R(self):
        self._R

    @property
    def vtip(self):
        return self._vtip

    @property
    def cl_slope(self):
        return self._cl_slope

    @property
    def sigma(self):
        return self._Nb*self._c/(np.pi*self._R)

    @Nb.setter
    def Nb(self, val):
        self._Nb = val
    
    @c.setter
    def c(self, val):
        self._c = val

    @R.setter
    def R(self, val):
        self._R = val

    @vtip.setter
    def vtip(self, val):
        self._vtip = val

    @cl_slope.setter
    def cl_slope(self, val):
        self._cl_slope = val