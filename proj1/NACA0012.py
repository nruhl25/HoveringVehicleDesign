# airfoil drag coefficient data for NACA0012 airfoil
# (Re=1.8x10^6, from NACA technical note 3361)

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

cd_list = np.array([0.01461567,
           0.0133942,
           0.01642706,
           0.01306375,
           0.03974506])

alphas_drag = np.array([-2.089976895,
               1.992422036,
               5.944140922,
               10.02673838,
               14.36712245])

alphas = np.deg2rad(alphas_drag)

def cd(alpha, cd0f, d1, d2):
    return cd0f + d1*alpha + d2*alpha**2

popt, pcov = curve_fit(cd, xdata=alphas, ydata=cd_list)
popt_deg, pcov_deg = curve_fit(cd, xdata=np.array(alphas_drag), ydata=cd_list)

print(popt)

plt.title("NACA0012 Airfoil Data (Re=1.8x10^6, from NACA technical note 3361)")
a_deg = np.linspace(alphas_drag[0], alphas_drag[-1], 100)
a = np.deg2rad(a_deg)
plt.plot(a_deg,cd(a_deg,*popt_deg), label="Quadratic Fit")
plt.plot(alphas_drag, cd_list, '.',label="Data")
plt.hlines(xmin=min(a_deg),xmax=max(a_deg),y=0.01,label=r"$c_{d,0}=0.01$ (hw3 assumption)", color='r')
plt.ylabel(f"Section drag coefficient, $c_d$")
plt.xlabel("Section angle of attack (deg)")
plt.grid()
plt.legend()
plt.show()
