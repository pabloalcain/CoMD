""" This script plots the interaction in QMD/CoMD model """
import numpy as np
import pylab as pl
from scipy import special
t0 = -356
t3 = 303
d0 = 0.16
u = 7.0/6.0
a = 32
sig_r = 1.3
cs = -0.33
e = np.sqrt(1.44)

r = np.linspace(0, 10, 1000)
dij = 1.0/(2*np.sqrt(np.pi)*sig_r)**3 * pl.exp(-(r/(2*sig_r))**2)


Vvol = t0/2 * dij/d0
V3 = 2*t3/(u + 1) * (dij/d0)**u
Vsym = a/2 * (dij/d0)
Vsurf = cs * (dij/d0) * (r**2-2*sig_r**2)/(4*sig_r**4)
Vcoul = e**2/r*special.erf(r/(2*sig_r))

pl.figure()
pl.plot(r, Vvol, label='Volume')
pl.plot(r, V3, label='3-body')
pl.plot(r, Vsym, label='Symmetry')
pl.plot(r, Vsurf, label='Surface')
pl.plot(r, Vcoul, label='Coulomb')
pl.legend()

pl.figure()
pl.plot(r, Vsurf + V3 + Vsym + Vvol + Vcoul, label='Proton-proton')
pl.plot(r, Vsurf + V3 - Vsym + Vvol, label='Proton-neutron')
pl.plot(r, Vsurf + V3 + Vsym + Vvol, label='Neutron-neutron')
pl.legend()

pl.show()

