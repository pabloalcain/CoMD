import numpy as np
import pylab as pl
import scipy as sp
from scipy.special import erf
from scipy.integrate import dblquad
sigma = 1.0


def function(z0, r, z):
  f = r/np.sqrt(z**2 + r**2) * np.exp(-(r**2 + (z - z0)**2)/(2*sigma**2))
  f *= erf(np.sqrt(r**2 + z**2)/(np.sqrt(2)*sigma))
  return f

def integrand(z0):
  """
  Return the integrand for z0
  """
  return lambda r, z: function(z0, r, z)


def Vinteg(r):
  V = np.zeros_like(r)
  err = np.zeros_like(r)
  for i, _r in enumerate(r):
    print i
    res = dblquad(integrand(_r), -np.inf, np.inf,
                   lambda _: 0, lambda _: np.inf)
    V[i] = res[0]
    print V[i]
    err[i] = res[1]
  return V/(sigma**3 * np.sqrt(2*np.pi)), err/(sigma**3 * np.sqrt(2*np.pi))

r = np.linspace(0, 10, 1000)
V = {}
V['QMD'] = 1/r*erf(r/(np.sqrt(2)*sigma))
V['Point particles'] = 1/r
res = Vinteg(r)
V['Integration'] = res[0]
err = res[1]
fig, ax = pl.subplots()
for key, val in zip(V.keys(), V.values()):
  ax.loglog(r, val, label=key)
ax.legend()
ax.set_ylabel('Potential')
ax.set_xlabel('Distance')
pl.show()
