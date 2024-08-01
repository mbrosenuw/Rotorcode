import numpy as np
from getPES import getPES
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.constants import physical_constants as pc

cmtoh = 100/pc['hartree-inverse meter relationship'][0]

aitoi = pc['atomic mass constant'][0] * pc['Angstrom star'][0]**2
jtocm = 1/(100*pc['inverse meter-joule relationship'][0])

hbar = pc['reduced Planck constant'][0]
def interPE(grid, PE, points):
    f = interp1d(grid, PE, kind='linear')
    x_new = np.linspace(grid.min(), grid.max(), points)
    # x_new = np.linspace(-115, 115, points)
    y_new = f(x_new)
    return x_new, y_new

def T(grid,A):
    n = len(grid)
    mat = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            coeff = A * (-1)**(i-j)
            if i == j:
                elem = (n*(n+1))/3
            else:
                num = np.cos(np.pi * (i-j)/(2*n + 1))
                denom = 2 * (np.sin(np.pi * (i-j)/ (2*n + 1)))**2
                elem = num/denom
            mat[i,j] = coeff * elem * cmtoh
    return mat

def CM_kinE(grid, mu):
    T = np.zeros((len(grid), len(grid)))
    for i in range(len(grid)):
       for j in range(len(grid)):
            coeff = mu * (-1)**(i-j)
            if i == j:
                val = np.pi**2/3
            else:
                val = 2/(i-j)**2
            T[i][j] = coeff*val
    return T

def V(grid, PE, units = True):
    n = len(grid)
    mat = np.zeros((n, n))

    for i in range(n):
        mat[i,i] = PE[i]
    if units:
        mat = mat * cmtoh
    return mat

def H(TE,VE):
    H = TE + VE
    energies, wfns = np.linalg.eigh(H)
    return energies, wfns

Iwater = aitoi * (1 * 0.98055084**2) * 2
Awater = jtocm * hbar ** 2 / (2*Iwater)
Iio = aitoi * (127 * 3.51321647**2)
Aio = jtocm * hbar ** 2 / (2*Iio)
A = Awater+Aio

grid, PE = getPES('iodsolidscan.csv')
grid, PE = interPE(grid, PE, 1000)
grid = (grid + 180)/360 * 2*np.pi
kinE = T(grid, A)
potE = V(grid, PE, False)
energies, wfns = H(kinE,potE)
print(energies[0:4]/cmtoh)
idx = 0
grid = grid / (2*np.pi)*360-180
for idx, wfn in enumerate((wfns.T)[5:8]):
    plt.plot(grid, (wfn**2 + energies[idx])/cmtoh)
plt.plot(grid, np.diag(potE)/cmtoh)
plt.show()

print('1st transition: ', (energies[1]-energies[0])/cmtoh)

