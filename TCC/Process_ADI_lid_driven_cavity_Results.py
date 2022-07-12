# -*- coding: utf-8 -*-
"""
Created on Wed May 18 08:21:43 2022

@author: Jose Carvalho
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join

#workdir = "C:/Users/Jose Carvalho/Documents/PROG/TCC/Trabalho 1 SNED Parte B/Iury MOD/"
workdir = "C:/Users/Jose Carvalho/source/repos/Navier_Stokes_Thomas_Parallel/Navier_Stokes_Thomas_Parallel/"
resdir = workdir + "Results/Navier-Stokes/"
figdir =  workdir + 'Figures/Navier-Stokes/'
tabledir = workdir + 'Tables/Navier-Stokes/'
maxstreamdir = workdir + 'Max Stream Function/Navier-Stokes/'

#%%
Nel = 64
Re = 100
which_h = 1

#%%
if which_h == 1:
    dtString = 'h2'
    dtStringPlot = 'h^2'
elif which_h == 2:
    dtString = 'h2s2'
    dtStringPlot = 'h^2/2'
elif which_h == 3:
    dtString = 'h2s4'
    dtStringPlot = 'h^2/4'
elif which_h == 4:
    dtString = 'h2s10'
    dtStringPlot = 'h^2/10'
elif which_h == 5:
    dtString = 'h3'
    dtStringPlot = 'h^3'
elif which_h == 6:
    dtString = 'h4'
    dtStringPlot = 'h^4'
else: print('Invalid which_h.\n')

#%%
files = listdir(resdir)
files = [file for file in files if isfile(join(resdir, file))]

psifilename = '_'.join(['psi', str(Nel), str(Re), 'dt' + str(which_h)])
omegafilename = '_'.join(['omega', str(Nel), str(Re), 'dt' + str(which_h)])
psifile =  [file for file in files if file.find(psifilename) != -1][0]
omegafile =  [file for file in files if file.find(omegafilename) != -1][0]

psi = pd.read_csv(resdir + psifile, header = None).values
omega = pd.read_csv(resdir + omegafile, header = None).values
    
#%%

a = 0.0
b = 1.0
h = (b-a)/Nel

x = [a + i*h for i in range(0, Nel + 1)]
y = [a + j*h for j in range(0, Nel + 1)]

plt.figure(1)
plt.title(r'Nel = ' + str(Nel) + ',   $\Delta t = ' + dtStringPlot + '$,   Re = ' + str(Re))
plt.suptitle('Função Corrente ' + r'$\psi$')
plt.contour(x, y, psi.T, [psi.min() - psi.min()/50.0*t for t in range(0, 50)] + [psi.max()/10.0*t for t in range(0, 11)])
plt.gcf().set_size_inches(7, 7)
plt.savefig(figdir + psifile[:-3] + 'pdf', bbox_inches = 'tight')
#plt.close()

plt.figure(2)
plt.title(r'Nel = ' + str(Nel) + ',   $\Delta t = ' + dtStringPlot + '$,   Re = ' + str(Re))
plt.suptitle('Vorticidade ' + r'$\omega$')
plt.contour(x, y, omega.T, 500)
plt.gcf().set_size_inches(7, 7)
plt.savefig(figdir + omegafile[:-3] + 'pdf', bbox_inches = 'tight')
#plt.close()

exit()
##############################################################################
#%%
NelGhia = 128
hghia = (b-a)/NelGhia

xghia = [a + i*hghia for i in range(0, NelGhia + 1)]
yghia = [a + j*hghia for j in range(0, NelGhia + 1)]

idxGhiaY = [0, 7, 8, 9, 13, 22, 36, 58, 64, 79, 94, 109, 122, 123, 124, 125, 128];
idxGhiaX = [0, 8, 9, 10, 12, 20, 29, 30, 64, 103, 110, 116, 121, 122, 123, 124, 128];

xghia = [xghia[i] for i in idxGhiaX];
yghia = [yghia[i] for i in idxGhiaY];

ReValuesGhia = [100, 400, 1000, 3200, 5000, 7500, 10000];

ughia = np.array([[0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000],
        [-0.03717, -0.08186, -0.18109, -0.32407, -0.41165, -0.43154, -0.42735],
        [-0.04192, -0.09266, -0.20196, -0.35344, -0.42901, -0.43590, -0.42537],
        [-0.04775, -0.10338, -0.22220, -0.37827, -0.43643, -0.43025, -0.41657],
        [-0.06434, -0.14612, -0.29730, -0.41933, -0.40435, -0.38324, -0.38000],
        [-0.10150, -0.24299, -0.38289, -0.34323, -0.33050, -0.32393, -0.32709],
        [-0.15662, -0.32726, -0.27805, -0.24427, -0.22855, -0.23176, -0.23186],
        [-0.21090, -0.17119, -0.10648, -0.086636, -0.07404, -0.07503, -0.07540],
        [-0.20581, -0.11477, -0.06080, -0.04272, -0.03039, -0.03800, -0.03111],
        [-0.13641, 0.02135, 0.05702, 0.07156, 0.08183, 0.08342, 0.08344],
        [ 0.00332, 0.16256, 0.18719, 0.19791, 0.20087, 0.20591, 0.20673],
        [ 0.23151, 0.29093, 0.33304, 0.34682, 0.33556, 0.34228, 0.34635],
        [ 0.68717, 0.55892, 0.46604, 0.46101, 0.46036, 0.47167, 0.47804],
        [ 0.73722, 0.61756, 0.51117, 0.46547, 0.45992, 0.47323, 0.48070],
        [ 0.78871, 0.68439, 0.57492, 0.48296, 0.46120, 0.47048, 0.47783],
        [ 0.84123, 0.75837, 0.65928, 0.53236, 0.48223, 0.47244, 0.47221],
        [ 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000]])

vghia = np.array([[0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000],
         [0.09233, 0.18360, 0.27485, 0.39560, 0.42447, 0.43979, 0.43983],
         [0.10091, 0.19713, 0.29012, 0.40917, 0.43329, 0.44030, 0.43733],
         [0.10890, 0.20920, 0.30353, 0.41906, 0.43648, 0.43564, 0.43124],
         [0.12317, 0.22965, 0.32627, 0.42768, 0.42951, 0.41824, 0.41487],
         [0.16077, 0.28124, 0.37095, 0.37119, 0.35368, 0.35060, 0.35070],
         [0.17507, 0.30203, 0.33075, 0.29030, 0.28066, 0.28117, 0.28003],
         [0.17527, 0.30174, 0.32235, 0.28188, 0.27280, 0.27348, 0.27224],
         [0.05454, 0.05186, 0.02526, 0.00999, 0.00945, 0.00824, 0.00831],
         [-0.24533, -0.38598, -0.31966, -0.31184, -0.30018, -0.30448, -0.30719],
         [-0.22445, -0.44993, -0.42665, -0.37401, -0.36214, -0.36213, -0.36737],
         [-0.16914, -0.23827, -0.51550, -0.44307, -0.41442, -0.41050, -0.41496],
         [-0.10313, -0.22847, -0.39188, -0.54053, -0.52876, -0.48590, -0.45863],
         [-0.08864, -0.19254, -0.33714, -0.52357, -0.55408, -0.52347, -0.49099],
         [-0.07391, -0.15663, -0.27669, -0.47425, -0.55069, -0.55216, -0.52987],
         [-0.05906, -0.12146, -0.21388, -0.39017, -0.49774, -0.53858, -0.54302],
         [0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000]])

#%%

I = Nel+1
J = I

u = np.empty((I, J))
v = np.empty((I, J))

u[0] = 0;
u[I-1] = 0;
u[:,0] = 0;
u[:,J-1] = 1;

v[0] = 0;
v[I-1] = 0;
v[:,0] = 0;
v[:,J-1] = 0;

for i in range(1, I-1):
    for j in range(1, J-1):
        u[i, j] = (psi[i][j+1] - psi[i][j-1])/(2*h)
        v[i, j] = -(psi[i+1][j] - psi[i-1][j])/(2*h)

umid = u[x.index(0.5), :]
vmid = v[:, y.index(0.5)]

#%%
plt.figure(3)
plt.plot(y, umid)
plt.xlabel('y, x')
plt.ylabel('u, v')
plt.grid()
plt.title('Nel = ' + str(Nel) + ',   ' + r'$\Delta t = ' + dtStringPlot + '$,   Re = ' + str(Re))
plt.suptitle(r'$u$ em $x = 0.5$, $v$ em $y = 0.5$')
plt.plot(yghia, ughia[:,ReValuesGhia.index(Re)], 'or', mfc = 'none')

plt.plot(x, vmid)
plt.plot(xghia, vghia[:,ReValuesGhia.index(Re)], 'xb', mfc = 'none')

plt.gcf().set_size_inches(7, 7)
plt.legend(['u - ADI', 'u - Ghia et. al', 'v - ADI', 'v - Ghia et. al'])
plt.savefig(figdir + 'uv_at_xy_0.5_' + '_'.join([str(Nel), str(Re), 'dt' + str(which_h)]) + '.pdf')
plt.close()

'''
plt.figure(4)
plt.plot(x, vmid)
plt.xlabel('x')
plt.ylabel('v')
plt.grid()
plt.title('Nel = ' + str(Nel) + ',   ' + r'$\Delta t = ' + dtStringPlot + '$,   Re = ' + str(Re))
plt.suptitle(r'$v$ em $y = 0.5$')
plt.plot(xghia, vghia[:,ReValuesGhia.index(Re)], 'or', mfc = 'none')
plt.gcf().set_size_inches(7, 7)
plt.legend(['ADI', 'Ghia et. al'])
plt.savefig(figdir + 'v_at_y_0.5_'  + '_'.join([str(Nel), str(Re), 'dt' + str(which_h)]) + '.pdf')
#plt.close()
'''
#%%
np.savetxt(tabledir + '_'.join(['umid', str(Nel), str(Re), 'dt' + str(which_h) + '.csv']), umid, delimiter = ',')
np.savetxt(tabledir + '_'.join(['vmid', str(Nel), str(Re), 'dt' + str(which_h) + '.csv']), vmid, delimiter = ',')

xgidx = [x.index(xg) for xg in xghia]
ygidx = [y.index(yg) for yg in yghia]

np.savetxt(tabledir + '_'.join(['umid_on_ghia', str(Nel), str(Re), 'dt' + str(which_h) + '.csv']), umid[ygidx], delimiter = ',', fmt = '%2.5f')
np.savetxt(tabledir + '_'.join(['vmid_on_ghia', str(Nel), str(Re), 'dt' + str(which_h) + '.csv']), vmid[xgidx], delimiter = ',', fmt = '%2.5f')

#%%
idxmaxpsi = np.where(abs(psi) == abs(psi).max())
rowmaxpsi = idxmaxpsi[0][0]
colmaxpsi = idxmaxpsi[1][0]

maxarray = np.array([psi[rowmaxpsi, colmaxpsi], omega[rowmaxpsi, colmaxpsi], x[rowmaxpsi], y[colmaxpsi]])
np.savetxt(maxstreamdir + '_'.join(['maxpsi', str(Nel), str(Re), 'dt' + str(which_h) + '.csv']), maxarray, delimiter = ',')

