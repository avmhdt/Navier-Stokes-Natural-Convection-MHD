# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 17:40:23 2021

@author: Jose Carvalho
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

workdir = "C:/Users/Jose Carvalho/Documents/PROG/TCC/Trabalho 1 SNED Parte B/Iury MOD/"
resdir = workdir + "Results/MHD/"

#%%
Nel = 128
Pr = 1.00
Ra = 1000000
Ha = 50
dt = 1
#it = 6190

N = 0.5

alpha_psi = 0.1#1.00
alpha_w = 0.01#1.00
alpha_T = 0.5#1.00

###########################################################################################
###########################################################################################
###########################################################################################
#%%

from os import listdir
from os.path import isfile, join

workdir = "C:/Users/Jose Carvalho/Documents/PROG/TCC/Trabalho 1 SNED Parte B/Iury MOD/"
resdir = workdir + "Results/MHD/"
figdir =  workdir + 'Figures/MHD/'
tabledir = workdir + 'Tables/MHD/'

#%%

which_h = dt

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

psifilename = '_'.join(['psi', str(Nel), '%.2f' % Pr, str(Ra), str(Ha), str(N), '%.2f' % alpha_psi, '%.2f' % alpha_w, '%.2f' % alpha_T, 'dt' + str(which_h)])
omegafilename = '_'.join(['omega', str(Nel), '%.2f' % Pr, str(Ra), str(Ha), str(N), '%.2f' % alpha_psi, '%.2f' % alpha_w, '%.2f' % alpha_T, 'dt' + str(which_h)])
Tfilename = '_'.join(['T', str(Nel), '%.2f' % Pr, str(Ra), str(Ha), str(N), '%.2f' % alpha_psi, '%.2f' % alpha_w, '%.2f' % alpha_T, 'dt' + str(which_h)])
psifile =  [file for file in files if file.find(psifilename) != -1][0]
omegafile =  [file for file in files if file.find(omegafilename) != -1][0]
Tfile =  [file for file in files if file.find(Tfilename) != -1][0]

psi = pd.read_csv(resdir + psifile, header = None).values
omega = pd.read_csv(resdir + omegafile, header = None).values
T = pd.read_csv(resdir + Tfile, header = None).values

#%%
a = 0.
b = 1.
h = (b-a)/Nel

x = np.linspace(a, b, Nel+1)
x, y = np.meshgrid(x, x)

plt.figure(1)
plt.title(r'Nel = ' + str(Nel) + ',   $\Delta t = ' + dtStringPlot + '$,   Pr = ' + str(Pr) + ',    Ra = ' + str(Ra))
plt.suptitle('Função Corrente ' + r'$\psi$')
plt.contour(x, y, psi.T, levels = 20, cmap = 'viridis')
plt.gcf().set_size_inches(7, 7)
plt.savefig(figdir + psifile[:-3] + 'pdf', bbox_inches = 'tight')
plt.close()

plt.figure(2)
plt.title(r'Nel = ' + str(Nel) + ',   $\Delta t = ' + dtStringPlot + '$,   Pr = ' + str(Pr) + ',    Ra = ' + str(Ra))
plt.suptitle('Vorticidade ' + r'$\omega$')
plt.contour(x, y, omega.T, levels = 20, cmap = 'viridis')
plt.gcf().set_size_inches(7, 7)
plt.savefig(figdir + omegafile[:-3] + 'pdf', bbox_inches = 'tight')
plt.close()

plt.figure(3)
plt.title(r'Nel = ' + str(Nel) + ',   $\Delta t = ' + dtStringPlot + '$,   Pr = ' + str(Pr) + ',    Ra = ' + str(Ra))
plt.suptitle('Temperatura ' + r'$T$')
cont = plt.contour(x, y, T.T, levels = 10, cmap = 'viridis')
plt.clabel(cont, inline = True, fontsize = 10)
plt.gcf().set_size_inches(7, 7)
plt.savefig(figdir + Tfile[:-3] + 'pdf', bbox_inches = 'tight')
plt.close()

#%%
I = J = Nel+1

#%%
u = np.empty((I, J))
v = np.empty((I, J))

u[0] = 0;
u[I-1] = 0;
u[:,0] = 0;
u[:,J-1] = 0;

v[0] = 0;
v[I-1] = 0;
v[:,0] = 0;
v[:,J-1] = 0;

for i in range(1, I-1):
    for j in range(1, J-1):
        u[i, j] = (psi[i][j+1] - psi[i][j-1])/(2.*h)
        v[i, j] = -(psi[i+1][j] - psi[i-1][j])/(2.*h)

#%%
plt.figure(4)
plt.plot(v[x[0, :] == 0.25][0], y[:, 0])
plt.plot(v[x[0, :] == 0.5][0], y[:, 0])
plt.plot(v[x[0, :] == 0.75][0], y[:, 0])
plt.grid()
plt.xlabel('v')
plt.ylabel('y')
plt.gcf().set_size_inches(7, 7)
plt.legend(['x = 0.25', 'x = 0.5', 'x = 0.75'])
plt.savefig(figdir + '_'.join(['v_at_x', str(Nel), '%.2f' % Pr, str(Ra), str(Ha), str(N), '%.2f' % alpha_psi, '%.2f' % alpha_w, '%.2f' % alpha_T, 'dt' + str(which_h)]) + '.pdf', bbox_inches = 'tight')
plt.close()

plt.figure(5)
plt.plot(T[x[0, :] == 0.25][0], y[:, 0])
plt.plot(T[x[0, :] == 0.5][0], y[:, 0])
plt.plot(T[x[0, :] == 0.75][0], y[:, 0])
plt.grid()
plt.xlabel('T')
plt.ylabel('y')
plt.gcf().set_size_inches(7, 7)
plt.legend(['x = 0.25', 'x = 0.5', 'x = 0.75'])
plt.savefig(figdir + '_'.join(['T_at_x', str(Nel), '%.2f' % Pr, str(Ra), str(Ha), str(N), '%.2f' % alpha_psi, '%.2f' % alpha_w, '%.2f' % alpha_T, 'dt' + str(which_h)]) + '.pdf', bbox_inches = 'tight')
plt.close()


#%%
Q = np.zeros((I,J))
for i in range(0, I):
    for j in range(0, J):
        if j == 0:
            Q[i, j] = v[i, j]*T[i, j] - (T[i, j+2] - T[i, j])/(2.*h)
        elif j == J-1:
            Q[i, j] = v[i, j]*T[i, j] - (T[i, j] - T[i, j-2])/(2.*h)
        else:
            Q[i, j] = v[i, j]*T[i, j] - (T[i, j+1] - T[i, j-1])/(2.*h)

Nux = np.zeros(J)
for j in range(J):
    for i in range(I-1):
        Nux[j] = Nux[j] + (Q[i, j] + Q[i+1, j])*h/2.

NuBar = 0.
for i in range(I-1):
    NuBar = NuBar + (Nux[i] + Nux[i+1])*h/2.

NuMid = Nux[x[0, :] == 0.5][0]
Nu0 = Nux[0]

print(NuBar)
print(Nu0)
print(NuMid)
#plt.plot(x[0,:], Nux)

NuArray = np.array([NuBar, NuMid, Nu0])

'''
Nux = np.zeros((I, J))
for j in range(J):
    for i in range(I):
        if j == 0:
            Nux[i, j] = -(T[i, j+2] - T[i, j])/(2.*h)
        elif j == J-1:
            Nux[i, j] = -(T[i, j] - T[i, j-2])/(2.*h)
        else:
            Nux[i, j] = -(T[i, j+1] - T[i, j-1])/(2.*h)

NuAvg = np.zeros(I)
for i in range(I):
    for j in range(J-1):
        NuAvg[i] = NuAvg[i] + (Nux[i, j] + Nux[i, j+1])*h/2.

NuBar = 0.
for i in range(I-1):
    NuBar = NuBar + (NuAvg[i] + NuAvg[i+1])*h/2.


print(NuBar)
plt.figure(4)
plt.plot(x[0, :], NuAvg)
plt.close()
'''
'''
np.savetxt(tabledir + '_'.join(['Nu', str(Nel), str(Pr), str(Ra), '%.2f' % alpha_psi, '%.2f' % alpha_w, '%.2f' % alpha_T, 'dt' + str(which_h) + '.csv']),
           NuArray, delimiter = ',', fmt = '%2.5f')
'''

#%%
'''
#dT = np.zeros(J)
#for j in range(J):
#    dT[j] = (T[2, j] - 2.*T[1, j] + T[0, j])/(h**2.)

#Nu0 = (dT[0] + 4.*dT[y[:, 0] == 0.5][0] + dT[J-1])/6.

Nu0 = 0.
for j in range(J-1):
    Nu0 = Nu0 + (dT[j] + dT[j+1])*h/2.
    
print(Nu0)
'''

