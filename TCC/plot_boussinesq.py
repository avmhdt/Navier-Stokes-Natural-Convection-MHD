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
resdir = workdir + "Results/Natural Convection Boussinesq/"

#%%
Nel = 80
Pr = 0.71
Ra = 1000000
dt = 1
#it = 6190

alpha_psi = 0.1#1.00
alpha_w = 0.01#1.00
alpha_T = 0.5#1.00

'''
omega = pd.read_csv('_'.join([resdir + "omega", str(Nel), str(Pr), str(Ra), "dt" + str(dt), str(it) + ".csv"]), header = None)
psi = pd.read_csv('_'.join([resdir + "psi", str(Nel), str(Pr), str(Ra), "dt" + str(dt), str(it) + ".csv"]), header = None)
T = pd.read_csv('_'.join([resdir + "T", str(Nel), str(Pr), str(Ra), "dt" + str(dt), str(it) + ".csv"]), header = None)

#%%

a = 0.
b = 1.
h = (b-a)/Nel

x = np.linspace(a, b, Nel+1)
x, y = np.meshgrid(x, x)

fig = plt.figure(1)
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(X = x, Y = y, Z = omega.values.T, rstride = 1, cstride = 1, cmap = 'viridis')
ax = fig.add_subplot(111)
ax.contour(x, y, omega.values.T, levels = 20, cmap = 'viridis')

fig = plt.figure(2)
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(X = x, Y = y, Z = psi.values.T, rstride = 1, cstride = 1, cmap = 'viridis')
ax = fig.add_subplot(111)
ax.contour(x, y, psi.values.T, levels = 20, cmap = 'viridis')

fig = plt.figure(3)
ax = fig.add_subplot(111)
cont = ax.contour(x, y, T.values.T, levels = 10, cmap = 'viridis')
ax.clabel(cont, inline = True, fontsize = 10)

'''
###########################################################################################
###########################################################################################
###########################################################################################
#%%

from os import listdir
from os.path import isfile, join

workdir = "C:/Users/Jose Carvalho/Documents/PROG/TCC/Trabalho 1 SNED Parte B/Iury MOD/"
resdir = workdir + "Results/Natural Convection Boussinesq/"
figdir =  workdir + 'Figures/Natural Convection Boussinesq/'
tabledir = workdir + 'Tables/Natural Convection Boussinesq/'

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

psifilename = '_'.join(['psi', str(Nel), str(Pr), str(Ra), '%.2f' % alpha_psi, '%.2f' % alpha_w, '%.2f' % alpha_T, 'dt' + str(which_h)])
omegafilename = '_'.join(['omega', str(Nel), str(Pr), str(Ra), '%.2f' % alpha_psi, '%.2f' % alpha_w, '%.2f' % alpha_T, 'dt' + str(which_h)])
Tfilename = '_'.join(['T', str(Nel), str(Pr), str(Ra), '%.2f' % alpha_psi, '%.2f' % alpha_w, '%.2f' % alpha_T, 'dt' + str(which_h)])
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
psimid = abs(psi[x[0, :] == 0.5, y[:, 0] == 0.5])
np.savetxt(tabledir + '_'.join(['psimid', str(Nel), str(Pr), str(Ra), '%.2f' % alpha_psi, '%.2f' % alpha_w, '%.2f' % alpha_T, 'dt' + str(which_h) + '.csv']),
           psimid, delimiter = ',', fmt = '%2.5f')

#%%
idxmaxpsi = np.where(abs(psi) == abs(psi).max())
rowmaxpsi = idxmaxpsi[0][0]
colmaxpsi = idxmaxpsi[1][0]

maxpsiarray = np.array([psi[rowmaxpsi, colmaxpsi], x[0, :][rowmaxpsi], y[:, 0][colmaxpsi]])
np.savetxt(tabledir + '_'.join(['maxpsi', str(Nel), str(Pr), str(Ra), '%.2f' % alpha_psi, '%.2f' % alpha_w, '%.2f' % alpha_T, 'dt' + str(which_h) + '.csv']),
           maxpsiarray, delimiter = ',', fmt = '%2.5f')

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

umid = u[x[0,:] == 0.5, :]
vmid = v[:, y[:,0] == 0.5]

umax = abs(umid).max()
vmax = abs(vmid).max()

yumax = y[:, 0][abs(umid)[0,:] == umax][0]
xvmax = x[0, :][abs(vmid)[:,0] == vmax][0]

uarray = np.array([umax, yumax])
varray = np.array([vmax, xvmax])

np.savetxt(tabledir + '_'.join(['umax', str(Nel), str(Pr), str(Ra), '%.2f' % alpha_psi, '%.2f' % alpha_w, '%.2f' % alpha_T, 'dt' + str(which_h) + '.csv']),
           uarray, delimiter = ',', fmt = '%2.5f')
np.savetxt(tabledir + '_'.join(['vmax', str(Nel), str(Pr), str(Ra), '%.2f' % alpha_psi, '%.2f' % alpha_w, '%.2f' % alpha_T, 'dt' + str(which_h) + '.csv']),
           varray, delimiter = ',', fmt = '%2.5f')

#%%
'''
dT = np.zeros(J)
for j in range(0, J):
    dT[j] = (T[2, j] - T[0, j])/(2.*h)

NuBar = 0.
for i in range(0, Nel):
    NuBar = NuBar + (dT[i] + dT[i+1])*h/2.

print(NuBar)

np.savetxt(tabledir + '_'.join(['NuBar', str(Nel), str(Pr), str(Ra), '%.2f' % alpha_psi, '%.2f' % alpha_w, '%.2f' % alpha_T, 'dt' + str(which_h) + '.csv']),
           np.array([NuBar]), delimiter = ',', fmt = '%2.5f')
'''
#%%
Q = np.zeros((I,J))
for i in range(0, I):
    for j in range(0, J):
        if i == 0:
            Q[i, j] = u[i, j]*T[i, j] - (T[i+2, j] - T[i, j])/(2.*h)
        elif i == I-1:
            Q[i, j] = u[i, j]*T[i, j] - (T[i, j] - T[i-2, j])/(2.*h)
        else:
            Q[i, j] = u[i, j]*T[i, j] - (T[i+1, j] - T[i-1, j])/(2.*h)
            
#Nux = np.zeros(I)
#for i in range(I):
#    Nux[i] = (Q[i, 0] + 4.*Q[i, y[:, 0] == 0.5] + Q[i, J-1])/6.

#NuBar = (Nux[0] + 4.*Nux[x[0, :] == 0.5] + Nux[I-1])/6.
#print(NuBar)

Nux = np.zeros(I)
for i in range(I):
    for j in range(J-1):
        Nux[i] = Nux[i] + (Q[i, j] + Q[i, j+1])*h/2.

NuBar = 0.
for i in range(I-1):
    NuBar = NuBar + (Nux[i] + Nux[i+1])*h/2.

NuMid = Nux[x[0, :] == 0.5][0]
Nu0 = Nux[0]

print(NuBar)
print(Nu0)
print(NuMid)
plt.plot(x[0,:], Nux)

NuArray = np.array([NuBar, NuMid, Nu0])

np.savetxt(tabledir + '_'.join(['Nu', str(Nel), str(Pr), str(Ra), '%.2f' % alpha_psi, '%.2f' % alpha_w, '%.2f' % alpha_T, 'dt' + str(which_h) + '.csv']),
           NuArray, delimiter = ',', fmt = '%2.5f')


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

