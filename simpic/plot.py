import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv('phi.dat',sep='\s+',header=None)
data = pd.DataFrame(data)

x = data[1]
y = data[2]
plt.plot(x, y,'r.',label='all time steps')

plt.xlabel('x[m]')
plt.ylabel('phi[V]')
plt.title('SIMPIC version ***: Potential calculation')
plt.legend()

plt.savefig("plot_phi.png")

plt.clf()

data2 = pd.read_csv('E.dat',sep='\s+',header=None)
data2 = pd.DataFrame(data2)

x2 = data2[1]
y2 = data2[2]
plt.plot(x2, y2,'b.',label='all time steps')

plt.xlabel('x[m]')
plt.ylabel('E[V/m]')
plt.title('SIMPIC version ***: Electric field calculation')
plt.legend()

plt.savefig("plot_E.png")

plt.clf()

data3 = pd.read_csv('nt.dat',sep='\s+',header=None)
data3 = pd.DataFrame(data3)

x3 = data3[0]
y3 = data3[1]
plt.plot(x3, y3,'g.',label='all time steps')

plt.xlabel('t[s]')
plt.ylabel('N[-]')
plt.title('SIMPIC version ***: Number of particles time evolution')
plt.legend()

plt.savefig("plot_nt.png")

plt.clf()

data4 = pd.read_csv('density.dat',sep='\s+',header=None)
data4 = pd.DataFrame(data4)

initial = data4.iloc[:500]

x4i = initial[1]
y4i = initial[2]
plt.plot(x4i, y4i,'b-',label='first time step')

plt.xlabel('x[m]')
plt.ylabel('density[$\mathrm{m}^{-3}$]')
plt.title('SIMPIC version ***: Density of particles calculation')
plt.legend()

plt.savefig("plot_density_first.png")

plt.clf()

last = data4.iloc[-500:]

x4l = last[1]
y4l = last[2]
plt.plot(x4l, y4l,'r-',label='last time step')

plt.xlabel('x[m]')
plt.ylabel('density[$\mathrm{m}^{-3}$]')
plt.title('SIMPIC version ***: Density of particles calculation')
plt.legend()

plt.savefig("plot_density_last.png")