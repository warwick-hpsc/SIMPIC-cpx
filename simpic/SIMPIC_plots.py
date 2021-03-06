# -*- coding: utf-8 -*-
"""SIMPIC_plots.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1nSTelqNPGnp-iCuypMBdUnqKzIVQkV9D

- This is a script that plots these four quantities for SIMPIC:

1. Density
2. Phase space (velocity and position)
3. Electric Field
4. Potential

- All data files are taken from the MPI version.
- The following three quantities are varied:
1. Number of particles
2. Number of time steps
3. Voltage on LHS

- Files "*_p##.dat" indicate the number of particles with np=100 being the default in the files with no suffix (e.g. "D.dat").
- Files "*_dt###.dat" indicate the number of time steps with dt=200 being the default in the files with no suffix (e.g. "D.dat").
- Files "*_V##k.dat" indicate the voltage on the LHS electrode with V=25k being the default in the files with no suffix (e.g. "D.dat").
"""

import matplotlib.pyplot as plt
import numpy as np

# Density data
D_dt100 = np.loadtxt('/content/D_dt100.dat')
D_dt300 = np.loadtxt('/content/D_dt300.dat')
D_p50 = np.loadtxt('/content/D_p50.dat')
D = np.loadtxt('/content/D.dat')
D_p200 = np.loadtxt('/content/D_p200.dat')

# Electric field data
E_dt100 = np.loadtxt('/content/E_dt100.dat')
E_dt300 = np.loadtxt('/content/E_dt300.dat')
E_p50 = np.loadtxt('/content/E_p50.dat')
E = np.loadtxt('/content/E.dat')
E_p200 = np.loadtxt('/content/E_p200.dat')

# Potential data
P_dt100 = np.loadtxt('/content/P_dt100.dat')
P_dt300 = np.loadtxt('/content/P_dt300.dat')
P_p50 = np.loadtxt('/content/P_p50.dat')
P = np.loadtxt('/content/P.dat')
P_p200 = np.loadtxt('/content/P_p200.dat')

# Phase space data
PS_dt100 = np.loadtxt('/content/PS_dt100.dat')
PS_dt300 = np.loadtxt('/content/PS_dt300.dat')
PS_p50 = np.loadtxt('/content/PS_p50.dat')
PS = np.loadtxt('/content/PS.dat')
PS_p200 = np.loadtxt('/content/PS_p200.dat')

#LHSV data
E_V20k = np.loadtxt('/content/E_V20k.dat')
E_V30k = np.loadtxt('/content/E_V30k.dat')
P_V20k = np.loadtxt('/content/P_V20k.dat')
P_V30k = np.loadtxt('/content/P_V30k.dat')
PS_V20k = np.loadtxt('/content/PS_V20k.dat')
PS_V30k = np.loadtxt('/content/PS_V30k.dat')
D_V20k = np.loadtxt('/content/D_V20k.dat')
D_V30k = np.loadtxt('/content/D_V30k.dat')

"""-----------------------------------------------------------------
PLOTS OF DENSITY
"""

plt.figure()
plt.suptitle('Plot of Density vs Position')
plt.title('(for varying number of particles)')
plt.xlabel('Position')
plt.ylabel('Particle Density')
plt.ylim(0.9e13,1.1e13)
plt.scatter(D[::100,1],D[::100,2], label='np=100')
plt.scatter(D_p200[::100,1],D_p200[::100,2], label='np=200')
plt.scatter(D_p50[::100,1],D_p50[::100,2], label='np=50')
plt.legend()
plt.show()

plt.figure()
plt.suptitle('Plot of Density vs Position')
plt.title('(for varying number of timesteps)')
plt.xlabel('Position')
plt.ylabel('Particle Density')
plt.ylim(0.9e13,1.1e13)
plt.scatter(D[::100,1],D[::100,2], label='dt=200')
plt.scatter(D_dt300[::100,1],D_dt300[::100,2], label='dt=300')
plt.scatter(D_dt100[::100,1],D_dt100[::100,2], label='dt=100')
plt.legend()
plt.show()

plt.figure()
plt.suptitle('Plot of Density vs Time')
plt.title('(for varying number of particles)')
plt.xlabel('Time')
plt.ylabel('Particle Density')
plt.ylim(0.9e13,1.1e13)
plt.scatter(D[::4*64,0],D[::4*64,2], label='np=100')
plt.scatter(D_p200[::4*64,0],D_p200[::4*64,2], label='np=200')
plt.scatter(D_p50[::4*64,0],D_p50[::4*64,2], label='np=50')
plt.legend()
plt.show()

plt.figure()
plt.suptitle('Plot of Density vs Time')
plt.title('(for varying number of timesteps)')
plt.xlabel('Time')
plt.ylabel('Particle Density')
plt.ylim(0.9e13,1.1e13)
plt.scatter(D[::100,0],D[::100,2], label='dt=200')
plt.scatter(D_dt300[::100,0],D_dt300[::100,2], label='dt=300')
plt.scatter(D_dt100[::100,0],D_dt100[::100,2], label='dt=100')
plt.legend()
plt.show()

"""- Particle density only has very small variations ~0.0025% over position and time. 
- There is no change in density for varying number of particles or number of timesteps (all datapoints overlapped).

-----------------------------------------------------------------
PLOTS OF PHASE SPACE
"""

plt.figure()
plt.suptitle('Plot of Velocity vs Position')
plt.title('(for varying number of particles)')
plt.xlabel('Position')
plt.ylabel('Particle Velocity')
plt.scatter(PS[::4*64,1],PS[::4*64,2], label='np=100')
plt.scatter(PS_p200[::4*64,1],PS_p200[::4*64,2], label='np=200')
plt.scatter(PS_p50[::4*64,1],PS_p50[::4*64,2], label='np=50')
plt.legend()
plt.show()

plt.figure()
plt.suptitle('Plot of Velocity vs Position')
plt.title('(for varying number of timesteps)')
plt.xlabel('Position')
plt.ylabel('Particle Velocity')
plt.scatter(PS[::4*64,1],PS[::4*64,2], label='dt=200')
plt.scatter(PS_dt300[::4*64,1],PS_dt300[::4*64,2], label='dt=300')
plt.scatter(PS_dt100[::4*64,1],PS_dt100[::4*64,2], label='dt=100')
plt.legend()
plt.show()

"""- No plots included for postition vs time as the phase space files only contain data for one timestep ==> time is constant wrt position.
- There is no change in velocity for varying number of particles or number of timesteps (all datapoints overlapped)

-----------------------------------------------------------------
PLOTS OF ELECTRIC FIELD
"""

plt.figure()
plt.suptitle('Plot of Electric Field vs Position')
plt.title('(for varying number of particles)')
plt.xlabel('Position')
plt.ylabel('Electric Field')
plt.scatter(E[::100,1],E[::100,2], label='np=100')
plt.scatter(E_p200[::100,1],E_p200[::100,2], label='np=200')
plt.scatter(E_p50[::100,1],E_p50[::100,2], label='np=50')
plt.legend()
plt.show()

plt.figure()
plt.suptitle('Plot of Electric Field vs Time on LHS')
plt.title('(for varying number of particles)')
plt.xlabel('Time')
plt.ylabel('Electric Field')
plt.scatter(E[::4*64+1,0],E[::4*64+1,2], label='np=100')
plt.scatter(E_p200[::4*64+1,0],E_p200[::4*64+1,2], label='np=200')
plt.scatter(E_p50[::4*64+1,0],E_p50[::4*64+1,2], label='np=50')
plt.legend()
plt.show()

plt.figure()
plt.suptitle('Plot of Electric Field vs Position')
plt.title('(for varying number of timesteps)')
plt.xlabel('Position')
plt.ylabel('Electric Field')
plt.scatter(E[::100,1],E[::100,2], label='dt=200')
plt.scatter(E_dt300[::100,1],E_dt300[::100,2], label='dt=300')
plt.scatter(E_dt100[::100,1],E_dt100[::100,2], label='dt=100')
plt.legend()
plt.show()

plt.figure()
plt.suptitle('Plot of Electric Field vs Time')
plt.title('(for varying number of timesteps)')
plt.xlabel('Time')
plt.ylabel('Electric Field')
plt.scatter(E[::4*64+1,0],E[::4*64+1,2], label='dt=200')
plt.scatter(E_dt300[::4*64+1,0],E_dt300[::4*64+1,2], label='dt=300')
plt.scatter(E_dt100[::4*64+1,0],E_dt100[::4*64+1,2], label='dt=100')
plt.legend()
plt.show()

"""- Electric field in the plots vs time is measured on the LHS plate.
- There is no change in electric field for varying number of particles or number of timesteps (all datapoints overlapped)

-----------------------------------------------------------------
PLOTS OF POTENTIAL
"""

plt.figure()
plt.suptitle('Plot of Potential vs Position')
plt.title('(for varying number of particles)')
plt.xlabel('Position')
plt.ylabel('Potential')
plt.scatter(P[::100,1],P[::100,2], label='np=100')
plt.scatter(P_p200[::100,1],P_p200[::100,2], label='np=200')
plt.scatter(P_p50[::100,1],P_p50[::100,2], label='np=50')
plt.legend()
plt.show()

plt.figure()
plt.suptitle('Plot of Potential vs Time')
plt.title('(for varying number of particles)')
plt.xlabel('Time')
plt.ylabel('Potential')
plt.scatter(P[::4*64+1,0],P[::4*64+1,2], label='np=100')
plt.scatter(P_p200[::4*64+1,0],P_p200[::4*64+1,2], label='np=200')
plt.scatter(P_p50[::4*64+1,0],P_p50[::4*64+1,2], label='np=50')
plt.legend()
plt.show()

plt.figure()
plt.suptitle('Plot of Potential vs Position')
plt.title('(for varying number of timesteps)')
plt.xlabel('Position')
plt.ylabel('Potential')
plt.scatter(P[::100,1],P[::100,2], label='dt=200')
plt.scatter(P_dt300[::100,1],P_dt300[::100,2], label='dt=300')
plt.scatter(P_dt100[::100,1],P_dt100[::100,2], label='dt=100')
plt.legend()
plt.show()

plt.figure()
plt.suptitle('Plot of Potential vs Time')
plt.title('(for varying number of timesteps)')
plt.xlabel('Time')
plt.ylabel('Potential')
plt.scatter(P[::4*64+1,0],P[::4*64+1,2], label='dt=200')
plt.scatter(P_dt300[::4*64+1,0],P_dt300[::4*64+1,2], label='dt=300')
plt.scatter(P_dt100[::4*64+1,0],P_dt100[::4*64+1,2], label='dt=100')
plt.legend()
plt.show()

"""- Potential in the plots vs time is measured on the LHS plate.
- There is no change in potential for varying number of particles or number of timesteps (all datapoints overlapped)

-----------------------------------------------------------------
VARYING LHS VOLTAGE
"""

plt.figure()
plt.suptitle('Plot of Density vs Position')
plt.title('(for varying LHS voltage)')
plt.xlabel('Position')
plt.ylabel('Density')
#plt.ylim(0.1e13,1.2e13)
plt.scatter(D_V30k[::100,1],D_V30k[::100,2], label='V=20k')
plt.scatter(D_V20k[::100,1],D_V20k[::100,2], label='V=30k')
plt.scatter(D[::100,1],D[::100,2], label='V=25k')
plt.legend()
plt.show()

plt.figure()
plt.suptitle('Plot of Density vs Time')
plt.title('(for varying LHS voltage)')
plt.xlabel('Time')
plt.ylabel('Density')
plt.scatter(D[::4*64,0],D[::4*64,2], label='V=25k')
plt.scatter(D_V20k[::4*64,0],D_V20k[::4*64,2], label='V=20k')
plt.scatter(D_V30k[::4*64,0],D_V30k[::4*64,2], label='V=30k')
plt.legend()
plt.show()

plt.figure()
plt.suptitle('Plot of Velocity vs Position')
plt.title('(for varying LHS voltage)')
plt.xlabel('Position')
plt.ylabel('Velocity')
plt.plot(PS_V30k[::100,1],PS_V30k[::100,2], label='V=20k')
plt.plot(PS_V20k[::100,1],PS_V20k[::100,2], label='V=30k')
plt.plot(PS[::100,1],PS[::100,2], label='V=25k')
plt.legend()
plt.show()

plt.figure()
plt.suptitle('Plot of Velocity vs Time')
plt.title('(for varying LHS voltage)')
plt.ylabel('Time')
plt.xlabel('Velocity')
plt.scatter(PS[::4*64+1,2],PS[::4*64+1,0], label='V=25k')
plt.scatter(PS_V20k[::4*64+1,2],PS_V20k[::4*64+1,0], label='V=20k')
plt.scatter(PS_V30k[::4*64+1,2],PS_V30k[::4*64+1,0], label='V=30k')
plt.legend()
plt.show()

plt.figure()
plt.suptitle('Plot of Electric Field vs Position')
plt.title('(for varying LHS voltage)')
plt.xlabel('Position')
plt.ylabel('Electric Field')
plt.plot(E_V30k[::100,1],E_V30k[::100,2], label='V=20k')
plt.plot(E_V20k[::100,1],E_V20k[::100,2], label='V=30k', color='green')
plt.plot(E[::100,1],E[::100,2], label='V=25k', color='orange')
plt.legend()
plt.show()

plt.figure()
plt.suptitle('Plot of Electric Field vs Time')
plt.title('(for varying LHS voltage)')
plt.xlabel('Time')
plt.ylabel('Electric Field')
plt.scatter(E[::4*64+1,0],E[::4*64+1,2], label='V=25k')
plt.scatter(E_V20k[::4*64+1,0],E_V20k[::4*64+1,2], label='V=20k')
plt.scatter(E_V30k[::4*64+1,0],E_V30k[::4*64+1,2], label='V=30k')
plt.legend()
plt.show()

plt.figure()
plt.suptitle('Plot of Potential vs Position')
plt.title('(for varying LHS voltage)')
plt.xlabel('Position')
plt.ylabel('Potential')
plt.scatter(P[::100,1],P[::100,2], label='V=25k')
plt.scatter(P_V30k[::100,1],P_V30k[::100,2], label='V=20k')
plt.scatter(P_V20k[::100,1],P_V20k[::100,2], label='V=30k')
plt.legend()
plt.show()

plt.figure()
plt.suptitle('Plot of Potential vs Time')
plt.title('(for varying LHS voltage)')
plt.xlabel('Time')
plt.ylabel('Potential')
plt.scatter(P[::4*64+1,0],P[::4*64+1,2], label='V=25k')
plt.scatter(P_V30k[::4*64+1,0],P_V30k[::4*64+1,2], label='V=20k')
plt.scatter(P_V20k[::4*64+1,0],P_V20k[::4*64+1,2], label='V=30k')
plt.legend()
plt.show()
