#!/usr/bin/env python3
import os
import sys
import glob

L = 1.0

def readValues(baseName, adjustPosition=False):
    values = {}
    files = glob.glob(f'{baseName}*')
    skipFiles = []
    for file in files:
        padding = len(file[len(baseName):-4])
        if padding == 0:
            skipFiles.append(file)
    [files.remove(_) for _ in skipFiles]
    N = len(files)

    padding = len(file[len(baseName):-4])

    for i in range(N):
        file = f"{baseName}{i:0{padding}.0f}.dat"
        average_boundaries = True
        prevT = '-1'
        with open(file, 'r') as f:
            for j, line in enumerate(f):
                if not line.strip():
                    continue
                split = line.split()
                t = split[0]
                if prevT != t:
                    average_boundaries = True
                    prevT = t
                data = split[1:]
                if adjustPosition:
                    data[0] = str(float(data[0]) + L * i/N)

                if t in values:
                    if i > 0 and average_boundaries:
                        average_boundaries = False
                        previousData = values[t][-1]
                        _N = len(data)
                        averagedData = [0.0 for _ in range(_N)]
                        for k in range(_N):
                            newValue = 0.5 * (float(data[k]) +
                                              float(previousData[k]))
                            averagedData[k] = str(newValue)
                        # Replace last entry with averaged
                        values[t][-1] = (*averagedData,)
                    else:
                        values[t].append((*data,))
                else:
                    values[t] = [(*data,)]
    return values

def writeValues(data, file):
    """Writes values to separate files.
    """
    DEL = ' '

    time = list(data.keys())
    time.sort()

    with open(file, 'w') as f:

        for key in time:
            for el in data[key]:
                line = f'{key}{DEL}'
                line += DEL.join(el) + '\n'
                f.write(line)
            f.write('\n\n')

E = readValues('E', adjustPosition=1)
phi = readValues('phi', adjustPosition=1)
density = readValues('density', adjustPosition=1)
vxx = readValues('vxx')

writeValues(E, 'E.dat')
writeValues(phi, 'phi.dat')
writeValues(density, 'density.dat')
writeValues(vxx, 'vxx.dat')


