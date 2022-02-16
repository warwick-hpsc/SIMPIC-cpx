# -*- coding: utf-8 -*-
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

#Set display Text
rc('text', usetex=True)
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('fontsize')

#Options
variables_list = ['E', 'density', 'phi']


E= np.loadtxt('E.dat')

def get_time_indexes(data_list):
    index_list = [0]
    previous_time = data_list[0,0]
    for index in range(len(data_list[:,0])):
        if data_list[index, 0] == previous_time:
            continue
        else:
            index_list.append(index)
            previous_time = data_list[index,0]
    #print(data_list[index_list, 0])
    #print(data_list[ [index+1 for index in index_list],0])
    #print(data_list[[index-1 for index in index_list],0])
    index_list.append(len(data_list[:,0]) + 1)
    return index_list
    
def plot_Var_vs_Space_forTime(xdata, ydata, time, output_name = '', ytag = '', title = '', suptit = ''):
    plt.figure()
    plt.title(title, fontsize = 15)
    plt.suptitle(suptit)
    plt.xlabel('Space')
    plt.ylabel(ytag)
    plt.plot(xdata, ydata)
    plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
    plt.savefig(output_name + '.png', dpi = 150, bbox_inches = 'tight')
    plt.close()

def plot_all_times(all_data_list, index_list, var,  suptit = '', output_prefix = ''):
    for item in range(len(index_list)-1):
        xdata = all_data_list[index_list[item]:(index_list[item+1]-1),1]
        ydata = all_data_list[index_list[item]:(index_list[item+1]-1),2]
        time  = "{:.5e}".format(all_data_list[index_list[item],0])
        output_name = output_prefix + var + '_' + time
        tit = var + ' plot for time = ' + time + ' s' 
        plot_Var_vs_Space_forTime(xdata, ydata, time , output_name, var, tit, suptit)

def autoplot(variables_list):
    control_get_index = True
    for variable in variables_list:
        print(variable)
        all_values_for_variable = np.loadtxt(variable + '.dat')
        if control_get_index == True:
            indexes = get_time_indexes(all_values_for_variable)
            control_get_indes = False
        plot_all_times(all_values_for_variable, indexes, variable)

autoplot(variables_list)
