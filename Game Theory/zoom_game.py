#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 00:29:30 2020

@author: Clement

Modeling Zoom Game
"""
import pandas
import numpy
import matplotlib.pyplot as plt
import os
from cycler import cycler

def choose_plan (current, list_t, list_student, list_pi, v_i, list_ei, name_set, name_trial): #Determine if set and trial exist already 
    list_subdirr1 = []
    for root, dirs, files in os.walk(current+'/Data set'):
        list_subdirr1.append(dirs)
    
    if name_set in list_subdirr1[0]:
        list_subdirr2 = []
        for root, dirs, files in os.walk(current+f'/Data set/{name_set}'):
            list_subdirr2.append(dirs)
    
        df_student = importing_student_set (current, name_set)
        
        if name_trial in list_subdirr2[0]:
            print(f'Replacing {name_trial} in existing {name_set}')
            
        else:
            print(f'Creation of {name_trial} in existing {name_set}')
            
    else:
        df_student = creating_student_set(list_t, list_student, list_pi, v_i, list_ei, name_set)
        print(f'Creation of {name_trial} in {name_set}')
    print('\n\n')           
    return df_student
         
def creation_folder (paths): #General function to create a folder
    current = os.path.dirname(os.path.realpath(__file__))
    list_directory = [os.path.normcase(current + directory) for directory in paths]
    
    for directory in list_directory:
        if os.path.exists(directory) == False:
            print('Directory created')
            os.makedirs(directory)
    list_return = [current + x + '/' for x in paths]      
    return list_return

def random_draw (mu, sigma, size): #Draw values and then filter out the one not between 0 and 1
    normal_draw = numpy.random.normal(loc=mu, scale=sigma, size=3*size)
    k=0
    
    normal_draw_mod = []
    
    for val in normal_draw:
        if val < 1 and val > 0 and k<size:
            normal_draw_mod.append(val)
            k+=1
            
    return normal_draw_mod

def normal_distrib (lim_x, mu, sigma): #Compute the theo value of normal distrib
    list_x = numpy.linspace(lim_x[0], lim_x[1], 100)
    list_y = [1/(sigma * numpy.sqrt(2 * numpy.pi)) * 
              numpy.exp( - (x - mu)**2 / (2 * sigma**2) ) for x in list_x]
    
    return list_x, list_y

def creating_student_set (list_t, list_student, list_pi, v_i, list_ei, name_set): #Create a new set of student
    student_index= pandas.Index(list_student, name='student')
    df_student = pandas.DataFrame(index = student_index, columns = ['p_i', 'v_i', 'e_i'])
    
    df_student.loc[:, 'p_i'] = list_pi
    df_student.loc[:, 'e_i'] = list_ei
    df_student.loc[:, 'v_i'] = v_i
    df_student.loc[:, 'I1'] = 2*df_student.loc[:, 'e_i']*(1-df_student.loc[:, 'v_i'])
    df_student.loc[:, 'I2'] = 2*df_student.loc[:, 'e_i']*(1+df_student.loc[:, 'v_i'])
    
    file_path = creation_folder([f'/Data set/{name_set}'])
    df_student.to_csv(os.path.normcase(file_path [0] + 'student_set.csv'), index=True)
    
    return df_student

def importing_student_set (current, name_set): #Import an existing set of studeny
    df_student = pandas.read_csv(os.path.normcase(f'{current}/Data set/{name_set}/student_set.csv'))
    df_student = df_student.set_index('student')
    
    return df_student
    
def creating_experience_df (list_student, list_t, df_full_reslts): #Create blank df used to save the game output
    df_camera = pandas.DataFrame(index = list_student, columns = list_t)
    df_payoff = pandas.DataFrame(index = list_student, columns = list_t)
    
    df_result_index = pandas.Index(list_t, name='time')
    df_result = pandas.DataFrame(index = df_result_index, columns = df_full_reslts.columns)
    
    return df_camera, df_payoff, df_result
    
def one_period (df_student, df_camera, df_payoff, df_result, coeff): #Detrmining payoff an camera status for one game
    list_t = df_camera.columns
    number_student = df_student.index.size
    
    for student in df_student.index:
        if coeff < df_student.loc[student, 'I1']:
            df_camera.loc[student, 0] = 0
            df_payoff.loc[student, 0] = df_student.loc[student, 'e_i']
            
        else:
            df_camera.loc[student, 0] = 1
            df_payoff.loc[student, 0] = coeff - df_student.loc[student, 'e_i']
    
    df_result.loc[0, 'N_C'] = df_camera.loc[:,0].sum()
    df_result.loc[0, 'PO(Te)'] = df_result.loc[0, 'N_C']/number_student
    df_result.loc[0, 'PO(St)'] = df_payoff.loc[:,0].sum()/number_student
    
    for prev_t, t in zip(list_t[:-1], list_t[1:]):
        df_camera_on = df_camera[(df_camera.loc[:,prev_t]==1)]
        for student in df_camera_on.index:
            if coeff <= df_student.loc[student, 'I2']:
                if df_result.loc[prev_t, 'N_C']/number_student < df_student.loc[student, 'p_i']:
                    df_camera.loc[student, t] = 0
                    df_payoff.loc[student, t] = df_student.loc[student, 'e_i']
        
        df_camera = df_camera.fillna(axis=1, method='ffill', limit=1)
        df_payoff = df_payoff.fillna(axis=1, method='ffill', limit=1)
        
        df_result.loc[t, 'N_C'] = df_camera.loc[:,t].sum()
        df_result.loc[t, 'PO(Te)'] = df_result.loc[t, 'N_C']/number_student
        df_result.loc[t, 'PO(St)'] = df_payoff.loc[:,t].sum()/number_student
         
    return df_camera, df_payoff, df_result

def several_periods (list_student, list_t, list_coeff): #Game for several coeff x
    full_result_index = pandas.MultiIndex.from_product([list_coeff, list_t], names=['coeff', 'time'])
    df_full_results = pandas.DataFrame(index=full_result_index, columns = ['N_C', 'PO(Te)', 'PO(St)'])
    
    for coeff in list_coeff:
        df_camera, df_payoff, df_result = creating_experience_df (list_student, list_t, df_full_results)
        df_camera, df_payoff, df_result = one_period (df_student, df_camera, df_payoff, df_result, coeff)
        df_full_results.loc[coeff] = df_result.values
        
    print(df_full_results)   
    return df_full_results

def plot (df_student, normal_distrib_ei, normal_distrib_pi, df_full_results, style_cycle, name_set, name_trial): #plotting results
    fig, axs = plt.subplots(2,2, figsize=(11.7, 8.3), num=f'Result Zoom Game {name_set} {name_trial}')
    fig.suptitle('Result of Zoom Game', fontsize=14)
    list_coeff = df_full_results.index.get_level_values('coeff').unique()
    
    axs[0][0].set_prop_cycle(style_cycle)
    for a_coeff in list_coeff :
        list_t = df_full_results.index.get_level_values('time').unique()
        axs[0][0].plot(list_t, df_full_results.loc[a_coeff, 'N_C'], label=f'$x=${round(a_coeff,2)}')
    
    axs[0][1].set_prop_cycle(style_cycle)
    axs[0][1].plot(list_coeff, df_full_results.groupby(level='coeff').tail(1).loc[:,'PO(Te)'], label='Teacher')
    axs[0][1].plot(list_coeff, df_full_results.groupby(level='coeff').tail(1).loc[:,'PO(St)'], label='Student')
    
    axs[1][0].hist(df_student.loc[:,'e_i'], bins=30, density=True, facecolor=color_c.by_key()['color'][0], edgecolor='black', linewidth=1, zorder=4)
    axs[1][0].plot(normal_distrib_ei[0], normal_distrib_ei[1], zorder=5, color=color_c.by_key()['color'][1], linewidth=1)
    
    axs[1][1].hist(df_student.loc[:,'p_i'], bins=30, density=True, facecolor=color_c.by_key()['color'][0], edgecolor='black', linewidth=1, zorder=4)
    axs[1][1].plot(normal_distrib_pi[0], normal_distrib_pi[1], zorder=5, color=color_c.by_key()['color'][1], linewidth=1)
    
    for axes, graph in zip(numpy.ravel(axs), df_plot.index):
        axes.set_xlabel(df_plot.loc[graph, 'x_title'])
        axes.set_ylabel(df_plot.loc[graph, 'y_title'])
        axes.set_title (df_plot.loc[graph, 'title'], fontsize=12)
        axes.grid(zorder=0)
        
        if df_plot.loc[graph, 'legend']:
            axes.legend(title=df_plot.loc[graph, 'legend_title'])
            
    fig.subplots_adjust(hspace = 0.3)
    
    file_paths = creation_folder ([f'/Data set/{name_set}/{name_trial}']) 
    fig.savefig(os.path.normcase(file_paths[0] + 'Graph_output.pdf'), format='pdf')
    fig.savefig(os.path.normcase(file_paths[0] + 'Graph_output.png'), format='png', dpi=300)
    fig.savefig(os.path.normcase(file_paths[0] + 'Graph_output.svg'), format='svg')

def export_save_para (df_full_results, list_student, para_ei, para_pi, v_i, list_t, list_x, name_set, name_trial):
    file_paths = creation_folder ([f'/Data set/{name_set}/{name_trial}']) 
    
    df_full_results.to_csv(os.path.normcase(file_paths[0] + 'results.csv'), index=True)
    
    list_lines = []
    list_lines.append(f'Parameters for {name_set}, {name_trial}\n')
    
    list_lines.append(f'Student set parameters ({name_set})')
    list_lines.append(f'Number of student: {len(list_student)}')
    list_lines.append(f'Effort values: {para_ei}')
    list_lines.append(f'Uncertainty on effort value: v_i={v_i}')
    list_lines.append(f'Minium portion of class with camera on: {para_pi}\n')
    
    list_lines.append(f'Game parameters ({name_trial})')
    list_lines.append(f'Number of period: {len(list_t)}')
    list_lines.append(f'List of coefficient: {list_x}')
    
    export_file = open(os.path.normcase(file_paths[0] + f'Parameters ({name_set}, {name_trial}).txt'), 'w')
    for line in list_lines:
        export_file.write(line+'\n')
    export_file.close()
    
    
# --- General parameters ---
current = os.path.dirname(os.path.realpath(__file__))
colors= ['#386cb0', '#D45050', '#7fc97f', '#9A46C4', '#F08328', '#a6cee3', 'k']
color_c = cycler(color=['#386cb0', '#D45050', '#7fc97f', '#9A46C4', '#F08328', '#a6cee3', 'k'])
marker_c = cycler(marker=['2', 4, 'o', 'x', 's', None])
line_c = cycler(linestyle=['solid', 'dashed'])
linewidth_c = cycler(linewidth=[1])
cycle_color = line_c[:1] * linewidth_c * marker_c[-1:] * color_c
style_cycle = cycle_color

# --- Setting up the parameters for the student set ---
name_set = 'Set 3' #Name of set for saving purposes

size = 100 #Students
list_student = numpy.arange(1, size+1,1)

mu1, sigma1 = 0.1, 0.3 #Normal distrib for the e_i values
para_ei = (f'Normal distrib with μ={mu1} and σ={sigma1}')
list_ei = random_draw (mu1, sigma1, size)

#list_pi = numpy.linspace(0.1, 1.0, 6)
mu2, sigma2 = 0.5, 0.3 #Normal distrib for the p_i values
para_pi = (f'Normal distrib with μ={mu2} and σ={sigma2}')
list_pi = random_draw (mu2, sigma2, size)

v_i = 0.2 #Constant value for the uncertainty

# --- Setting up the game parameter ---
name_trial = 'Trial 1' #Name for saving
list_t = numpy.arange(0, 7, 1) #Number of period
list_x = numpy.linspace(0.1, 0.3, 3) #Coeff

#Title for the plot
df_plot = pandas.DataFrame(index= ['Plot1', 'Plot2', 'Plot3', 'Plot4'],
                           columns = ['title', 'x_title', 'y_title', 'legend', 'legend_title'],
                           data =[['Number of students with camera on\nfor several coefficient $x$', 'Subperiod $t_k$', 'Number of student with camera on $N_{k, x}$', True, 'Coefficient'],
                                    ["Teacher's and averege student final payoff\nfor several coefficient $x$", 'Coefficient $x$', r"Payoff $PO$", True, 'Payoff'],
                                    ['Histogram of the effort values $e_i$', 'Value of $e_i$', 'Density', False,''],
                                    ['Histogram of the $p_i$ values', 'Value of $p_i$', 'Density', False,'']])


# --- Main program ---
df_student = choose_plan (current, list_t, list_student, list_pi, v_i, list_ei, name_set, name_trial) #Create a new set or reuse the same set
df_full_results = several_periods (list_student, list_t, list_x) #Compute the outcome of the game

normal_distrib_ei = normal_distrib ([df_student.loc[:,'e_i'].min(), df_student.loc[:,'e_i'].max()], mu1, sigma1) #Calculating normal dstrib for e_i
normal_distrib_pi = normal_distrib ([df_student.loc[:,'p_i'].min(), df_student.loc[:,'p_i'].max()], mu2, sigma2) #Calculating normal dstrib for p_i

plot (df_student, normal_distrib_ei, normal_distrib_pi, df_full_results, style_cycle, name_set, name_trial) #Plot results
export_save_para (df_full_results, list_student, para_ei, para_pi, v_i, list_t, list_x, name_set, name_trial) #Save parameters used for game and save results
