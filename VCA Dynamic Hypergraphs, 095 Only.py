#!/usr/bin/env python
# coding: utf-8

# In[16]:


import sys
import os
import pandas
from pandas import ExcelWriter
import matplotlib
from matplotlib import pyplot as plt
import statistics
import math

# Authors: Ashti M. Shah, University of Pittsburgh School of Medicine 
# Mentors: Dr. Yoram Vodovotz and Dr. Ruben Zamora
# September, 2022

# This specific code is written to create the spatio-temporal mapping of inflammatory mediators across 3 different
# tissues over three day time intervals. However, the code itself is easily modifiable to incorporate multiple
# tissues or other types of data that are apt for for spatio-temporal modeling. 


# Begin by importing data: the input spreadsheet should be organized such that the first column specifies 
# the tissue (ex. muscle, plasma, skin). The second column should specify the time point of the sample (ex. d0
# d3, d5). Each remaining column is a specific inflammatory mediator and the row contains the amount of mediator
# (ex. pg/mL, or some other unit) within the specified tissue at the indicated time point.


"Import data"
raw_data = pandas.DataFrame()
raw_data = pandas.read_csv('A3_A4_Input.csv') #CHANGE INPUT FILE NAME
cytokines = list(raw_data.columns.values)[2:] #list of all cytokines

"Sort data by tissue type"
muscle_raw_data = raw_data.loc[raw_data['Tissue'] == 'Muscle']
skin_raw_data = raw_data.loc[raw_data['Tissue'] == 'Skin']
plasma_raw_data = raw_data.loc[raw_data['Tissue'] == 'Plasma']

# For our computational analysis, we use the median value of a cytokine at a given time point across all samples

"Function to get median value for each cytokine and time point by tissue"
def get_median_for_tissue(tissue_data, cytokines):
    "tissue_data: a pandas data frame of of all the data from a specific tissue"
    "cytokines: A list of all the cytokines in the data table"
    results = pandas.DataFrame()
    all_times_str = ['d0', 'd3','d5','d7','d9','d11','d20','d23','d25','d27','d31']
    results['Day']= all_times_str
    for c in cytokines:
        "cur_df: data frame where the first column is all the time points in a tissue and the second"
        "column is the cytokine quantification of cytokine, c"
        cur_df = pandas.DataFrame()
        cur_df['Time'] = tissue_data['Day']
        cur_df[c] = tissue_data [c]
        list_cur_median_cyt = []
        for day in all_times_str:
            cur_day = cur_df.loc[cur_df['Time'] == day] #data frame of all samples at unique time point
                                                        #and cytokine quantification
            list_cur_median_cyt.append(cur_day[c].median()) #calculate the median for a given time point
        results[c]=list_cur_median_cyt
    return results

#Get the median cytokine value across all samples for each cytokine at each time point

MUSCLE = get_median_for_tissue(muscle_raw_data, cytokines)
SKIN = get_median_for_tissue(skin_raw_data, cytokines)
PLASMA = get_median_for_tissue(plasma_raw_data, cytokines)

#Get a focused version of the data frame MUSCLE, SKIN, or PLASMA, that includes only three consecutive
#time points

"Function to get a Pandas DF of relevant data for correlation matrix:"
"Rows are 3 consecutive time points, columns are cytokine values (pg/mL)"
def get_rel_data (tissue_data, time):
    "tissue_data: a pandas data frame of of all the data from a specific tissue"
    "time: an array of the strings of time names, ex: ['d0', 'd3', 'd5']"
    "rel_data: a df containing data for a given tissue for just three consecutive time points"
    rel_data = tissue_data.loc[tissue_data['Day'] == time[0]]
    rel_data = rel_data.append(tissue_data.loc[tissue_data['Day'] == time[1]])
    rel_data = rel_data.append(tissue_data.loc[tissue_data['Day'] == time[2]])
    return rel_data

#Get a table of the pearson's correlation coefficient across three time points for each cytokine within a tissue

"Function to calculate significant cytokines within a given tissue"
def get_significant_cytokines (dynamic_data, time, cytokines, tissue_name):
    "dynamic_data: pandas data frame of cytokines values in one tissue from three consecutive time points"
    "time: An array of three integers specifying the time points, ex: [0, 3, 5]"
    "cytokines: an array of cytokine names"
    "tissue name: a string specifying the name of the tissue, ex: 'Muscle' or 'Skin'"
    "significant data: a table including the cytokines where |r| > 0.7"
    significant_data = pandas.DataFrame()
    list_corr_coeff = []
    for c in cytokines:
        y_val = dynamic_data[c] #a single cytokine across three time points
        df = pandas.DataFrame(y_val)
        df['time'] = time
        df = df[['time', c]] #df is a data frame with three time points in column 1 and the median cytokine value
                             #for a given cytokine in at that time in column 2
        corr_matrix = df.corr(method = 'pearson')
        r = (corr_matrix.iloc[0,1])
        "If r > 0.7 and < 0.95, then r = 0.7; if r > 0.95, then r = 0.95"
        if r>0.7 and r <0.95:
            r = 0.7
        elif r >= 0.95:
            r = 0.95
        elif r < -0.7 and r > -0.95:
            r = -0.7
        elif r <= -0.95:
            r = -0.95
        list_corr_coeff.append(r)
    significant_data['Cytokines in %s'%tissue_name] = cytokines
    significant_data['Pearsons R - %s'%tissue_name] = list_corr_coeff    
    significant_data = significant_data[abs(significant_data['Pearsons R - %s'%tissue_name]) >= 0.7]
    return significant_data

"Function to get all Pearsons correlations for a given time interval in all three tissues"
def get_dynamic_interval_data (Muscle, Skin, Plasma, time_str, time_int, cytokines):
    "Muscle: Muscle_Mouse1; All data from muscle for a specific mouse"
    "Skin: Skin_Mouse1; All data from skin for a specific mouse"
    "Plasma: Plasma_Mouse1; All data from plasma for a specific mouse"
    "time_str: list of strings that represent time, ex: ['d0','d3','d5']"
    "time_int: list of times as integers, ex: [0, 3, 5]"
    "cytokines: list of cytokines"
    Muscle_data = get_rel_data(Muscle, time_str)
    rel_muscle = get_significant_cytokines(Muscle_data, time_int, cytokines,'Muscle')
    Skin_data = get_rel_data(Skin, time_str)
    rel_skin = get_significant_cytokines (Skin_data, time_int, cytokines, 'Skin')
    Plasma_data = get_rel_data(Plasma, time_str)
    rel_plasma = get_significant_cytokines(Plasma_data, time_int, cytokines, 'Plasma')
    results = pandas.concat([rel_muscle, rel_skin, rel_plasma], axis=1)
    return results

# Function to sort cytokines into groups of tissue. For ex: if IL-17A appears in the dictionary key 'muscle'
# and the dictionary key 'plasma', then we move IL-17A to the dictionary key 'muscle and plasma'
# This function is written to be used in hypergraph_grouped_edges

def get_hypergraph_combined_edges(dict_nodes_edges):
    "dict_nodes_edges: a dictionary where the keys are nodes and the definitions are a list of cytokines"
    "found within that node. Only the keys 'muscle', 'skin', and 'plasma' have definitions at the beginning"
    "of this function. By the end, edges that surround multiple nodes are sorted into the appropriate key"
    "this new dictionary is returned"
    for m in dict_nodes_edges['muscle']:
        if m in dict_nodes_edges['skin'] and m in dict_nodes_edges['plasma']:
            dict_nodes_edges['muscle, skin, and plasma'].append(m)
        elif m in dict_nodes_edges['skin']:
            dict_nodes_edges['muscle and skin'].append(m)
        elif m in dict_nodes_edges['plasma']:
            dict_nodes_edges['muscle and plasma'].append(m)     
    for s in dict_nodes_edges['skin']:
        if s in dict_nodes_edges['plasma'] and s not in dict_nodes_edges['muscle']:
            dict_nodes_edges['skin and plasma'].append(s)
    for msp in dict_nodes_edges['muscle, skin, and plasma']:
        dict_nodes_edges['muscle'].remove(msp)
        dict_nodes_edges['skin'].remove(msp)
        dict_nodes_edges['plasma'].remove(msp)
    for mp in dict_nodes_edges['muscle and plasma']:
        dict_nodes_edges['muscle'].remove(mp)
        dict_nodes_edges['plasma'].remove(mp)
    for ms in dict_nodes_edges['muscle and skin']:
        dict_nodes_edges['muscle'].remove(ms)
        dict_nodes_edges['skin'].remove (ms)
    for sp in dict_nodes_edges['skin and plasma']:
        dict_nodes_edges['skin'].remove(sp)
        dict_nodes_edges['plasma'].remove(sp)
    return dict_nodes_edges

# Function to group all of the edges into subsets of nodes (ex. muscle and skin, muscle, and plasma,
# muscle alone)

def hypergraphs_grouped_edges(cur_graph):
    "cur_graph: a data frame containing 6 columns organized as --> cytokines in muscle, pearson's r muscle, etc."
    "results: a data frame of cytokines organized by groups of tissues that they appear in as well as "
    "pearson's r"
    "Create dictionaries to sort cytokines that appear in multiple tissues"
    edge_neg_095 ={'muscle':[], 'skin':[], 'plasma':[], 'muscle and skin':[], 'muscle and plasma':[], 'skin and plasma':[], 'muscle, skin, and plasma':[]}
    edge_neg_07 = {'muscle':[], 'skin':[], 'plasma':[], 'muscle and skin':[], 'muscle and plasma':[], 'skin and plasma':[],'muscle, skin, and plasma':[]}
    edge_pos_095 = {'muscle':[], 'skin':[], 'plasma':[], 'muscle and skin':[], 'muscle and plasma':[], 'skin and plasma':[],'muscle, skin, and plasma':[]}
    edge_pos_07 = {'muscle':[], 'skin':[], 'plasma':[], 'muscle and skin':[], 'muscle and plasma':[], 'skin and plasma':[], 'muscle, skin, and plasma':[]}
    
    # iterate through each row of cur_graph     
    for index, row in cur_graph.iterrows():
        muscle_sig = row['Pearsons R - Muscle']
        skin_sig = row ['Pearsons R - Skin']
        plasma_sig = row['Pearsons R - Plasma']
        # Dependeing on the value of r, append the cytokine in muscle, skin, or plasma into its
        # respective dictionary of edges
        if math.isnan(muscle_sig) == False:
            if muscle_sig == 0.7:
                edge_pos_07['muscle'].append(row['Cytokines in Muscle'])
            elif muscle_sig == 0.95:
                edge_pos_095['muscle'].append(row['Cytokines in Muscle'])
            elif muscle_sig == -0.7:
                edge_neg_07['muscle'].append(row['Cytokines in Muscle'])
            elif muscle_sig == -0.95:
                edge_neg_095['muscle'].append(row['Cytokines in Muscle'])
        if math.isnan(skin_sig) == False:
            if skin_sig == 0.7:
                edge_pos_07['skin'].append(row['Cytokines in Skin'])
            elif skin_sig == 0.95:
                edge_pos_095['skin'].append(row['Cytokines in Skin'])
            elif skin_sig == -0.7:
                edge_neg_07['skin'].append(row['Cytokines in Skin'])
            elif skin_sig == -0.95:
                edge_neg_095['skin'].append(row['Cytokines in Skin'])
        if math.isnan(plasma_sig) == False:
            if plasma_sig == 0.7:
                edge_pos_07['plasma'].append(row['Cytokines in Plasma'])
            elif plasma_sig == 0.95:
                edge_pos_095['plasma'].append(row['Cytokines in Plasma'])
            elif plasma_sig == -0.7:
                edge_neg_07['plasma'].append(row['Cytokines in Plasma'])
            elif plasma_sig == -0.95:
                edge_neg_095['plasma'].append(row['Cytokines in Plasma'])

    edge_neg_095_sorted = pandas.DataFrame.from_dict(get_hypergraph_combined_edges(edge_neg_095), orient = 'index')
    edge_neg_07_sorted = pandas.DataFrame.from_dict(get_hypergraph_combined_edges(edge_neg_07), orient = 'index')
    edge_pos_095_sorted = pandas.DataFrame.from_dict(get_hypergraph_combined_edges(edge_pos_095), orient = 'index')
    edge_pos_07_sorted = pandas.DataFrame.from_dict(get_hypergraph_combined_edges(edge_pos_07), orient = 'index')
    keys = ['Edge = 0.95', 'Edge = 0.7', 'Edge = -0.95', 'Edge = -0.7']
    results = pandas.concat([edge_pos_095_sorted, edge_pos_07_sorted, edge_neg_095_sorted, edge_neg_07_sorted], axis = 1, keys = keys)   
    return results

#Function to plot a singular dynamic hypergraph
def plot_hypergraph(cur_dict, color, figure, grid_spec, panel, title, thickness):
    "cur_dict: a dictionary where the keys are nodes or groups of nodes and the definitions are edges"
    "ex: edge_neg_095_sorted"
    "color: the color of the graph"
    "figure: the actual figure panel"
    "grid_spec: to space the graphs"
    "panel: A, B, C, or D, to specify the location of the subplot on the overall figure"
    "title: title of the subgraph"
    "thickness: an integer, indicating the thickness of the lines of the graph"
    
    if panel == 'A':
        ax = figure.add_subplot(grid_spec[0,0])
    elif panel == 'B':
        ax = figure.add_subplot(grid_spec[0,1])
    elif panel == 'C':
        ax = figure.add_subplot(grid_spec[1,0])
    else:
        ax = figure.add_subplot(grid_spec[1,1])
        
    plt.title(title, size=36)
    x_coords_edge = [0.75,6]
    # coordinates for the location of groups of edges depending on the edge group (ex: muscle and plasma located
    # at (10, 8.34))
    muscle_only_y = [10, 10]
    skin_only_y =[5, 5]
    plasma_only_y = [0.1, 0.1]
    muscle_skin_y1 = [10, 6.67]
    muscle_skin_y2 = [5, 6.67]
    muscle_plasma_y1 = [10, 8.34]
    muscle_plasma_y2 = [0.1, 8.34]
    plasma_skin_y1 = [0.1, 3.34]
    plasma_skin_y2 = [5, 3.34]
    msp_y1 = [10, 1.67]
    msp_y2 = [5, 1.67]
    msp_y3 = [0.1, 1.67]
    
    muscle_label = 'MUSCLE'
    skin_label = 'SKIN'
    plasma_label = 'PLASMA'
    plt.text(-2,10.2, muscle_label, fontsize=36, verticalalignment='top')
    plt.text(-2,5.2, skin_label, fontsize=36, verticalalignment='top')
    plt.text(-2,0.2, plasma_label, fontsize=36,verticalalignment='top')
    if (cur_dict['muscle']) != []:
        plt.plot(x_coords_edge, muscle_only_y, color, linewidth = thickness)
        if len(cur_dict['muscle']) > 5:
            plt.text(6.5, 10.3, "           ".join(cur_dict['muscle'][0:6]),fontsize=24, verticalalignment='top')
            plt.text(6.5, 10, "           ".join(cur_dict['muscle'][6:]),fontsize=24, verticalalignment='top')
        else:
            plt.text(6.5, 10.2, "           ".join(cur_dict['muscle']),fontsize=24, verticalalignment='top')
    if (cur_dict['skin']) != []:
        plt.plot(x_coords_edge, skin_only_y, color, linewidth = thickness)
        if len(cur_dict['skin']) > 5:
            plt.text(6.5, 5.3, "           ".join(cur_dict['skin'][0:6]),fontsize=24, verticalalignment='top')
            plt.text(6.5, 4.6, "           ".join(cur_dict['skin'][6:]),fontsize=24, verticalalignment='top')
        else:
            plt.text(6.5, 5.2, "           ".join(cur_dict['skin']),fontsize=24, verticalalignment='top')
    if (cur_dict['plasma']) != []:
        plt.plot(x_coords_edge, plasma_only_y, color, linewidth = thickness)
        if len(cur_dict['plasma']) > 5:
            plt.text(6.5, 0.7, "           ".join(cur_dict['plasma'][0:6]),fontsize=24, verticalalignment='top')
            plt.text(6.5, 0.1, "           ".join(cur_dict['plasma'][6:]),fontsize=24, verticalalignment='top')
        else:
            plt.text(6.5, 0.2, "           ".join(cur_dict['plasma']),fontsize= 24, verticalalignment='top')
    if (cur_dict['muscle and skin']) != []:
        plt.plot(x_coords_edge, muscle_skin_y1, color, linewidth = thickness)
        plt.plot(x_coords_edge, muscle_skin_y2, color, linewidth = thickness)
        if len(cur_dict['muscle and skin']) > 5:
            plt.text(6.5, 6.97, "           ".join(cur_dict['muscle and skin'][0:6]),fontsize=24, verticalalignment='top')
            plt.text(6.5, 6.67, "           ".join(cur_dict['muslce and skin'][6:]),fontsize=24, verticalalignment='top')
        else:
            plt.text(6.5, 6.87, "           ".join(cur_dict['muscle and skin']),fontsize=24, verticalalignment='top')
    if (cur_dict['muscle and plasma']) != []:
        plt.plot(x_coords_edge, muscle_plasma_y1, color, linewidth = thickness)
        plt.plot(x_coords_edge, muscle_plasma_y2, color, linewidth = thickness)
        if len(cur_dict['muscle and plasma']) > 5:
            plt.text(6.5, 8.64, "           ".join(cur_dict['muscle and plasma'][0:6]),fontsize=24, verticalalignment='top')
            plt.text(6.5, 8.34, "           ".join(cur_dict['muscle and plasma'][6:]),fontsize=24, verticalalignment='top')
        else: 
            plt.text(6.5, 8.54, "           ".join(cur_dict['muscle and plasma']),fontsize=24, verticalalignment='top')
    if (cur_dict['skin and plasma']) != []:
        plt.plot(x_coords_edge, plasma_skin_y1, color, linewidth = thickness)
        plt.plot(x_coords_edge, plasma_skin_y2, color, linewidth = thickness)
        if len(cur_dict['skin and plasma']) > 5:
            plt.text(6.5, 3.64, "           ".join(cur_dict['skin and plasma'][0:6]),fontsize=24, verticalalignment='top')
            plt.text(6.5, 2.90, "           ".join(cur_dict['skin and plasma'][6:]),fontsize=24, verticalalignment='top')
        else:
            plt.text(6.5, 3.54, "           ".join(cur_dict['skin and plasma']),fontsize=24, verticalalignment='top')
    if (cur_dict['muscle, skin, and plasma']) != []:
        plt.plot(x_coords_edge, msp_y1, color, linewidth = thickness)
        plt.plot(x_coords_edge, msp_y2, color, linewidth = thickness)
        plt.plot(x_coords_edge, msp_y3, color, linewidth = thickness)
        if len(cur_dict['muscle and plasma']) > 5:
            plt.text(6.5, 1.97, "           ".join(cur_dict['muscle, skin, and plasma'][0:6]),fontsize=24, verticalalignment='top')
            plt.text(6.5, 1.67, "           ".join(cur_dict['muscle, skin, and plasma'][6:]),fontsize=24, verticalalignment='top')
        else:
            plt.text(6.5, 1.77, "           ".join(cur_dict['muscle, skin, and plasma']),fontsize=24, verticalalignment='top')
        

    plt.xlim(0,10)
    plt.ylim(-1,12)
    plt.axis('off')
    return figure


# Function to generate hypergraph image for a dynamic time interval
def generate_dynamic_hypergraphs(cur_graph, title):
    "cur_graph: a data frame containing 6 columns organized as --> cytokines in muscle, pearson's r muscle, etc."
    "title: title for dynamic hypergraph image"
    edge_neg_095 ={'muscle':[], 'skin':[], 'plasma':[], 'muscle and skin':[], 'muscle and plasma':[], 'skin and plasma':[], 'muscle, skin, and plasma':[]}
    edge_neg_07 = {'muscle':[], 'skin':[], 'plasma':[], 'muscle and skin':[], 'muscle and plasma':[], 'skin and plasma':[],'muscle, skin, and plasma':[]}
    edge_pos_095 = {'muscle':[], 'skin':[], 'plasma':[], 'muscle and skin':[], 'muscle and plasma':[], 'skin and plasma':[],'muscle, skin, and plasma':[]}
    edge_pos_07 = {'muscle':[], 'skin':[], 'plasma':[], 'muscle and skin':[], 'muscle and plasma':[], 'skin and plasma':[], 'muscle, skin, and plasma':[]}
    
    # iterate through each row of cur_graph     
    for index, row in cur_graph.iterrows():
        muscle_sig = row['Pearsons R - Muscle']
        skin_sig = row ['Pearsons R - Skin']
        plasma_sig = row['Pearsons R - Plasma']
        # Dependeing on the value of r, append the cytokine in muscle, skin, or plasma into its
        # respective dictionary of edges
        if math.isnan(muscle_sig) == False:
            if muscle_sig == 0.7:
                edge_pos_07['muscle'].append(row['Cytokines in Muscle'])
            elif muscle_sig == 0.95:
                edge_pos_095['muscle'].append(row['Cytokines in Muscle'])
            elif muscle_sig == -0.7:
                edge_neg_07['muscle'].append(row['Cytokines in Muscle'])
            elif muscle_sig == -0.95:
                edge_neg_095['muscle'].append(row['Cytokines in Muscle'])
        if math.isnan(skin_sig) == False:
            if skin_sig == 0.7:
                edge_pos_07['skin'].append(row['Cytokines in Skin'])
            elif skin_sig == 0.95:
                edge_pos_095['skin'].append(row['Cytokines in Skin'])
            elif skin_sig == -0.7:
                edge_neg_07['skin'].append(row['Cytokines in Skin'])
            elif skin_sig == -0.95:
                edge_neg_095['skin'].append(row['Cytokines in Skin'])
        if math.isnan(plasma_sig) == False:
            if plasma_sig == 0.7:
                edge_pos_07['plasma'].append(row['Cytokines in Plasma'])
            elif plasma_sig == 0.95:
                edge_pos_095['plasma'].append(row['Cytokines in Plasma'])
            elif plasma_sig == -0.7:
                edge_neg_07['plasma'].append(row['Cytokines in Plasma'])
            elif plasma_sig == -0.95:
                edge_neg_095['plasma'].append(row['Cytokines in Plasma'])

    edge_neg_095_sorted = get_hypergraph_combined_edges(edge_neg_095)
    edge_neg_07_sorted = get_hypergraph_combined_edges(edge_neg_07)
    edge_pos_095_sorted = get_hypergraph_combined_edges(edge_pos_095)
    edge_pos_07_sorted = get_hypergraph_combined_edges(edge_pos_07)
    
    fig = plt.figure(figsize=(18, 18))
    gs = fig.add_gridspec(nrows=2, ncols=1, hspace= 0.5, wspace=1.5)
    
    fig = plot_hypergraph(edge_pos_095_sorted, 'k-', fig, gs, 'A', "Pearson's r > +0.95", 8)
#     fig = plot_hypergraph(edge_pos_07_sorted, 'k-', fig, gs, 'B', "Pearson's r > +0.7", 2)
    fig = plot_hypergraph(edge_neg_095_sorted, 'r-', fig, gs, 'C', "Pearson's r < - 0.95", 8)
#     fig = plot_hypergraph(edge_neg_07_sorted, 'r-', fig, gs, 'D', "Pearson's r < -0.7", 2)
    fig.suptitle(title, size=40)
    fig.savefig('A3A4_%s.png'%title, bbox_inches="tight") #CHANGE FILE NAME

#Function to get dynamic hypergraphs, in tabular form, across all dynamic time intervals
#Saves dynamic hypergraphs to a multi-tab excel file
def get_all_dynamic_hypergraphs_EXCEL(MUSCLE, SKIN, PLASMA, cytokines):
    "Inputs: median cytokine values at each time point in muscle, skin, and plasma"
    "cytokines: a list of cytokines"
    all_times_str = ['d0', 'd3','d5','d7','d9','d11','d20','d23','d25','d27','d31']
    all_times_int = [0, 3, 5, 7, 9, 11, 20, 23, 25, 27, 31]
    writer = pandas.ExcelWriter('Dynamic_Hypergraphs_Grouped_edges_A9_095.xlsx') #CHANGE FILE NAME
    for n in range (0,9):
        cur_times_str = all_times_str[n:n+3]
        cur_times_int = all_times_int [n:n+3]
        cur_dynamic_hypergraph = get_dynamic_interval_data(MUSCLE, SKIN, PLASMA, cur_times_str, cur_times_int, cytokines)   
        cur_dynamic_hypergraph = cur_dynamic_hypergraph.reset_index()
        cur_dynamic_hypergraph = hypergraphs_grouped_edges(cur_dynamic_hypergraph)
        cur_dynamic_hypergraph.to_excel(writer, '_'.join(cur_times_str))
    writer.save()
    
#Function to get dynamic hypergraphs, in image form, across all dynamic time intervals
#Saves dynamic hypergraph images to local folder
def get_all_dynamic_hypergraphs_IMGS(MUSCLE, SKIN, PLASMA, cytokines):
    all_times_str = ['d0', 'd3','d5','d7','d9','d11','d20','d23','d25','d27','d31']
    all_times_int = [0, 3, 5, 7, 9, 11, 20, 23, 25, 27, 31]
    for n in range (0,9):
        cur_times_str = all_times_str[n:n+3]
        cur_times_int = all_times_int [n:n+3]
        cur_dynamic_hypergraph = get_dynamic_interval_data(MUSCLE, SKIN, PLASMA, cur_times_str, cur_times_int, cytokines)   
        cur_dynamic_hypergraph = cur_dynamic_hypergraph.reset_index()
        cur_dynamic_hypergraph = generate_dynamic_hypergraphs(cur_dynamic_hypergraph, ', '.join(cur_times_str))
    
get_all_dynamic_hypergraphs_IMGS(MUSCLE, SKIN, PLASMA, cytokines)








# In[ ]:





# In[ ]:




