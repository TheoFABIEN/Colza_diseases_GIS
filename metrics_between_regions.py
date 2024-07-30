import os
import tkinter as tk
from tkinter import filedialog

import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import shutil

import utils 


folders_path = filedialog.askdirectory()

if os.path.exists(folders_path + '/0general_plots'):
    shutil.rmtree(folders_path + '/0general_plots')

def find_files(directory):
    metrics_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.startswith('metrics') and not file.endswith('.png'):
                #yield os.path.join(root, file)
                metrics_files.append(os.path.join(root, file))
    return metrics_files


def list_subdirectories(directory):
    return [
        os.path.join(directory, d) 
        for d in os.listdir(directory) 
        if os.path.isdir(os.path.join(directory, d))
    ]


total_num_colza = []
total_surf_colza = []
total_prop_colza = []
total_surf_ann = []
total_surf_per = []
total_surf_unknown = []



sub_folders = list_subdirectories(folders_path)
region_list = [0]* len(sub_folders)

for idx, region_dir in enumerate(sub_folders):

    num_colza_region = []
    surf_colza_region = []
    prop_colza_region = []
    surf_ann_region = []
    surf_per_region = []
    surf_unknown_region = []
    

    for file_path in find_files(region_dir):

        if file_path is not None:
            region_list[idx] = os.path.basename(region_dir)

        metrics = pd.read_csv(file_path, index_col = 0)

        num_colza = metrics.loc['Number_colza_fields']
        surf_colza = metrics.loc['Total_surface_colza']
        prop_colza = metrics.loc['Proportion of fields that contain colza']
        surf_ann = metrics.loc['Surface_Ann']
        surf_per = metrics.loc['Surface_Per']
        surf_unknown = metrics.loc['Surface_unknown']
        

        if num_colza is not None:
            num_colza_region.extend(num_colza.tolist())
        if surf_colza is not None:
            surf_colza_region.extend(surf_colza.tolist())
        if prop_colza is not None:
            prop_colza_region.extend(prop_colza.tolist())
        if surf_ann is not None:
            surf_ann_region.extend(surf_ann.tolist())
        if surf_per is not None:
            surf_per_region.extend(surf_per.tolist())
        if surf_unknown is not None:
            surf_unknown_region.extend(surf_unknown.tolist())



    total_num_colza.append(num_colza_region)
    total_surf_colza.append(surf_colza_region)
    total_prop_colza.append(prop_colza_region)
    total_surf_ann.append(surf_ann_region)
    total_surf_per.append(surf_per_region)
    total_surf_unknown.append(surf_unknown_region)
    
region_list = [s.replace('_', ' ') for s in region_list]



####### Create and save figure ##############

os.mkdir(folders_path+ '/0general_plots')

plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False

fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(
    2, 3, figsize  = (12, 7)
)
fig.subplots_adjust(bottom= .25, top = .93, hspace = .45, wspace = .4) 
ax1.boxplot(total_num_colza)
ax2.boxplot(total_surf_colza)
ax3.boxplot(total_prop_colza)
ax4.boxplot(total_surf_ann)
ax5.boxplot(total_surf_per)
ax6.boxplot(total_surf_unknown)

ax1.set_xticklabels('', rotation = 45, ha = 'right')
ax2.set_xticklabels('', rotation = 45, ha = 'right')
ax3.set_xticklabels('', rotation = 45, ha = 'right')
ax4.set_xticklabels(region_list, rotation = 45, ha = 'right')
ax5.set_xticklabels(region_list, rotation = 45, ha = 'right') 
ax6.set_xticklabels(region_list, rotation = 45, ha = 'right') 
ax1.set_title('Number of colza fields')
ax2.set_title('Total surface of colza')
ax3.set_title('Proportion of fields that contain colza (or mustard)')
ax4.set_title('Total surface of annual')
ax5.set_title('Total surface of permanent')
ax6.set_title('Total surface of unknown')


fig.savefig(folders_path + '/0general_plots/metrics_between_regions.png')




########### Time instead of regions #########


num_years = len(metrics.columns) # Use the last dataframe we opened
all_metrics_files = find_files(folders_path)

total_num_colza = np.zeros((len(all_metrics_files), num_years))
total_surf_colza = np.zeros((len(all_metrics_files), num_years))
total_surf_ann = np.zeros((len(all_metrics_files), num_years))
total_surf_per = np.zeros((len(all_metrics_files), num_years))
total_surf_unknown = np.zeros((len(all_metrics_files), num_years))
total_area_bound = np.zeros((len(all_metrics_files), num_years))


for idx, metrics_file in enumerate(all_metrics_files):

    metrics_df = pd.read_csv(metrics_file, index_col = 0)

    for t, year in enumerate(metrics_df.columns):

        numcol = metrics_df.loc['Number_colza_fields', year]
        surfcol = metrics_df.loc['Total_surface_colza', year]
        surfann = metrics_df.loc['Surface_Ann', year]
        surfper = metrics_df.loc['Surface_Per', year]
        surfunknown = metrics_df.loc['Surface_unknown', year]
        boundarea = metrics_df.loc['Area_boundary_zones_with_prev_year', year]

        if numcol is not None:
            total_num_colza[idx, t] = numcol
            total_surf_colza[idx, t] = surfcol
            total_surf_ann[idx, t] = surfann
            total_surf_per[idx, t] = surfper
            total_surf_unknown[idx, t] = surfunknown
            total_area_bound[idx, t] = boundarea 



fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(
    2, 3, figsize  = (12, 7)
)
fig.subplots_adjust(bottom= .1, top = .93, hspace = .45, wspace = .4) 
ax1.boxplot(total_num_colza)
ax2.boxplot(total_surf_colza)
ax3.boxplot(total_surf_ann)
ax4.boxplot(total_surf_per)
ax5.boxplot(total_surf_unknown)
ax6.boxplot(total_area_bound)

year_names = [f"Year {x}" for x in range(1, num_years + 1)]

ax1.set_xticklabels('')
ax2.set_xticklabels('')
ax3.set_xticklabels('')
ax4.set_xticklabels(year_names, rotation = 45, ha = 'right')
ax5.set_xticklabels(year_names, rotation = 45, ha = 'right') 
ax6.set_xticklabels(year_names, rotation = 45, ha = 'right') 
ax1.set_title('Number of colza fields')
ax2.set_title('Total surface of colza')
ax3.set_title('Total surface of annual')
ax4.set_title('Total surface of permanent')
ax5.set_title('Total surface of unknown')
ax6.set_title('Total surface of meeting zones between t and t+1')

fig.savefig(folders_path + '/0general_plots/metrics_between_years.png')



















