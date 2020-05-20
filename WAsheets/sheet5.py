# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 17:10:28 2020

@author: ntr002
"""
import os
import sys
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import datetime
import calendar
from matplotlib.colors import LinearSegmentedColormap
from dateutil.relativedelta import relativedelta
from . import GIS_functions as gis
from . import calculate_flux as cf
from . import get_dictionaries as gd
from . import hydroloop as hl

def main(BASIN,unit_conversion=1):
    '''
    unit_conversion: 1 for TCM, 1000 for MCM, 1e6 for BCM or km3)
    '''
    ### check requirements
    requirements=gd.get_requirements_for_sheet(5)    
    if not cf.check_requirement_sheet(BASIN,requirements):
        sys.exit("ERROR: Data requirements for Sheet 3 are not fulfilled")
    #create folder to save intermetidate data
    folder=os.path.join(BASIN['output_folder'],
                        'csv','timeseries')        
    if not os.path.exists(folder):
        os.makedirs(folder)
    output_file=os.path.join(folder,'sheet5_{0}.csv')
    #create folder to save sheet 3 csv
    sheet_folder=os.path.join(BASIN['output_folder'],'csv','sheet5') 
    if not os.path.exists(folder):
        os.makedirs(sheet_folder) 
    ### Get LU classes dictionary
    lu_dict=gd.get_sheet1_classes()  
    ### Calculate fraction to split outflow
    # non-utilizable fraction
    # non-recoerable fraction
    # commited fraction  
    
    ### Calulate yearly data to fill in Sheet 5  
    data={'basin':dict()}
    # Full basin
    data['basin']['sro']=cf.calc_flux_per_basin(
             BASIN['data_cube']['sro_yearly'], 
             BASIN['gis_data']['lu_map'], 
             BASIN['gis_data']['basin_mask'],
             chunksize=BASIN['chunksize'], 
             output=output_file.format('Basin_SRO'), 
             lu_dictionary=lu_dict, #calc for LU cat
             quantity='volume')
    data['basin']['bf']=cf.calc_flux_per_basin(
             BASIN['data_cube']['sro_yearly'], 
             BASIN['gis_data']['lu_map'], 
             BASIN['gis_data']['basin_mask'],
             chunksize=BASIN['chunksize'], 
             output=output_file.format('Basin_BF'), 
             lu_dictionary=lu_dict, #calc for LU cat
             quantity='volume')
    data['basin']['supply']=cf.calc_flux_per_basin(
             BASIN['data_cube']['sro_yearly'], 
             BASIN['gis_data']['lu_map'], 
             BASIN['gis_data']['basin_mask'],
             chunksize=BASIN['chunksize'], 
             output=output_file.format('Basin_Supply'),
             quantity='volume')
    data['basin']['sroincr']=cf.calc_flux_per_basin(
             BASIN['data_cube']['sro_yearly'], 
             BASIN['gis_data']['lu_map'], 
             BASIN['gis_data']['basin_mask'],
             chunksize=BASIN['chunksize'], 
             output=output_file.format('Basin_SROincr'),
             quantity='volume')
    data['basin']['percincr']=cf.calc_flux_per_basin(
             BASIN['data_cube']['sro_yearly'], 
             BASIN['gis_data']['lu_map'], 
             BASIN['gis_data']['basin_mask'],
             chunksize=BASIN['chunksize'], 
             output=output_file.format('Basin_PERCincr'),
             quantity='volume')
 
    #Sub basins
    for ID in BASIN['gis_data']['subbasin_mask'].keys():
        data[ID]=dict()
        data[ID]['sro']=cf.calc_flux_per_basin(
             BASIN['data_cube']['sro_yearly'], 
             BASIN['gis_data']['lu_map'], 
             BASIN['gis_data']['subbasin_mask'][ID],
             chunksize=BASIN['chunksize'], 
             output=output_file.format(
                     'Subbasin_{0}_SRO'.format(ID)), 
             lu_dictionary=lu_dict, #calc for LU cat
             quantity='volume')
        data[ID]['bf']=cf.calc_flux_per_basin(
             BASIN['data_cube']['sro_yearly'], 
             BASIN['gis_data']['lu_map'], 
             BASIN['gis_data']['subbasin_mask'][ID],
             chunksize=BASIN['chunksize'], 
             output=output_file.format(
                     'Subbasin_{0}_BF'.format(ID)),  
             lu_dictionary=lu_dict, #calc for LU cat
             quantity='volume')
        data[ID]['supply']=cf.calc_flux_per_basin(
             BASIN['data_cube']['sro_yearly'], 
             BASIN['gis_data']['lu_map'], 
             BASIN['gis_data']['subbasin_mask'][ID],
             chunksize=BASIN['chunksize'], 
             output=output_file.format(
                     'Subbasin_{0}_Supply'.format(ID)), 
             quantity='volume')
        data[ID]['sroincr']=cf.calc_flux_per_basin(
             BASIN['data_cube']['sro_yearly'], 
             BASIN['gis_data']['lu_map'], 
             BASIN['gis_data']['subbasin_mask'][ID],
             chunksize=BASIN['chunksize'], 
             output=output_file.format(
                     'Subbasin_{0}_SROincr'.format(ID)), 
             quantity='volume')
        data[ID]['percincr']=cf.calc_flux_per_basin(
             BASIN['data_cube']['sro_yearly'], 
             BASIN['gis_data']['lu_map'], 
             BASIN['gis_data']['subbasin_mask'][ID],
             chunksize=BASIN['chunksize'], 
             output=output_file.format(
                     'Subbasin_{0}_PERCincr'.format(ID)), 
             quantity='volume')  
    #connection between subbasin
    dico_in = BASIN['param']['dico_in']
    dico_out = BASIN['param']['dico_out']
    ### Fill data in Sheet 5 csv
    for i in range(data['basin']['sro']):
        year=data['basin']['sro'].index[i]
        results=Vividict()
        for subbasin in data.keys():
            for j in len(lu_dict.keys()):
                results['surf_runoff'][subbasin][lu_dict.keys()[j]]=\
                data[subbasin]['sro'][data[subbasin]['sro'].columns[j]].values[i] 
                results['base_runoff'][subbasin][lu_dict.keys()[j]]=\
                data[subbasin]['bf'][data[subbasin]['bf'].columns[j]].values[i]
            results['base_runoff'][subbasin]                 
                

        #write sheet 5 csv
        output_fh=os.path.join(sheet_folder,'sheet5_{0}.csv'.format(year))
        create_sheet5_csv(year,results,output_fh)    
        
    return True

class Vividict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value        

def create_sheet5_csv(dresults, output_fh):
    """
    Create the csv-file for sheet 5.

    Parameters
    ----------
    results : dict
        Dictionary of results generated in sheet5_run.py
    output_fh : str
        Filehandle to store the csv-file.
    """
    first_row = ['SUBBASIN', 'VARIABLE', 'VALUE', 'UNITS']
    if not os.path.exists(os.path.split(output_fh)[0]):
        os.makedirs(os.path.split(output_fh)[0])
    csv_file = open(output_fh, 'w')
    writer = csv.writer(csv_file, delimiter=';', lineterminator = '\n')
    writer.writerow(first_row)
    lu_classes = ['PROTECTED', 'UTILIZED', 'MODIFIED', 'MANAGED']
    for sb in list(dresults['surf_runoff'].keys()):
        writer.writerow([sb, 'Inflow', '{0}'.format(dresults['inflows'][sb]), 'km3'])
        for lu_class in lu_classes:
            writer.writerow([sb, 'Fast Runoff: '+lu_class, '{0}'.format(dresults['surf_runoff'][sb][lu_class]), 'km3'])
            writer.writerow([sb, 'Slow Runoff: ' +lu_class, '{0}'.format(dresults['base_runoff'][sb][lu_class]), 'km3'])
        writer.writerow([sb, 'Total Runoff', '{0}'.format(dresults['total_runoff'][sb]), 'km3'])
        writer.writerow([sb, 'SW withdr. manmade', '{0}'.format(dresults['withdrawls'][sb]['man']), 'km3'])
        writer.writerow([sb, 'SW withdr. natural', '{0}'.format(dresults['withdrawls'][sb]['natural']), 'km3'])
        writer.writerow([sb, 'SW withdr. total', '{0}'.format(dresults['withdrawls'][sb]['man']+dresults['withdrawls'][sb]['natural']), 'km3'])
        writer.writerow([sb, 'Return Flow SW', '{0}'.format(dresults['return_sw_sw'][sb]), 'km3'])
        writer.writerow([sb, 'Return Flow GW', '{0}'.format(dresults['return_gw_sw'][sb]), 'km3'])
        writer.writerow([sb, 'Total Return Flow', '{0}'.format(dresults['return_sw_sw'][sb]+dresults['return_gw_sw'][sb]), 'km3'])
        writer.writerow([sb, 'Outflow: Total', '{0}'.format(dresults['total_outflow'][sb]), 'km3'])
        writer.writerow([sb, 'Outflow: Committed', 'nan', 'km3'])
        writer.writerow([sb, 'Outflow: Non Recoverable', 'nan', 'km3'])
        writer.writerow([sb, 'Outflow: Non Utilizable', 'nan', 'km3'])
        writer.writerow([sb, 'Outflow: Utilizable', 'nan', 'km3'])
        writer.writerow([sb,'Interbasin Transfer','{0}'.format(dresults['interbasin_transfers'][sb]),'km3'])
        writer.writerow([sb, 'SW storage change', '{0}'.format(dresults['deltaS'][sb]), 'km3'])
    csv_file.close()
    return

