# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 17:06:08 2020

@author: ntr002
"""
import os
import csv
import pandas as pd
from . import calculate_flux as cf
from . import get_dictionaries as gd
from . import hydroloop as hl

def main(BASIN,unit_conversion=1000):
    '''
    unit_conversion: 1 for TCM, 1000 for MCM, 1e6 for BCM or km3)
    '''
    #check requirements
#    requirements=gd.get_sheet_requirements(1)    
#    if not cf.check_requirement_sheet(BASIN,requirements):
#        print("ERROR: Data requirements for Sheet 1 are not fulfilled")
#        return None
    #get sheet 1 dictionary
    lu_dictionary=gd.get_sheet1_classes()  
    
    folder=os.path.join(BASIN['output_folder'],'csv','timeseries')        
    if not os.path.exists(folder):
        os.makedirs(folder)
    output_file=os.path.join(folder,'sheet1_{0}.csv')
    
    #Calulate yearly data to fill in Sheet 1
    df_P=cf.calc_flux_per_basin(BASIN['data_cube']['monthly']['p'], 
                                BASIN['gis_data']['basin_mask'],
                                chunksize=BASIN['chunksize'],
                                output=output_file.format('basin_p_monthly'),
                                quantity='volume')
    df_ET=cf.calc_flux_per_basin(BASIN['data_cube']['monthly']['et'], 
                                BASIN['gis_data']['basin_mask'], 
                                chunksize=BASIN['chunksize'],
                                output=output_file.format('basin_p_monthly'),
                                quantity='volume')
    df_ETrain=cf.calc_flux_per_LU_class(BASIN['data_cube']['monthly']['etrain'], 
                             BASIN['gis_data']['lu_map'], 
                             BASIN['gis_data']['basin_mask'],
                     chunksize=BASIN['chunksize'], #option to process in chunks
                     output=output_file.format('basin_etrain_monthly'), 
                     #option to save output as csv                     
                     lu_dictionary=lu_dictionary, #calc for LU categories
                     quantity='volume')
    
    df_ETincr=cf.calc_flux_per_LU_class(BASIN['data_cube']['monthly']['etincr'], 
                             BASIN['gis_data']['lu_map'], 
                             BASIN['gis_data']['basin_mask'],
                     chunksize=BASIN['chunksize'], #option to process in chunks
                     output=output_file.format('basin_etincr_monthly'), 
                     #option to save output as csv                     
                     lu_dictionary=lu_dictionary, #calc for LU categories
                     quantity='volume')
    
    sheet_folder=os.path.join(BASIN['output_folder'],'csv','sheet1') 
    if not os.path.exists(sheet_folder):
        os.makedirs(sheet_folder) #create sheet1 folder       
    #Fill in Sheet 1 csv
    monthly_csvs=[]
    for i in range(len(df_P)):
        results=dict()
        results['p_advection']=df_P[df_P.columns[0]].values[i]/unit_conversion
        results['landscape_et_plu']=df_ETrain[df_ETrain.columns[0]].values[i]/unit_conversion
        results['landscape_et_ulu']=df_ETrain[df_ETrain.columns[1]].values[i]/unit_conversion
        results['landscape_et_mlu']=df_ETrain[df_ETrain.columns[2]].values[i]/unit_conversion
        results['landscape_et_mwu']=df_ETrain[df_ETrain.columns[3]].values[i]/unit_conversion
        
        results['uf_plu']=df_ETincr[df_ETincr.columns[0]].values[i]/unit_conversion
        results['uf_ulu']=df_ETincr[df_ETincr.columns[1]].values[i]/unit_conversion
        results['uf_mlu']=df_ETincr[df_ETincr.columns[2]].values[i]/unit_conversion
        results['uf_mwu']=df_ETincr[df_ETincr.columns[3]].values[i]/unit_conversion
        
        results['manmade']=df_ETincr[df_ETincr.columns[3]].values[i]/unit_conversion
        results['natural']=(df_ETincr[df_ETincr.columns[0]].values[i]
        +df_ETincr[df_ETincr.columns[1]].values[i]
        +df_ETincr[df_ETincr.columns[2]].values[i])/unit_conversion
        #year value
        year=df_P.index[i].year
        month=df_P.index[i].month
        #Read time-series input
        ts_data_sheet1=['q_in_sw', 'q_in_gw', 'q_in_desal',
                        'q_outflow','q_out_sw','q_out_gw']
        for key in ts_data_sheet1:
            if BASIN['ts_data'][key]['basin'] is not None:
                df=pd.read_csv(BASIN['ts_data']['basin'][key],sep=';',index_col=0)
                results[key]=df[df.columns[0]].values[i]
            else:
                results[key]=0.
        #calulate water balance
       
        P=results['p_advection']
        ET=df_ET[df_ET.columns[0]].values[i]/unit_conversion
        Qin=results['q_in_sw']+results['q_in_gw']+results['q_in_desal']
        Qout=results['q_outflow']+results['q_out_sw']+results['q_out_gw']
        #In case yearly dS is not available, calculate dS
        if 'dS' not in BASIN['ts_data'].keys():
            results['dS']=calc_water_balance_residual(P,ET,Qin=Qin,Qout=Qout)
        #in case yearly dS is available, calculate q_outflow
        else:
            df=pd.read_csv(BASIN['ts_data']['dS'],sep=';',index_col=0)
            results['dS']=df[df.columns[0]].values[i]            
            results['q_outflow']=calc_water_balance_residual(P,ET,
                   Qin=Qin,dS=results['dS'])
        #write sheet 1 csv
        output_fh=os.path.join(sheet_folder,'sheet1_{0}_{1}.csv'.format(year,month))
        create_sheet1_csv(results, output_fh) #write results to sheet1 csv
        monthly_csvs.append(output_fh)
    ##calculate yearly sheets
    yearly_folder=os.path.join(sheet_folder,'yearly') 
    if not os.path.exists(yearly_folder):
        os.makedirs(yearly_folder) #create sheet1 folder  
    yearly_csvs=hl.calc_yearly_sheet(monthly_csvs,
                                     yearly_folder,
                                     hydroyear=BASIN['hydroyear'])
    return yearly_csvs

def calc_water_balance_residual(P,ET,dS=None,Qin=None,Qout=None):
    if Qin is not None:
        gross_inflow=P+Qin
    else:
        gross_inflow=P
    if dS is not None:
        Qout=gross_inflow-dS
        return Qout
    if Qout is not None:
        dS=gross_inflow-Qout
        return dS
        
def create_sheet1_csv(results, output_fh):
    """
    Create the csv-file needed to plot sheet 1.
    
    Parameters
    ----------
    results : dict
        Dictionary generated by calc_sheet1.
    output_fh : str
        Filehandle to store the csv-file.
    """
    first_row = ['CLASS', 'SUBCLASS', 'VARIABLE', 'VALUE']
    
    if not os.path.exists(os.path.split(output_fh)[0]):
        os.makedirs(os.path.split(output_fh)[0])
    
    csv_file = open(output_fh, 'w')
    writer = csv.writer(csv_file, delimiter=';', lineterminator = '\n')
    writer.writerow(first_row)

    writer.writerow(['INFLOW', 'PRECIPITATION', 'Rainfall',
                     '{0}'.format(results['p_advection'])])
    writer.writerow(['INFLOW', 'PRECIPITATION', 'Snowfall',
                     0.])
    writer.writerow(['INFLOW', 'PRECIPITATION', 'Precipitation recycling',
                     0.])
    writer.writerow(['INFLOW', 'SURFACE WATER', 'Main riverstem',
                     '{0}'.format(results['q_in_sw'])])
    writer.writerow(['INFLOW', 'SURFACE WATER', 'Tributaries',
                     0.])
    writer.writerow(['INFLOW', 'SURFACE WATER', 'Utilized surface water',
                     0.])
    writer.writerow(['INFLOW', 'SURFACE WATER', 'Flood', 
                     0.])
    writer.writerow(['INFLOW', 'GROUNDWATER', 'Natural', 
                     '{0}'.format(results['q_in_gw'])])
    writer.writerow(['INFLOW', 'GROUNDWATER', 'Utilized',
                     0.])
    writer.writerow(['INFLOW', 'OTHER', 'Desalinized', 
                     '{0}'.format(results['q_in_desal'])])
    writer.writerow(['STORAGE', 'CHANGE', 'Surface storage', 
                     '{0}'.format(results['dS'])])
    writer.writerow(['STORAGE', 'CHANGE', 'Storage in sinks',
                     0.])
    writer.writerow(['OUTFLOW', 'ET LANDSCAPE', 'Protected',
                     '{0}'.format(results['landscape_et_plu'])])
    writer.writerow(['OUTFLOW', 'ET LANDSCAPE', 'Utilized',
                     '{0}'.format(results['landscape_et_ulu'])])
    writer.writerow(['OUTFLOW', 'ET LANDSCAPE', 'Modified',
                     '{0}'.format(results['landscape_et_mlu'])])
    writer.writerow(['OUTFLOW', 'ET LANDSCAPE', 'Managed',
                     '{0}'.format(results['landscape_et_mwu'])])
    writer.writerow(['OUTFLOW', 'ET UTILIZED FLOW', 'Protected',
                     '{0}'.format(results['uf_plu'])])
    writer.writerow(['OUTFLOW', 'ET UTILIZED FLOW', 'Utilized',
                     '{0}'.format(results['uf_ulu'])])
    writer.writerow(['OUTFLOW', 'ET UTILIZED FLOW', 'Modified',
                     '{0}'.format(results['uf_mlu'])])
    writer.writerow(['OUTFLOW', 'ET UTILIZED FLOW', 'Managed',
                     '{0}'.format(results['uf_mwu'])])        
    writer.writerow(['OUTFLOW', 'ET INCREMENTAL', 'Manmade',
                     '{0}'.format(results['manmade'])])  
    writer.writerow(['OUTFLOW', 'ET INCREMENTAL', 'Natural',
                     '{0}'.format(results['natural'])])
    writer.writerow(['OUTFLOW', 'SURFACE WATER', 'Main riverstem',
                     '{0}'.format(results['q_outflow'])])
    writer.writerow(['OUTFLOW', 'SURFACE WATER', 'Tributaries', 
                     0.])  
    writer.writerow(['OUTFLOW', 'SURFACE WATER', 'Utilized surface water',
                     0.])  
    writer.writerow(['OUTFLOW', 'SURFACE WATER', 'Flood',
                     0.])
    writer.writerow(['OUTFLOW', 'SURFACE WATER', 'Interbasin transfer',
                     '{0}'.format(results['q_out_sw'])])
    writer.writerow(['OUTFLOW', 'GROUNDWATER', 'Natural', 
                     '{0}'.format(results['q_out_gw'])]) 
    writer.writerow(['OUTFLOW', 'GROUNDWATER', 'Utilized', 
                     0.])
    writer.writerow(['OUTFLOW', 'OTHER', 'Non-utilizable',
                     0.])
    writer.writerow(['OUTFLOW', 'OTHER', 'Other', 
                     0.])
    writer.writerow(['OUTFLOW', 'RESERVED', 'Commited',
                     0.]) 
    writer.writerow(['OUTFLOW', 'RESERVED', 'Navigational',
                     0.]) 
    writer.writerow(['OUTFLOW', 'RESERVED', 'Environmental',
                     0.]) 
    
    csv_file.close()
    
