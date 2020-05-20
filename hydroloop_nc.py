# -*- coding: utf-8 -*-
"""

"""
import os
from . import calculate_flux as cf
from . import hydroloop as hl

wp_folder=r"D:\Projects\ADB\analysis\Karnataka\Results_Karnataka\SM_balance"
BASIN={
       'name': 'Kanartaka',
       'hydroyear':'A-DEC', #Water year end month
       'chunksize':None,
       'output_folder':r"D:\Projects\ADB\analysis\Karnataka\Results_Karnataka\Hydroloop",
       'gis_data':{
               'basin_mask': os.path.join(wp_folder,
                                          'K2_K3_K4','Tot_Mask_WGS84.tif'),
               'subbasin_mask':{
                       1: r".tif",
                       2: r".tif"
                       },
               'lu_map': os.path.join(wp_folder,
                                          'K2_K3_K4','K2K3K4_LU_WA.nc'),
               'dem':None,
               },
       'data_cube':{
               'p_monthly':os.path.join(wp_folder,
                                          'K2_K3_K4','K2K3K4_P_GPM.nc'),
               'et_monthly':os.path.join(wp_folder,
                                          'K2_K3_K4','K2K3K4_ETa_SSEBop.nc'),
               'etincr_monthly':os.path.join(wp_folder,
                                          'K2_K3_K4','etincr_monthly.nc'),
               'etrain_monthly':os.path.join(wp_folder,
                                          'K2_K3_K4','etrain_monthly.nc'),
               'lai_monthly':os.path.join(wp_folder,
                                          'K2_K3_K4','K2K3K4_LAI_MOD15.nc'),
               
               
               },
        'ts_data':{
                'q_in_sw':{
                        'basin':r"total_inflow.csv",
                        1:r"inflow_subbasin.csv",
                        2:r"inflow_subbasin.csv",
                        },
                'q_in_gw':None,
                'q_in_desal':None,
                'q_outflow':None,
                'q_out_sw':{
                        'basin':r"interbasin_transfer.csv",
                        1:r"interbasin_transfer_subbasin.csv",
                        2:r"interbasin_transfer_subbasin.csv",
                        },
                                
                },
        'params':{
            'crops':[
                    (r'seasons_Kharif.csv', 
                     'Kharif', 'Rainfed', 'Cereals', 35.0) ,
                    (r'seasons_Kharif.csv', 
                     'Kharif', 'Irrigated', 'Cereals', 54.0),                                       
                    (r'seasons_Rabi.csv', 
                     'Rabi', 'Rainfed', 'Root/tuber', 36.0),                                       
                    (r'seasons_Rabi.csv', 
                     'Rabi', 'Irrigated', 'Root/tuber', 55.0),                                                            
                    (r'seasons_Zaid.csv', 
                     'Zaid', 'Rainfed', 'Leguminous', 37.0),                                       
                    (r'seasons_Zaid.csv', 
                     'Zaid', 'Irrigated', 'Leguminous', 56.0),                                       
                    (r'seasons_double_triple.csv', 
                     'Double/Triple Crop', 'Irrigated', 'Sugar', 57.0),                                       
                    (r'seasons_33.csv', 
                     'Forest plantation', 'Rainfed', '-', 33.0),                                       
                    (r'seasons_52.csv', 
                     'Forest plantation', 'Irrigated', '-', 52.0),                                                            
                    ],
            'wcpc':110, #Water consumption per capita per day in [liter/person/day]
            'dico_in':{1:[], 2:[1]},
            'dico_out':{1:[], 2:[0]},
        }
        
       }
       
sheet1(BASIN,unit_conversion=1000)
sheet2(BASIN,unit_conversion=1000)


#%% Test discharge

bf_nc=r"D:\Projects\ADB\analysis\Karnataka\Results_Karnataka\New_SMBalance\K4\bf_monthly.nc"
sroincr_nc=r"D:\Projects\ADB\analysis\Karnataka\Results_Karnataka\New_SMBalance\K4\d_sro_monthly.nc"
sro_nc=r"D:\Projects\ADB\analysis\Karnataka\Results_Karnataka\New_SMBalance\K4\sro_monthly.nc"
supply_nc=r"D:\Projects\ADB\analysis\Karnataka\Results_Karnataka\New_SMBalance\K4\supply_monthly.nc"

basin_mask=r"D:\Projects\ADB\analysis\Karnataka\Results_Karnataka\New_SMBalance\K4\Shape\K4.tif"

df_bf=calc_flux_per_basin(bf_nc, basin_mask,
                          chunksize=None,output=None,quantity='volume')
df_sro=calc_flux_per_basin(sro_nc, basin_mask,
                          chunksize=None,output=None,quantity='volume')
df_sroincr=calc_flux_per_basin(sroincr_nc, basin_mask,
                          chunksize=None,output=None,quantity='volume')
df_supply=calc_flux_per_basin(supply_nc, basin_mask,
                          chunksize=None,output=None,quantity='volume')
df_supply=df_supply*0.5
discharge,dS_sw=calc_sw_from_wp(df_sro,df_sroincr,df_bf,df_supply,plot=True)

plt.plot(discharge,'k',label='discharge')
plt.plot(dS_sw,'orange',label='dS_sw')
plt.legend()

