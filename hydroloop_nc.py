# -*- coding: utf-8 -*-
"""

"""
from WAsheets import *
from WAsheets.sheet1 import main as sheet1
from WAsheets.sheet2 import main as sheet2

BASIN={
       'name': 'Kanartaka',
       'hydroyear':'A-DEC', #Water year end month
       'chunksize':None,
       'output_folder':r"D:\Projects\ADB\analysis\Karnataka\Results_Karnataka\Hydroloop",
       'gis_data':{
               'basin_mask': r"D:\Projects\ADB\analysis\Karnataka\Results_Karnataka\SM_balance\K2_K3_K4\Tot_Mask_WGS84.tif",
               'subbasin_mask':{
                       1: r".tif",
                       2: r".tif"
                       },
               'lu_map': r"D:\Projects\ADB\analysis\Karnataka\Results_Karnataka\SM_balance\K2_K3_K4\K2K3K4_LU_WA.nc",
               'dem':None,
               },
       'data_cube':{
               'p_monthly':r"D:\Projects\ADB\analysis\Karnataka\Results_Karnataka\SM_balance\K2_K3_K4\K2K3K4_P_GPM.nc",
               'et_monthly':r"D:\Projects\ADB\analysis\Karnataka\Results_Karnataka\SM_balance\K2_K3_K4\K2K3K4_ETa_SSEBop.nc",
               'etincr_monthly':r"D:\Projects\ADB\analysis\Karnataka\Results_Karnataka\SM_balance\K2_K3_K4\etincr_monthly.nc",
               'etrain_monthly':r"D:\Projects\ADB\analysis\Karnataka\Results_Karnataka\SM_balance\K2_K3_K4\etrain_monthly.nc",
               'lai_monthly':r"D:\Projects\ADB\analysis\Karnataka\Results_Karnataka\SM_balance\K2_K3_K4\K2K3K4_LAI_MOD15.nc",
               
               },
        'ts_data':{
                'q_in_sw':None,
                'q_in_gw':None,
                'q_in_desal':None,
                'q_outflow':None,
                'q_out_sw':None,
                                
                },
        'param':{
            'crops':[],
            'non_crop':[],
            'dico_in':{1:[], 2:[1]},
            'dico_out':{1:[], 2:[0]},
            'fraction_xs':[10, 50, 10, 50],
        }
        
       }
       
sheet1(BASIN,unit_conversion=1000)
sheet2(BASIN,unit_conversion=1000)
       