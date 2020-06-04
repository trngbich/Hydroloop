# -*- coding: utf-8 -*-
"""

"""
import os
import sys
sys.path.append(r'D:\Projects\ADB\analysis\Karnataka\Scripts')
import pandas as pd

from WAsheets import calculate_flux as cf
from WAsheets import hydroloop as hl

wp_folder=r"D:\Projects\ADB\analysis\Karnataka\SMbalance_results\K2"
BASIN={
       'name': 'K2',
       'hydroyear':'A-MAY', #Water year end month
       'chunksize':None,
       'unit_conversion':1e6, #1e3 for input data in MCM, 1e6 for input data in km3
       'output_folder':r"D:\Projects\ADB\analysis\Karnataka\Hydroloop_results\K2",
       'gis_data':{
               'basin_mask': os.path.join(wp_folder,
                                          'Shape','K2.tif'),
               'subbasin_mask':{
                       1: r"D:\Projects\ADB\analysis\Karnataka\Hydroloop_results\K2\1_subbasin.tif",
                       2: r"D:\Projects\ADB\analysis\Karnataka\Hydroloop_results\K2\2_subbasin.tif"
                       },
               'dem':None,
               'aeisw':r"D:\Projects\ADB\analysis\Karnataka\Hydroloop_results\K2\gmia_v5_aeisw_pct_aei.tif", #area equipped with surface water irrigation percentage
               'population':r"D:\Projects\ADB\analysis\Karnataka\Hydroloop_results\K2\Pop_K2_K3_K4_corrected.tif",
               },
       'data_cube':{
           'monthly':{
               'p':os.path.join(wp_folder,
                                          'K2_P_CHIRPS.nc'),
                'etref':os.path.join(wp_folder,
                                          'K2_ETref.nc'),
               'et':os.path.join(wp_folder,
                                          'K2_ETa_SSEBop.nc'),
               'i':os.path.join(wp_folder,
                                          'i_monthly.nc'),
               't':None,
               'e':None,
               'nrd':os.path.join(wp_folder,
                                          'nRD_monthly.nc'),
               'etincr':os.path.join(wp_folder,
                                          'etincr_monthly.nc'),
               'etrain':os.path.join(wp_folder,
                                          'etrain_monthly.nc'),
               'lai':os.path.join(wp_folder,
                                          'K2_LAI_MOD15.nc'),
              'ndm':os.path.join(wp_folder,
                                          'K2_DPM_PROBA-V.nc'),
             'sro':os.path.join(wp_folder,
                                         'sro_monthly.nc'),
             'sroincr':os.path.join(wp_folder,
                                          'd_sro_monthly.nc'),
             'perc':os.path.join(wp_folder,
                                          'perco_monthly.nc'),
             'percincr':os.path.join(wp_folder,
                                          'd_perco_monthly.nc'),
             'bf':os.path.join(wp_folder,
                                          'bf_monthly.nc'),
            'supply':os.path.join(wp_folder,
                                          'supply_monthly.nc')
               },
           'yearly':{
                'lu':os.path.join(wp_folder,
                                          'K2_LU_WA.nc'),
                   }      
                     },
        'ts_data':{
                'q_in_sw':{
                        'basin':None,
                        1:None, #unit MCM
                        2:None,
                        },
                'q_in_gw':{
                        'basin':None,
                        1:None, #unit MCM
                        2:None,
                        },
                'q_in_desal':{
                        'basin':None,
                        1:None, #unit MCM
                        2:None,
                        },
                'q_outflow':{ #river flow
                        'basin':None,
                        1:None,
                        2:None
                        },
                'q_out_sw':{ #interbasin transfer
                        'basin':None,
                        1:None,
                        2:None,
                        },
                'q_out_gw':{ #interbasin transfer
                        'basin':None,
                        1:None,
                        2:None,
                        },
                'dS_sw':{
                        'basin':None,
                        1:None, #unit MCM
                        2:None,
                        },
                                
                },
        'params':{
            'crops':{#lu_code: [r'seasons_dates.csv','Crop_type']
                    35.0: [r'crop_season\seasons_Kharif.csv','N/A'],
                    54.0: [r'crop_season\seasons_Kharif.csv','N/A'],
                    36.0: [r'crop_season\seasons_Rabi.csv','N/A'],
                    55.0: [r'crop_season\seasons_Rabi.csv','N/A'],
                    37.0: [r'crop_season\seasons_Zaid.csv','N/A'],
                    56.0: [r'crop_season\seasons_Zaid.csv','N/A'],
                    57.0: [r'crop_season\seasons_double_triple.csv','N/A'],
                    33.0: [r'crop_season\seasons_plantation.csv','N/A'],                    
                    52.0: [r'crop_season\seasons_plantation.csv','N/A'],                      
                    },
            'dico_in':{1:[], 2:[1]},
            'dico_out':{1:[], 2:[0]},
            'residential_sw_supply_fraction':0.6,
            'wcpc':110, #Water consumption per capita per day in [liter/person/day]
            'wcpc_min':100 #minimum demand
        }
        
       }

#%% Calculate hydroloop datacube

### Resample yearly LU to monthly netCDF
yearly_nc=BASIN['data_cube']['yearly']['lu']
sample_nc=BASIN['data_cube']['monthly']['p']
monthly_nc=cf.resample_to_monthly_dataset(yearly_nc, sample_nc,
                                start_month=0,
                                output=None,
                                chunksize=None) 
BASIN['data_cube']['monthly']['lu']=monthly_nc
### Split ETI
e_nc,i_nc,t_nc=hl.split_ETI(et_nc=BASIN['data_cube']['monthly']['et'],
                            i_nc=BASIN['data_cube']['monthly']['i'],
                            t_nc=BASIN['data_cube']['monthly']['t'],
                              p_nc=BASIN['data_cube']['monthly']['p'],
                              lai_nc=BASIN['data_cube']['monthly']['lai'],
                              nrd_nc=BASIN['data_cube']['monthly']['nrd'],
                              ndm_nc=BASIN['data_cube']['monthly']['ndm']
              )        
BASIN['data_cube']['monthly']['i']=i_nc
BASIN['data_cube']['monthly']['t']=t_nc
BASIN['data_cube']['monthly']['e']=e_nc
### Calculate Water Productivity
wp_nc=hl.flow_ratio(BASIN['data_cube']['monthly']['ndm'],
                    BASIN['data_cube']['monthly']['et'],
                    name='dry_mass_water_productivity',
                    attrs={
                            'units':'0.1kg/m3',
                            'source':'Dry-Mass Productivity (kg/ha/month) divided by \
                            Actual evapotranspiration (mm/month)',
                            'quantity':'Dry-Mass Water Productivity'
                            }
                    )
BASIN['data_cube']['monthly']['ndm_wp']=wp_nc
### split supply
sw_supply_fraction_nc=hl.calc_sw_supply_fraction_by_LU(BASIN['data_cube']['monthly']['lu'],
                                                       BASIN['gis_data']['aeisw'])
sw_supply_nc,gw_supply_nc=hl.split_flow(BASIN['data_cube']['monthly']['supply'],
              fraction_nc=sw_supply_fraction_nc)
### demand of land surface
demand_nc=hl.calc_land_surface_water_demand(BASIN['data_cube']['monthly']['lai'],
                                  BASIN['data_cube']['monthly']['etref'],
                                  BASIN['data_cube']['monthly']['p'],
                                  BASIN['data_cube']['monthly']['lu'])
### non-consumed supply or return flow
return_nc=hl.substract_flow(BASIN['data_cube']['monthly']['supply'],
                                  BASIN['data_cube']['monthly']['etincr'],
                                  name='return')
### split return by sroincr/total_incremental ratio
sw_return_fraction_nc=hl.calc_sw_return_fraction(
        BASIN['data_cube']['monthly']['sroincr'],
        BASIN['data_cube']['monthly']['percincr'])
sw_return_nc,gw_return_nc=hl.split_flow(return_nc,fraction_nc=sw_return_fraction_nc)

### residential supply and demand
residential_supply_nc=hl.calc_residential_water_consumption(
        BASIN['gis_data']['population'],
        BASIN['gis_data']['basin_mask'],
        BASIN['data_cube']['monthly']['lu'],
        wcpc=110,
        flow_type='supply')
residential_demand_nc=hl.calc_residential_water_consumption(
        BASIN['gis_data']['population'],
        BASIN['gis_data']['basin_mask'],
        BASIN['data_cube']['monthly']['lu'],
        wcpc=100,
        flow_type='demand')

### split return flow by source sw/gw
return_sw_from_sw_nc,return_sw_from_gw_nc=hl.split_flow(
        sw_return_nc,fraction_nc=sw_supply_fraction_nc)
return_gw_from_sw_nc,return_gw_from_gw_nc=hl.split_flow(
        gw_return_nc,fraction_nc=sw_supply_fraction_nc)

### split residential supply by sw/gw fraction
f=BASIN['params']['residential_sw_supply_fraction']
sw_residential_supply_nc,gw_residential_supply_nc=hl.split_flow(
        residential_supply_nc,fraction=f)
### add residential sw/gw supply to sw/gw supply and sw/gw return
BASIN['data_cube']['monthly']['supply_sw']=hl.add_flow(
        sw_supply_nc,sw_residential_supply_nc,name='total_sw_supply')
BASIN['data_cube']['monthly']['supply_gw']=hl.add_flow(
        gw_supply_nc,gw_residential_supply_nc,name='total_gw_supply')

#assume that residential supply from sw return to sw, from gw return to gw
BASIN['data_cube']['monthly']['return_sw_from_sw']=hl.add_flow(
        return_sw_from_sw_nc,sw_residential_supply_nc,name='total_return_sw_from_sw')
BASIN['data_cube']['monthly']['return_gw_from_gw']=hl.add_flow(
        return_gw_from_gw_nc,gw_residential_supply_nc,name='total_return_gw_from_gw')

BASIN['data_cube']['monthly']['return_sw_from_gw']=return_sw_from_gw_nc
BASIN['data_cube']['monthly']['return_gw_from_sw']=return_gw_from_sw_nc
### add residential demand to total demand
BASIN['data_cube']['monthly']['demand']=hl.add_flow(
        demand_nc,residential_demand_nc,name='total_demand')
### total return and supply
BASIN['data_cube']['monthly']['return_sw']=hl.add_flow(
        BASIN['data_cube']['monthly']['return_sw_from_gw'],
        BASIN['data_cube']['monthly']['return_sw_from_sw'],
        name='return_sw')
BASIN['data_cube']['monthly']['return_gw']=hl.add_flow(
        BASIN['data_cube']['monthly']['return_gw_from_gw'],
        BASIN['data_cube']['monthly']['return_gw_from_sw'],
        name='return_gw')
BASIN['data_cube']['monthly']['supply']=hl.add_flow(
        BASIN['data_cube']['monthly']['supply_sw'],
        BASIN['data_cube']['monthly']['supply_gw'],
        name='total_supply')
### calculate recharge
BASIN['data_cube']['monthly']['recharge']=BASIN['data_cube']['monthly']['perc']
#not sure if perc from SM_balance include percincr or not
#BASIN['data_cube']['monthly']['recharge']=hl.substract_flow(
#        BASIN['data_cube']['monthly']['perc'],
#        BASIN['data_cube']['monthly']['percincr'],
#        name='recharge')

#%% Calculate monthly discharge timeseries
### Calculate subbasin-wide timeseries

for sb in BASIN['gis_data']['subbasin_mask']:
    subbasin={}
    for key in ['sro','return_sw','bf','supply_sw']:
        output=os.path.join(BASIN['output_folder'],
                            'subbasin_{0}_{1}.csv'.format(sb,key))
        df=cf.calc_flux_per_basin(BASIN['data_cube']['monthly'][key],
                               BASIN['gis_data']['subbasin_mask'][sb],
                               output=output)
        subbasin[key]=df        
    ## read subbasin inflow
    if len(BASIN['params']['dico_in'][sb])==0: #no inflow
        inflow=None 
    else: #1 or more inflows
        for i in range(len(BASIN['params']['dico_in'][sb])):            
            if BASIN['params']['dico_in'][sb][i] == 0: #inflow from outside
                df_inflow=pd.read_csv(BASIN['ts_data']['q_in_sw'][sb],
                                      sep=';',index_col=0)
            else: #inflow from upstream subbasin                  
                subbasin_in=BASIN['params']['dico_in'][sb][i]
                df_inflow=pd.read_csv(BASIN['ts_data']['q_outflow'][subbasin_in],
                                      sep=';',index_col=0) 
                df_inflow=df_inflow*BASIN['unit_conversion']
                #assuming that outflow of upstream subbasin was calculated before
            if i == 0:
                inflow=df_inflow
            else:
                inflow=inflow+df_inflow

    ## calculate sw discharge and dS from pixel-based model results
    output=os.path.join(BASIN['output_folder'],
                        'subbasin_{0}_{1}.csv'.format(sb,'{0}'))
    discharge,dS_sw=hl.calc_sw_from_wp(subbasin['sro'],
                                       subbasin['return_sw'],
                                       subbasin['bf'],
                                       subbasin['supply_sw'],
                                       inflow=inflow,
                                       output=output,
                                       outflow=True, #not endorheic basin
                                       unit_conversion=BASIN['unit_conversion'] #MCM
                                       )
    BASIN['ts_data']['q_outflow'][sb]=discharge
    BASIN['ts_data']['dS_sw'][sb]=dS_sw
    
# outflow of basin is outflow of downstream subbasin   
for sb in BASIN['params']['dico_out']:
    if 0 in BASIN['params']['dico_out'][sb]: #if subbasin outflow is basin outflow
        BASIN['ts_data']['q_outflow']['basin']=BASIN['ts_data']['q_outflow'][sb]

# dS_sw of basin is sum of dS_sw of all subbasins
for i in range(len(BASIN['gis_data']['subbasin_mask'])):
    sb=list(BASIN['gis_data']['subbasin_mask'].keys())[i]
    df=pd.read_csv(BASIN['ts_data']['dS_sw'][sb],sep=';',index_col=0)
    if i==0:
        dS_sw=df
    else:
        dS_sw=dS_sw+df
dS_sw.to_csv(os.path.join(BASIN['output_folder'],
                        'basin_dS_sw.csv'),sep=';')
BASIN['ts_data']['dS_sw']['basin']=os.path.join(BASIN['output_folder'],
                        'basin_dS_sw.csv')
#%% yearly datacube
for key in BASIN['data_cube']['monthly']:
    if key != 'lu':
        BASIN['data_cube']['yearly'][key]=cf.create_yearly_dataset(
                BASIN['data_cube']['monthly'][key], hydroyear=BASIN['hydroyear'])
#%% Calculate yearly and (intermediate) monthly sheet csvs
from WAsheets import sheet1
from WAsheets import sheet2
from WAsheets import sheet3
from WAsheets import sheet4
from WAsheets import sheet5
from WAsheets import sheet6

sheet1_yearly_csvs=sheet1.main(BASIN,unit_conversion=BASIN['unit_conversion'])
sheet2_yearly_csvs=sheet2.main(BASIN,unit_conversion=BASIN['unit_conversion'])
sheet3_yearly_csvs=sheet3.main(BASIN,unit_conversion=BASIN['unit_conversion'])
sheet4_yearly_csvs=sheet4.main(BASIN,unit_conversion=BASIN['unit_conversion'])
sheet5_yearly_csvs=sheet5.main(BASIN,unit_conversion=BASIN['unit_conversion'])
sheet6_yearly_csvs=sheet6.main(BASIN,unit_conversion=BASIN['unit_conversion'])
#%% Print hydro-yearly sheet csv
from WAsheets import print_sheet as ps

if BASIN['unit_conversion'] == 1e6:
    str_unit='km3/year'
elif BASIN['unit_conversion']==1e3:
    str_unit='MCM/year'
    
for sheet1_csv in sheet1_yearly_csvs:
    period=os.path.basename(sheet1_csv).split('.')[0].split('_')[-1]
    output=sheet1_csv.replace('.csv','.pdf')
    ps.print_sheet1(BASIN['name'],period=period,
                    output=output,units=str_unit,data=sheet1_csv)

for sheet2_csv in sheet2_yearly_csvs:
    period=os.path.basename(sheet2_csv).split('.')[0].split('_')[-1]
    output=sheet2_csv.replace('.csv','.pdf')
    ps.print_sheet2(BASIN['name'],period=period,output=output,
                    units=str_unit,data=sheet2_csv)
    
for sheet3_csv in sheet3_yearly_csvs:
    period=os.path.basename(sheet3_csv).split('.')[0].split('_')[-1]
    output=sheet3_csv.replace('.csv','.pdf')
    ps.print_sheet3(BASIN['name'],period=period,output=output,units=[str_unit,'kg/ha','kg/m3'],data=sheet3_csv)
    
for sheet4_csv in sheet4_yearly_csvs:
    period=os.path.basename(sheet4_csv).split('.')[0].split('_')[-1]
    output=[sheet4_csv.replace('.csv','_part1.pdf'),sheet4_csv.replace('.csv','_part2.pdf')]
    ps.print_sheet4(BASIN['name'],period=period,output=output,
                    units=[str_unit,str_unit],data=[sheet4_csv,sheet4_csv])
    
for sheet5_csv in sheet5_yearly_csvs:
    period=os.path.basename(sheet5_csv).split('.')[0].split('_')[-1]
    output=sheet5_csv.replace('.csv','.pdf')
    ps.print_sheet5(BASIN['name'],sb_codes=[1,2],
                    period=period,output=output,units=str_unit,
                    data=sheet5_csv)
    
for sheet6_csv in sheet6_yearly_csvs:
    period=os.path.basename(sheet6_csv).split('.')[0].split('_')[-1]
    output=sheet6_csv.replace('.csv','.pdf')
    ps.print_sheet6(BASIN['name'],period=period,output=output,units=str_unit,
                    data=sheet6_csv)

#%% Pickle BASIN
import pickle
import datetime

time=str(datetime.datetime.now()).replace(' ','').replace(':','-')
filename=os.path.join(BASIN['output_folder'],
                      'BASIN_{0}.pickle'.format(time))
with open(filename, 'wb') as handle:
    pickle.dump(BASIN, handle)
    
#%% Test discharge

#bf_nc=r"D:\Projects\ADB\analysis\Karnataka\Results_Karnataka\New_SMBalance\K4\bf_monthly.nc"
#sroincr_nc=r"D:\Projects\ADB\analysis\Karnataka\Results_Karnataka\New_SMBalance\K4\d_sro_monthly.nc"
#sro_nc=r"D:\Projects\ADB\analysis\Karnataka\Results_Karnataka\New_SMBalance\K4\sro_monthly.nc"
#supply_nc=r"D:\Projects\ADB\analysis\Karnataka\Results_Karnataka\New_SMBalance\K4\supply_monthly.nc"
#
#basin_mask=r"D:\Projects\ADB\analysis\Karnataka\Results_Karnataka\New_SMBalance\K4\Shape\K4.tif"
#
#df_bf=cf.calc_flux_per_basin(bf_nc, basin_mask,
#                          chunksize=None,output=None,quantity='volume')
#df_sro=cf.calc_flux_per_basin(sro_nc, basin_mask,
#                          chunksize=None,output=None,quantity='volume')
#df_sroincr=cf.calc_flux_per_basin(sroincr_nc, basin_mask,
#                          chunksize=None,output=None,quantity='volume')
#df_supply=cf.calc_flux_per_basin(supply_nc, basin_mask,
#                          chunksize=None,output=None,quantity='volume')
#df_supply=df_supply*0.5
#discharge,dS_sw=hl.calc_sw_from_wp(df_sro,df_sroincr,df_bf,df_supply,plot=True)
#
#plt.plot(discharge,'k',label='discharge')
#plt.plot(dS_sw,'orange',label='dS_sw')
#plt.legend()

