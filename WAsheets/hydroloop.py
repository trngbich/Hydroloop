# -*- coding: utf-8 -*-
"""
Created on Fri May 1 15:09:33 2020

@author: ntr002

Functions to calculate extra data needed for WA+ sheets from pixel-based model

Input:
    
Output:

"""
import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import calendar
from itertools import groupby
from operator import itemgetter
from . import calculate_flux as cf
from . import get_dictionaries as gd
from . import GIS_functions as gis


#%% General
def check_requirement_sheet(BASIN,requirements,
                             gis_data=True,
                             data_cube=True,
                             param=True): 
    #check parameters requirements
    if len(requirements['param']-BASIN['gis_data'].keys())>0:
        print('Required parameters for Sheet {0} are missing: '.format(
                requirements['sheet']))
        print(requirements['gis_data']-BASIN['gis_data'].keys())
        gis_data=False        
    #check GIS data requirements
    if len(requirements['gis_data']-BASIN['gis_data'].keys())>0:
        print('Required GIS data for Sheet {0} are missing: '.format(
                requirements['sheet']))
        print(requirements['gis_data']-BASIN['gis_data'].keys())
        gis_data=False
    #check data cube requirements
    if len(requirements['data_cube']-BASIN['data_cube'].keys())>0:
        missing_data=list(requirements['data_cube']-BASIN['data_cube'].keys())
        print('Required spatial data for Sheet {0} are missing: '.format(
                requirements['sheet']))
        print(missing_data)
        for data in missing_data:            
            if data.replace('yearly','monthly') in BASIN['data_cube'].keys():
                print('Convert {0} from {1}'.format(data,
                      data.replace('yearly','monthly')))
                output=cf.create_yearly_dataset(
                        BASIN['data_cube'][data.replace('yearly',
                                         'monthly')],
                                        hydroyear=BASIN['hydroyear'],
                                        chunksize=BASIN['chunksize']
                                         )
                BASIN['data_cube'][data]=output
            else:
                print('Missing mothly data to convert to {0}'.format(data))
                data_cube=False
    
    return gis_data*data_cube*param
#%% Sheet 2 functions
def split_ETI(et_nc,lai_nc,output_folder):
    '''
    split actual evapotranspiration into 3 components
    Evaporation, Transpiration, and Interception
    '''
    return

#%% Sheet 3 functions
def calc_GBP(npp_nc,output_folder):
    '''
    calculate Gross Biomass Production
    '''
    return

def calc_GBWP(gbp_nc, et_nc, output_folder):
    '''
    calculate Gross Biomass Water Productivity
    '''
    return

def calc_crop_val_per_season(val_nc,lu_map,crop_lucl,
                           begin='2009-01-01',end='2009-06-12'):
    '''
    val_nc
    '''
    return
#%% Sheet 4 functions
def calc_sw_supply_fraction_by_LU(lu_nc,aeisw_nc,
                                  output=None,
                                  chunksize=None):
    '''
    calculate SW supply fraction by LU map
    lu_nc: string
        path to Land Use datacube netCDF file
    aeisw_nc: string
        Area equipped with surface water irrigation map (netCDF)
        
    '''
    lucs = gd.get_sheet4_6_classes() 
    fractions = gd.get_sheet4_6_fractions() #SW supply fractions
    LU=cf.open_nc(lu_nc,chunksize=chunksize,layer=0)   
    
    sw_supply_fraction=LU*0
    # fraction from dictionary
    for key in list(fractions.keys()):
        classes = lucs[key]
        fraction = fractions[key]
        print('{0} {1} {2}'.format(key,classes,fraction))
        sw_supply_fraction = xr.where(
                LU.isin(classes),fraction,sw_supply_fraction)        
    # update fraction with aeisw (GMIA)
    aeisw=cf.open_nc(aeisw_nc,chunksize=chunksize,layer=0)
    aeisw=aeisw/100 #convert percentage to fraction
    aeisw = xr.where(xr.isnan(aeisw),
                     aeisw.mean(skipna=True),
                     aeisw)  #assume average aeisw where NaN
    sw_supply_fraction = xr.where(
                LU.isin(lucs['Irrigated crops']),# for irrigated crops
                aeisw,#sw_supply_fraction = aeisw
                sw_supply_fraction)          
    #save output    
    sw_supply_fraction=sw_supply_fraction.transpose(
            'time','latitude','longitude')               

    sw_supply_fraction.attrs={'units':'fraction',
                              'source':'LUWA, GMIA_aeisw',
                              'quantity':'Fraction of SW supply'
                                  }
    sw_supply_fraction.name = 'sw_supply_fraction'
    
    comp = dict(zlib=True, 
                complevel=9, 
                least_significant_digit=2, 
                chunksizes=chunksize)
    encoding = {'sw_supply_fraction': comp}    
    if output is None:
        output=os.path.join(os.path.dirname(lu_nc),
                            'sw_supply_fraction.nc')
    sw_supply_fraction.load().to_netcdf(output,encoding=encoding)
    print('Save monthly LU datacube as {0}'.format(output))
    del LU
    return sw_supply_fraction

def calc_sw_return_fraction():
    '''
    DTOT = DSRO + DPERC
    SWRETFRAC = LU * 0
    SWRETFRAC[DTOT > 0] = (DSRO/(DTOT))[DTOT > 0]
    '''
    return
def split_flow_by_fraction(flow_nc,fraction_nc,
               output=None,chunksize=None,
               sub_names = ['sw','gw'] ):
    '''
    split flow_nc into 2 flows using a fraction map
    for example, used for splitting supply into sw and gw supply, 
    or for splitting return flow into sw and gw return flow
    '''
    flow=cf.open_nc(flow_nc,chunksize=chunksize,layer=0)
    fraction=cf.open_nc(fraction_nc,chunksize=chunksize,layer=0)
    
    sw_flow=flow*fraction
    gw_flow=flow-sw_flow
    
    #modify attributes of splitted flow datasets
    sw_flow=sw_flow.transpose('time','latitude','longitude')
    gw_flow=gw_flow.transpose('time','latitude','longitude')
    name=flow.name
    sw_flow.name='{0}_{1}'.format(name,sub_names[0])
    gw_flow.name='{0}_{1}'.format(name,sub_names[1])
    sw_flow.attrs={'units':flow.attrs['units'],
                              'source':'LUWA, GMIA_aeisw',
                              'quantity':flow.attrs['quantity']
                                  }
    gw_flow.attrs={'units':flow.attrs['units'],
                              'source':'LUWA, GMIA_aeisw',
                              'quantity':flow.attrs['quantity']
                                  }
    #save output
    if output is None:
        output=flow_nc.replace('.nc','_{0}.nc')
    comp = dict(zlib=True, 
                complevel=9, 
                least_significant_digit=2, 
                chunksizes=chunksize)
    sw_output=output.format(sub_names[0])
    gw_output=output.format(sub_names[1])
    print('Save splitted {0} as {1} and {2}'.format(name,sw_output,gw_output))
    sw_flow.load().to_netcdf(sw_output,encoding={sw_flow.name:comp})
    gw_flow.load().to_netcdf(gw_output,encoding={gw_flow.name:comp})    
    return sw_flow,gw_flow

def calc_nonconsumed_supply(supply_nc,etincr_nc,
                            output=None, chunksize=None):
    '''
    calculate non-consumed supply by substracting Incremental ET from supply
    
    supply_nc: str
        path to Supply netCDF file
    etincr_nc: str
        path to Incremental ET netCDF file
    '''
    supply=cf.open_nc(supply_nc,chunksize=chunksize,layer=0)
    etincr=cf.open_nc(etincr_nc,chunksize=chunksize,layer=0)
    
    non_consumed=supply-etincr
    
    non_consumed.name='non_consumed_supply'
    non_consumed=non_consumed.transpose('time','latitude','longitude')
    non_consumed.attrs={'units':supply.attrs['units'],
                              'source':'Supply - ET Incremental',
                              'quantity':'non_consumed_supply'
                                  }
    if output is None:
        output=os.path.join(os.path.dirname(supply_nc),
                            'non_consumed_supply.nc')
    comp = dict(zlib=True, 
                complevel=9, 
                least_significant_digit=2, 
                chunksizes=chunksize)
    non_consumed.load().to_netcdf(output,
                     encoding={'non_consumed_supply':comp})       
    return non_consumed

def calc_land_surface_water_demand(lai_nc, etref_nc, p_nc, lu_nc, 
                output=None,chunksize=None):
    '''
    calculate water demand of land surface based on LAI
    lai_nc: str
        path to LAI netcdf file
    etref_nc: str
        path to ET reference netcdf file
    p_nc: str
        path to Precipiation netcdf file
    lu_nc: str
        path to LU netcdf file (invariant or monthly)
        
    return 
    demand: xr.DataArray
        water demand datacube
    '''
    #read netcdf data
    lai=cf.open_nc(lai_nc,chunksize=chunksize,layer=0)
    et_reference=cf.open_nc(etref_nc,chunksize=chunksize,layer=0)
    lu=cf.open_nc(lu_nc,chunksize=chunksize,layer=0)
    p=cf.open_nc(p_nc,chunksize=chunksize,layer=0)
    
    #calculate Potential ET using LAI-based crop coefficient Kc
    kc = xr.where(lu.isin([4, 5, 30, 23, 24, 63]), #mask water classes
                  1.4, #water KC
                  (1 - xr.ufuncs.exp(-0.5 * lai)) / 0.76) #non-water KC   
    
    et_potential = kc*et_reference
    #calculate land surface water demand as the gap between potential ET and effective rainfall
    phi = et_potential / p
    
    effective_rainfall = xr.ufuncs.sqrt(phi*xr.ufuncs.tanh(1/phi)\
                                        *(1-xr.ufuncs.exp(-phi))) * p
    
    demand = et_potential - effective_rainfall
    #add attributes to demand data
    demand.name='land_surface_water_demand'
    demand=demand.transpose('time','latitude','longitude')
    demand.attrs={'units':p.attrs['units'],
                              'source':'Potential ET - Effective Rainfall',
                              'quantity':'land_surface_water_demand'
                                  }   
    #save results
    if output is None:
        output=os.path.join(os.path.dirname(lai_nc),'water_demand.nc')
    comp = dict(zlib=True, 
                complevel=9, 
                least_significant_digit=2, 
                chunksizes=chunksize)
    demand.load().to_netcdf(output,
                     encoding={'land_surface_water_demand':comp})           
    return demand

def add_residential_water_consumption(population_tif,area_tif,
                               lu_nc,blue_water_nc,
                               wcpc=100,flow_type='demand',
                               chunksize=None,
                               output=None):
    '''
    calculate domestic water demand or supply based on population map and
    water consumption per capita and add values to demand and supply dataset.
    
    wcpc: number
        water consumption per capita (L/person/day)
        use minimal wcpc to calculate residential demand 
        and actual wcpc to calculate residential supply
    lu_nc: str
        path to monthly LU NetCDF file
    blue_water_nc: str
        path to demand or supply NetCDF file
    flow_type: str
        'demand' or 'supply'
    '''
    # Read input data
    population=gis.OpenAsArray(population_tif,nan_values=True)
    area=gis.OpenAsArray(area_tif,nan_values=True)
    lu=cf.open_nc(lu_nc,chunksize=chunksize,layer=0)
    blue_water=cf.open_nc(blue_water_nc,chunksize=chunksize,layer=0)
    #get residential classes
    sheet4_lucs=gd.get_sheet4_6_classes() 
    classes = sheet4_lucs['Residential']    
    # Convert [pop/ha] to [pop/pixel]
    population *= area * 100    
    # Calculate WC per pixel in [mm/day]
    residential_wc = wcpc * population * 10**-6 / area
    
    # Calculate number days per month
    months=pd.to_datetime(np.array(blue_water.time))
    ndays=np.array(
            [calendar.monthrange(date.year, 
                                 date.month)[1] for date in months])
    ndays_da = xr.DataArray(ndays, dims=["time"])
    ndays_da['time']=blue_water['time']
    
    # Calculate WC per pixel in [mm/month]
    mask=xr.where(lu.isin(classes),1,0)
    ndays_mask=mask*ndays_da
    residential_wc*=ndays_mask
    # Add residential water to total blue water
    mean=residential_wc.mean(dim=['latitude','longitude'],skipna=True)
    residential_wc=xr.where(xr.ufuncs.isnan(residential_wc),
                            -999,#fill nan temporarily with finite value
                            residential_wc)
    nan_residential_wc=xr.where( #residential pixels that has nan wc
        lu.isin(classes) * (residential_wc==-999),
                    1.0,#True
                    0.0) #False
    residential_wc+=nan_residential_wc*mean #fill nan residential_wc w/ mean 
    total_blue_water=blue_water+residential_wc #add residential_wc
   
    # Save total_blue_water
    total_blue_water=total_blue_water.transpose('time',
                                                 'latitude','longitude')
    total_blue_water.name=blue_water.name
    total_blue_water.attrs={'units':blue_water.attrs['units'],
                              'source':blue_water.attrs['source']\
                        + ' , added residential water {0}'.format(flow_type),
                              'quantity':blue_water.attrs['quantity']
                                  }
    if output is None:
        output=blue_water_nc.replace('.nc',
                                     '_with_residential_{0}.nc'.format(
                                             flow_type))
    comp = dict(zlib=True, 
                complevel=9, 
                least_significant_digit=2, 
                chunksizes=chunksize)
    print('Save {0} added residential water as {1}'.format(
            blue_water.name,output))
    total_blue_water.load().to_netcdf(output,
                         encoding={total_blue_water.name:comp})    
    #close files and return results
    blue_water.close()
    lu.close()
    del population
    del area
    return output

#%% Sheet 5 functions
def calc_sw_from_wp(sro,sroincr,bf,supply_sw,
                         inflow=None,outflow=True,plot=False):
    '''
    calculate surface water outflow and storage from pixel-based model results
    In DataFrame format
    TF: Total flow = SRO + SROincr + BF
    SRO: Surface runoff
    SROincr: Incremental surface runoff
    BF: Baseflow
    SUPsw: Surface water Supply
    
    
    '''
    total_flow=sro[sro.columns[0]]\
                +sroincr[sroincr.columns[0]]\
                +bf[bf.columns[0]]
    available_sw=total_flow-supply_sw[supply_sw.columns[0]]
    if inflow is not None:
        available_sw=available_sw+inflow[inflow.columns[0]]
    #plot flows
    if plot:
        plt.plot(sro,'g',label='SRO')
        plt.plot(sroincr,'b',label='SROincr')
        plt.plot(bf,'r',label='BF')
        plt.plot(supply_sw,'lightblue',label='SupplySW')
        plt.legend()
    discharge,dS_sw=distribute_available_sw(available_sw,outflow=outflow)
    index=available_sw.index
    discharge=pd.DataFrame(data=discharge,index=index,columns=['Discharge'])
    dS_sw=pd.DataFrame(data=dS_sw,index=index,columns=['dS_sw'])
    return discharge,dS_sw

def distribute_available_sw(available_sw,outflow=True):
    '''
    Distribute monthly available surface water into discharge and surface
    water storage change.
    available_sw= total_runoff + other_inflow - interbasin_outflow - sw_supply
    total_runoff= surface_runoff + incremental_surface_runoff + base_flow
    
    Assumption: The total negative available surface water in dry period is 
    compensated by the storing water to SW storage from the previous month
    with positive available surface water. The amount of stored water each 
    month is proportional to the positive available water each during the wet 
    period.
    ex.
    available_sw=[8 4 0 -2 -2 -2]
    return
    dS=[4 2 0 -2 -2 -2]
    discharge=[4 2 0 0 0 0]
    ----
    available_sw:np.array or pandas.Series
    outflow: True (default)
        True for basin with outflow
        False for endorheic basin
    available_sw=np.array([8,4,0,-2,-2,-2])
    '''
  
    ### if basin/subbasin has surface water outflow
    if outflow is False:
        dS=available_sw
        discharge = np.zeros(len(available_sw))
        
    ### if basin/subbasin has no surface water outflow
    else:
        dS = np.zeros(len(available_sw))
        for i in np.where(available_sw < 0):
            dS[i] = available_sw[i] 
            # DeltaS = RO(includes_returns) + otherInflows - W  where negative
    
        discharge = np.copy(available_sw )
        discharge[discharge < 0] = 0
        # Discharge = RO(includes_returns) + otherInflows - W where positive 
        
        ##backwards estimate SW storage increase in previous months based on negative dS        
        dS_temp=np.copy(dS)
        data = list(np.where(dS_temp < 0)[0])
        prev_end = 0
        for k, g in groupby(enumerate(data), lambda i_x: i_x[0]-i_x[1]):
            decr_st = list(map(itemgetter(1), g))
            start = decr_st[0]
            end = decr_st[-1]
            deltaS_season = np.sum(dS_temp[start:end+1])
            avail_per_month = available_sw[prev_end:start]
            avail_prev_season = np.sum(available_sw[prev_end:start])
            weight_per_month = avail_per_month/avail_prev_season
            if start > 0:
                if avail_prev_season < abs(deltaS_season):
                    print('Warning, insufficient runoff and inflows in previous months to meet DeltaS')
                    dS_temp[prev_end:start] = np.min((avail_per_month, -deltaS_season * weight_per_month), axis=0)
                    discharge[prev_end:start] = discharge[prev_end:start] \
                                             - dS_temp[prev_end:start]
                else:
                    dS_temp[prev_end:start] = dS_temp[prev_end:start] \
                                - deltaS_season * weight_per_month
                    discharge[prev_end:start] = discharge[prev_end:start]\
                                                    - dS_temp[prev_end:start]
            prev_end = end + 1
        dS = dS_temp        
    return discharge,dS

#%% Sheet 6 functions

def calc_recharge():
    '''
    recharge = percolation - incremental percolation
    '''
    return