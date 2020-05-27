# -*- coding: utf-8 -*-
"""
Created on Fri May 1 15:09:33 2020

@author: ntr002

Functions to calculate extra data needed for WA+ sheets from pixel-based model results (NetCDF dataset)

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
                             params=True): 
    #check parameters requirements
    if len(requirements['params']-BASIN['gis_data'].keys())>0:
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
    
    return gis_data*data_cube*params
#%% Sheet 2 functions
    
def split_ETI(et_nc,i_nc=None,t_nc=None,
              chunksize=None,**input_data):
    '''
    split actual evapotranspiration into 3 components
    Evaporation, Transpiration, and Interception
    
    '''
    keys=input_data.keys()
    #Interception
    if i_nc is None:
        #check input data        
        list_rq=['p_nc','lai_nc','nrd_nc']
        check=[item in keys for item in list_rq]
        if not all(check): #not enough input_data
            print('ERROR:Missing data to calculate interception')
            print('p_nc: path to precipitation dataset')
            print('lai_nc: path to LAI dataset')
            print('nrd_nc: path to number of rainy days dataset')            
            return None
        else: #input_data is sufficient to calculate i
            p_nc=input_data['p_nc']
            lai_nc=input_data['lai_nc']
            nrd_nc=input_data['nrd_nc']
            i_nc=calc_interception(p_nc,lai_nc,nrd_nc,
                              chunksize=chunksize)
            input_data['i_nc']=i_nc
            return
    #Transpiration
    if t_nc is None:
        #check input data        
        list_rq=['et_nc','i_nc','ndm_nc']
        check=[item in keys for item in list_rq]
        if not all(check): #not enough input_data
            print('ERROR:Missing data to calculate             transpiration')
            print('et_nc: path to actual evapotranspiration                  dataset')
            print('ndm_nc: path to Net-Dry-Matter dataset')
            print('i_nc: path to interception dataset')
            return None
        else: #input_data is sufficient to calculate t        
            if 'ndm_max_original' in keys:#check ndm_max method
                ndm_max_original=input_data['ndm_max_original']
            else: #if no method specified, keep False
                ndm_max_original=False
            et_nc=input_data['et_nc'] #
            i_nc=input_data['i_nc']
            ndm_nc=input_data['ndm_nc']
            t_nc=calc_transpiration(et_nc,i_nc,ndm_nc,
                       chunksize=chunksize,
                       ndm_max_original=ndm_max_original)            
    #Evaporation
    e_nc=calc_evaporation(et_nc,i_nc,t_nc,chunksize=chunksize)            
    return e_nc,i_nc,t_nc

def calc_interception(p_nc,lai_nc,nrd_nc,
                      output=None,
                      chunksize=None):
    '''
    calculate interception from LAI and precipitation
    '''
    #open dataset
    p=cf.open_nc(p_nc,chunksize=chunksize,layer=0)
    lai=cf.open_nc(lai_nc,chunksize=chunksize,layer=0)
    nrd=cf.open_nc(nrd_nc,chunksize=chunksize,layer=0)
    #fill nan
    
    #calculate Interception
    i = lai * (1 - (1 + (p/nrd) \
                    * (1 - xr.ufuncs.exp(-0.5 * lai)) \
                    * (1/lai))**-1) * nrd
    #fill in nan value due to missing lai value    
    i = xr.where(xr.ufuncs.isnan(i),0,i)+p*0
    #save dataset
    i=i.transpose('time','latitude','longitude') 
    i.name='Interception'
    i.attrs={'units':p.attrs['units'],
    'source':'LAI* (1-(1+(P/nRD)*(1-exp(-0.5*LAI))*(1/LAI))**-1)*nRD',
    'quantity':'interception'
    }
    comp = dict(zlib=True, 
                complevel=9, 
                least_significant_digit=2, 
                chunksizes=chunksize)
    encoding = {i.name: comp}    
    if output is None:
        output=os.path.join(os.path.dirname(p_nc),
                            'interception.nc')
    print('Save Interception datacube as {0}'.format(output)) 
    i.load().to_netcdf(output,encoding=encoding)
    #close dataset
    p.close()
    lai.close()
    nrd.close()
    return output

def calc_transpiration(et_nc,i_nc,ndm_nc,
                       chunksize=None,
                       output=None,
                       ndm_max_original=False):
    '''
    calculate transpiration from net-dry-matter, given actual 
    evapotranspiration and interception
    '''
    #open dataset
    et=cf.open_nc(et_nc,chunksize=chunksize,layer=0)
    i=cf.open_nc(i_nc,chunksize=chunksize,layer=0)
    ndm=cf.open_nc(ndm_nc,chunksize=chunksize,layer=0)
    #calculate maximum NDM
    if not ndm_max_original:
        ndm_month_avg=ndm.groupby('time.month').mean(
                dim=['time'],skipna=True) #calculate monthly mean
        ndm_month_std=ndm.groupby('time.month').std(
                dim=['time'],skipna=True) #calculate monthly std
        ndm_month_max=ndm_month_avg+2*ndm_month_std
        #spread monthly max to all time steps in datacube
        ndm_max=ndm*0        
        for i in range(len(ndm_max)):
           m=pd.to_datetime(ndm.time[0].data).month #month index
           ndm_max[i]=ndm_month_max[m-1]
        
    if ndm_max_original:
        ndm_month_max=ndm.groupby('time.month').max(
                dim=['time'],skipna=True) #calculate monthly max      
        #spread monthly max to all time steps in datacube
        ndm_max=ndm*0  
        for i in range(len(ndm_max)):
           m=pd.to_datetime(ndm.time[0].data).month #month index
           ndm_max[i]=ndm_month_max[m-1]/0.95
    #calculate Transpiration
    t=xr.ufuncs.minimum((ndm/ndm_max),0.95)*(et-i)
    #fill in nan value due to missing ndm value    
    t = xr.where(xr.ufuncs.isnan(t),0,t)+et*0
    #save dataset
    t=t.transpose('time','latitude','longitude') 
    t.name='transpiration'
    t.attrs={'units':et.attrs['units'],
    'source':'minimum(0.95,NDM/NDM_max)*(ETa-I)',
    'quantity':'Transpiration'
    }
    comp = dict(zlib=True, 
                complevel=9, 
                least_significant_digit=2, 
                chunksizes=chunksize)
    encoding = {t.name: comp}    
    if output is None:
        output=os.path.join(os.path.dirname(et_nc),
                            'transpiration.nc')
    print('Save Transpiration datacube as {0}'.format(output)) 
    t.load().to_netcdf(output,encoding=encoding)
    #close dataset
    et.close()
    i.close()
    ndm.close()
    return output

def calc_evaporation(et_nc,i_nc,t_nc,
                     chunksize=None,
                     output=None):
    '''
    calculate evaporation from actual evapotranspiration, given interception
    and transpiration portions
    
    '''
    #open dataset
    et=cf.open_nc(et_nc,chunksize=chunksize,layer=0)
    i=cf.open_nc(i_nc,chunksize=chunksize,layer=0)
    t=cf.open_nc(t_nc,chunksize=chunksize,layer=0)
    #calculate evaporation
    e=et-i-t
    #save dataset
    e=e.transpose('time','latitude','longitude') 
    e.name='evaporation'
    e.attrs={'units':et.attrs['units'],
    'source':'ETa - I - T',
    'quantity':'Evaporation'
    }
    comp = dict(zlib=True, 
                complevel=9, 
                least_significant_digit=2, 
                chunksizes=chunksize)
    encoding = {e.name: comp}    
    if output is None:
        output=os.path.join(os.path.dirname(et_nc),
                            'evaporation.nc')
    print('Save Evaporation datacube as {0}'.format(output)) 
    e.load().to_netcdf(output,encoding=encoding)
    #close dataset
    et.close()
    i.close()
    t.close()    
    
    return output

def calc_ndm(npp_nc):
    '''
    Calculate Net Dry Matter based on monthly Net Primary Production
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
def add_flow(flow_nc,additional_flow_nc,
             flow_name='total_flow',
             output=None,chunksize=None):
    '''
    Add additional_flow_nc to flow_nc and save to new dataset
    '''
    flow=cf.open_nc(flow_nc,chunksize=chunksize,layer=0)
    additional_flow=cf.open_nc(additional_flow_nc,
                               chunksize=chunksize,layer=0)
    source_info='sum of '+flow.name+' and '+additional_flow.name
    total_flow=flow+additional_flow
    #Add attributes
    total_flow.name=flow_name
    total_flow=total_flow.transpose('time',
                                    'latitude','longitude')
    total_flow.attrs={'units':flow.attrs['units'],
                    'source':source_info,
                    'quantity':flow_name
                                  }    
    if output is None:
        output=os.path.join(
                os.path.dirname(flow_nc),
                '{0}.nc'.format(flow_name))
    comp = dict(zlib=True, 
                complevel=9, 
                least_significant_digit=2, 
                chunksizes=chunksize)
    print('Save summed {0} as {1}'.format(flow_name,output))
    total_flow.load().to_netcdf(output,
                 encoding={total_flow.name:comp})  
    return output  

def split_flow(flow_nc,
               fraction=0.5,
               fraction_nc=None,                           
               output=None,
               chunksize=None,
               sub_names = ['sw','gw'] ):
    '''
    split flow_nc into 2 flows using a fraction map or value.
    For example, used to split supply into sw and gw supply, 
    or for split return flow into sw and gw return flow
    '''
    flow=cf.open_nc(flow_nc,chunksize=chunksize,layer=0)
    source_info=flow.name
    if fraction_nc is not None: #use fraction map
        fraction=cf.open_nc(fraction_nc,
                            chunksize=chunksize,layer=0)
        source_info+= ' multiplied by ' +fraction.name
    else:
        source_info+= ' multiplied by {0}'.format(fraction)

    flow_one=flow*fraction 
    flow_two=flow-flow_one 
    
    #modify attributes of splitted flow datasets
    flow_one=flow_one.transpose('time','latitude','longitude')
    flow_two=flow_two.transpose('time','latitude','longitude')
    name=flow.name
    flow_one.name='{0}_{1}'.format(name,sub_names[0])
    flow_two.name='{0}_{1}'.format(name,sub_names[1])
    flow_one.attrs={'units':flow.attrs['units'],
                    'source':source_info,
                    'quantity':flow.attrs['quantity']
                                  }
    flow_two.attrs={'units':flow.attrs['units'],
                    'source':flow.name+' minus ' +source_info,
                    'quantity':flow.attrs['quantity']
                                  }
    if output is None:
        output=flow_nc.replace('.nc','_{0}.nc')
    comp = dict(zlib=True, 
                complevel=9, 
                least_significant_digit=2, 
                chunksizes=chunksize)
    one_output=output.format(sub_names[0])
    two_output=output.format(sub_names[1])
    print('Save splitted {0} as {1} and {2}'.format(
            name,one_output,two_output))
    flow_one.load().to_netcdf(one_output,
                 encoding={one_output.name:comp})
    flow_two.load().to_netcdf(two_output,
                 encoding={two_output.name:comp})    
    return one_output,two_output
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

def calc_residential_water_consumption(population_tif,area_tif,
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
    residential_wc+=blue_water*0 #replace nan in blue_water     
    # Save total_blue_water
    residential_wc=residential_wc.transpose('time',
                                                 'latitude','longitude')
    residential_wc.name=blue_water.name
    residential_wc.attrs={'units':blue_water.attrs['units'],
                      'source':'WCPC*Population/area',
                      'quantity':'residential water {0}'.format(flow_type)
                                  }
    if output is None:
        output=os.path.join(os.path.dirname(blue_water_nc),
                                     'residential_{0}.nc'.format(
                                             flow_type))
    comp = dict(zlib=True, 
                complevel=9, 
                least_significant_digit=2, 
                chunksizes=chunksize)
    print('Save residential water {0} as {1}'.format(
            flow_type,output))
    residential_wc.load().to_netcdf(output,
                         encoding={residential_wc.name:comp})    
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