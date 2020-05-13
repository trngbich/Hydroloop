# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 16:50:07 2020

@author: ntr002
"""
import numpy as np
import pandas as pd
import xarray as xr
import os
import GIS_functions as gis

def open_nc(input_nc, chunksize=None, layer=None):
    if chunksize is None:
        dts=xr.open_dataset(input_nc)
    else:
        dts=xr.open_dataset(input_nc,
                            chunks={'time':chunksize[0],
                                    'longitude':chunksize[1],
                                    'latitude':chunksize[2]})
    if type(layer) is int: #select dataarray by index
        layer_name=list(dts.keys())[layer]
        return dts[layer_name]
    elif type(layer) is str: #select dataarray by name
        return dts[layer]
    else:
        return dts
    
def create_yearly_dataset(monthly_nc, output=None, 
                          hydroyear='A-DEC',chunksize=None):
    '''
     create yearly dataset from monthly dataset
     monthly_nc: str 
         path to monthly NetCDF file
     output: str 
         path to output yearly NetCDF file
     hydroyear: str 
         End month of hydrological year. default is 'A-DEC' for December
    chunksize: list of 3 int
        time chunk, x and y chunks
        default is None mean not using chunks
    '''
    #check output path
    if output is None:
        if monthly_nc[-10:] == 'monthly.nc':
            output=monthly_nc.replace('monthly.nc','yearly.nc')
        else:
            output=monthly_nc.replace('.nc','_yearly.nc')
    else:
        if not os.path.exists(os.path.dirname(output)):
            os.makedirs(os.path.dirname(output))
    #open monthly nc
    dts=open_nc(monthly_nc, chunksize=chunksize)
    #resample to hydrological years
    dts_y=dts.resample(time=hydroyear).sum(dim=['time'],skipna=False)
    dts_y.to_netcdf(output)
    
    #close netcdf file
    dts.close()
    dts_y.close()
    
    return output

def calc_flux_per_basin(dts_nc, basin_mask, chunksize=None,output=None):
    '''
    calculate flux per basin/sub-basin
    input_nc: str
        path to dataset (NetCDF) (in mm)
    basin_mask: str OR np.array/xr.DataArray
        str 
            path to basin mask (GeoTIFF)
        np.array/xr.DataArray 
            pixel area of basin (in km2)
    output: str
        path to output (csv)
        default is None 
    
    return
        dataframe (in TCM)
    '''
    dts=open_nc(dts_nc,chunksize=chunksize)
    #read area mask
    if type(basin_mask) is str:
        basin=gis.OpenAsArray(basin_mask,nan_values=True)
        area_map=gis.MapPixelAreakm(basin_mask)
        area_mask=area_map*basin
    else:
        area_mask=basin_mask
    #calculate flux map
    dts_m=dts*area_mask #flux = depth*area                
    df=dts_m.sum(dim=['latitude','longitude']).to_dataframe() #export data
    if output is not None:
        df.to_csv(output,sep=';') #save data as csv
        print('Save basin flux as {0}'.format(output))
    dts.close()
    return df


def aggregate_by_lu_unique(dts,LU,how='sum'):
    '''aggregate dataset by unique LU classes in LU map(s)
    '''
    unique_LU=np.unique(LU) #get unique landuse classes
    unique_LU=unique_LU[~np.isnan(unique_LU)] #exclude nan
    data=[] #create empty data list
    for lucl in unique_LU: #agrregate total fluxes per each lu class
        dts_lu=dts.where(LU==lucl,np.nan) #mask only lu class
        if how=='sum':
            df_lu=dts_lu.sum(dim=[
                    'latitude',
                    'longitude'
                    ]).to_dataframe() #sum of all pixels in lu class
        elif how=='mean':
            df_lu=dts_lu.mean(dim=[
                    'latitude',
                    'longitude'
                    ]).to_dataframe() #mean of all pixels in lu class          
        df_lu=df_lu.drop(columns='time')        
        df_lu.columns=['{0}-{1}'.format(lucl,
                       col) for col in df_lu.columns] #rename column        

        data.append(df_lu) #append data list by lu class
    df=pd.concat(data, axis=1) #merge all results into 1 dataframe
    return df

def aggregate_by_lu_dictionary(dts,LU,lu_dictionary,how='sum'):
    '''aggregate dataset by LU classes categories 
    '''
    data=[] #create empty data list
    for key in lu_dictionary: #agrregate total fluxes per each lu class
        classes=lu_dictionary[key]
        dts_lu=dts.where(LU.isin(classes),np.nan) #mask only lu class
        if how=='sum':
            df_lu=dts_lu.sum(dim=[
                    'latitude',
                    'longitude'
                    ]).to_dataframe() #sum of all pixels in lu class
        elif how=='mean':
            df_lu=dts_lu.mean(dim=[
                    'latitude',
                    'longitude'
                    ]).to_dataframe() #mean of all pixels in lu class            
        df_lu=df_lu.drop(columns='time')        
        df_lu.columns=['{0}-{1}'.format(key,
                       col) for col in df_lu.columns] #rename column
        df_lu=df_lu.reset_index(drop=True, inplace=True)
        data.append(df_lu) #append data list by lu class
    df=pd.concat(data, axis=1) #merge all results into 1 dataframe
    return df

def calc_flux_per_LU_class(dts_nc, lu_nc, basin_mask,
                     chunksize=None, #option to process in chunks
                     output=None, #option to save output as csv                     
                     lu_dictionary=None, #option to calculate by LU category
                     quantity='volume'): #option to calculate L or L^3
    '''
    calculate flux per LU class in WA+ LU map
    input_nc: str
        path to yearly dataset (NetCDF) (in mm)
    lu_nc: str
        path to NetCDF of LULC map
        
    basin_mask: str OR np.array/xr.DataArray
        str 
            path to basin mask (GeoTIFF)
        np.array/xr.DataArray 
            pixel area of basin (in km2)
    output: str
        path to output (csv)
        default is None 
    quantity: str
        volume or depth
    
    return
        dataframe (in TCM)
    '''
    dts=open_nc(dts_nc,chunksize=chunksize)
    lu=open_nc(lu_nc,chunksize=chunksize,layer=0)
    
    #read basin mask
    if type(basin_mask) is str:
        basin=gis.OpenAsArray(basin_mask,nan_values=True)
        
    if quantity=='volume':
        #get area mask
        if type(basin_mask) is str:
            area_map=gis.MapPixelAreakm(basin_mask)
            area_mask=area_map*basin
        else:
            area_mask=basin_mask 
        dts_m=dts*area_mask #flux = depth*area   
        dts_y=dts_m.groupby('time.year').sum(dim=['time'],skipna=False)
        method='sum'
        
    elif quantity=='depth':
        dts_m=dts*basin
        dts_y=dts_m.groupby('time.year').sum(dim=['time'],skipna=False)
        method='mean'
        
    n_lu=len(lu.time) #number of landuse map
    
    if n_lu==1: #single landuse map
        LU=lu[0] #get single landuse map            
    elif n_lu==len(dts_y.time): #dynamic yearly landuse  maps    
        LU=lu.groupby('time.year').first()
    
    if lu_dictionary is None:
        df=aggregate_by_lu_unique(dts_y,LU,how=method)
    elif type(lu_dictionary) is dict:
        df=aggregate_by_lu_dictionary(dts_y,LU,lu_dictionary,
                                      how=method)
        
    if output is not None: #export result if output path is defined
        df.to_csv(output,sep=';')
        print('Save LU flux as {0}'.format(output))
    lu.close()
    dts.close()
    return df

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
                output=create_yearly_dataset(
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
