# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 16:50:07 2020

@author: ntr002
"""
import numpy as np
import pandas as pd
import xarray as xr
import os
from . import GIS_functions as gis

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

def resample_to_monthly_dataset(yearly_nc, sample_nc,
                                start_month=0,
                                output=None,
                                chunksize=None):
    '''
    yearly_nc: a yearly dataset to resample to monthly
    sample_nc: a monthly dataset to sample 'time' dimension
    start_month: the index of start month. 
        default start_month = 0 means yearly value resample from the first
        month in sample_nc
                
    Resample a yearly netCDF dataset to monthly netCDF dataset
    Where the value of each month is the same with the value of the year
    '''
    dts1=open_nc(yearly_nc,chunksize=chunksize,layer=0)
    dts2=open_nc(sample_nc,chunksize=chunksize,layer=0)  
   
    for i in range(len(dts1.time)):
        for t in range(i*12-start_month,i*12+12-start_month):
            LU=dts1.isel(time=i)
            if t==0:
                dts=LU
            else:
                dts = xr.concat([dts, LU], dim='time')  
    dts['time']=dts2['time']
    #change coordinates order to [time,latitude,longitude]
    dts=dts.transpose('time','latitude','longitude')               

    dts.attrs=dts1.attrs
    dts.name = dts1.name
    
    comp = dict(zlib=True, 
                complevel=9, 
                least_significant_digit=2, 
                chunksizes=chunksize)
    if output is None:
        output=yearly_nc.replace('.nc','_resampled_monthly.nc')
    encoding = {dts.name: comp}
    dts.load().to_netcdf(output,encoding=encoding)
    dts1.close()
    dts2.close()
    print('Save monthly LU datacube as {0}'.format(output))
    return output

def calc_flux_per_basin(dts_nc, basin_mask, 
                        chunksize=None,output=None,quantity='volume'):
    '''
    calculate flux per basin/sub-basin
    input_nc: str
        path to dataset (NetCDF) (in mm)
    basin_mask: str OR np.array/xr.DataArray
        str 
            path to basin mask (GeoTIFF)
        np.array/xr.DataArray 
            pixel area of basin (in km2)
            or mask of basin
    quantity: str
        'volume' OR 'depth'
        
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
    else: #basin_mask is 2D array
        area_mask=basin_mask
    #calculate flux 
    if quantity=='volume':
        dts_m=dts*area_mask #flux = depth*area                
        df=dts_m.sum(dim=['latitude','longitude']).to_dataframe() #export data
    elif quantity=='depth':
        dts_m=dts*basin_mask
        df=dts_m.mean(dim=['latitude','longitude']).to_dataframe() #export data
    if output is not None:
        df.to_csv(output,sep=';') #save data as csv
        print('Save basin flux as {0}'.format(output))
    dts.close()
    return df


def calc_flux_per_LU_class(dts_nc, lu_nc, basin_mask,
                     chunksize=None, 
                     output=None,            
                     lu_dictionary=None, 
                     quantity='volume'):
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
    chunksize: list
        [t,x,y] chunksize of datacube
    output: str
        path to output (csv)
        default is None 
    quantity: str
        'volume': multiplied with pixel area (km2). 
                 if variable unit is mm, 
                 the result will be in 10^3 m3 or 10^-6 km3
                 if the variable unit is kg/ha,
                 the result will be in 10^2 kg
        'depth': keep the variable unit
    
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
        method='sum'
        
    elif quantity=='depth':
        dts_m=dts*basin
        method='mean'
        
#    n_lu=len(lu.time) #number of landuse map
    
#    if n_lu==1: #single landuse map
#        LU=lu[0] #get single landuse map            
    
    if lu_dictionary is None:
        df=aggregate_by_lu_unique(dts_m,lu,how=method)
    elif type(lu_dictionary) is dict:
        df=aggregate_by_lu_dictionary(dts_m,lu,lu_dictionary,
                                      how=method)
        
    if output is not None: #export result if output path is defined
        df.to_csv(output,sep=';')
        print('Save LU flux as {0}'.format(output))
    lu.close()
    dts.close()
    return df

def aggregate_by_lu_unique(dts,LU,how='sum'):
    '''aggregate dataset by unique LU classes in LU map(s)
    '''
    unique_LU=np.unique(LU) #get unique landuse classes
    unique_LU=unique_LU[~np.isnan(unique_LU)] #exclude nan
    data=[] #create empty data list
    
    LU=dts*0+LU #Trick: to keep same time dimension
    dts=LU*0+dts #Trick: to keep same time dimension
    for lucl in unique_LU: #agrregate total fluxes per each lu class
        dts_lu=xr.where(LU==lucl,dts,np.nan) #mask only lu class
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
        if len(df_lu.columns)>1:      
            df_lu.columns=['{0}-{1}'.format(lucl,
                           col) for col in df_lu.columns] #rename column with variable   
        else:
            df_lu.columns=['{0}'.format(lucl) for col in df_lu.columns] #rename column             
        data.append(df_lu) #append data list by lu class
    df=pd.concat(data, axis=1) #merge all results into 1 dataframe
    return df

def aggregate_by_lu_dictionary(dts,LU,lu_dictionary,how='sum'):
    '''aggregate dataset by LU classes categories 
    '''
    data=[] #create empty data list
    LU=dts*0+LU #Trick: to keep same time dimension
    dts=LU*0+dts #Trick: to keep same time dimension
    for key in lu_dictionary: #agrregate total fluxes per each lu class
        classes=lu_dictionary[key]
        dts_lu=xr.where(LU.isin(classes),dts,np.nan) #mask only lu class
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
        if len(df_lu.columns)>1:      
            df_lu.columns=['{0}-{1}'.format(key,
                           col) for col in df_lu.columns] #rename column with variable   
        else:
            df_lu.columns=['{0}'.format(key) for col in df_lu.columns] #rename column             
        data.append(df_lu) #append data list by lu class
    df=pd.concat(data, axis=1) #merge all results into 1 dataframe
    return df

