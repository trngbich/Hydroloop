# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 17:10:28 2020

@author: ntr002
"""
import os
import csv
import numpy as np
from . import calculate_flux as cf
from . import get_dictionaries as gd
from . import hydroloop as hl

def main(BASIN,unit_conversion=1000):
    '''
    unit_conversion: 1 for TCM, 1000 for MCM, 1e6 for BCM or km3)
    '''
    #check requirements
#    requirements=gd.get_sheet_requirements(2)    
#    if not cf.check_requirement_sheet(BASIN,requirements):
#        print("ERROR: Data requirements for Sheet 2 are not fulfilled")
#        return None
    #create folder to save intermetidate data
    folder=os.path.join(BASIN['output_folder'],
                        'csv','timeseries')        
    if not os.path.exists(folder):
        os.makedirs(folder)
    output_file=os.path.join(folder,'sheet2_{0}.csv')
    
    #Calulate yearly data to fill in Sheet 2    
    ET=cf.calc_flux_per_LU_class(BASIN['data_cube']['monthly']['et'], 
                         BASIN['data_cube']['monthly']['lu'], 
                         BASIN['gis_data']['basin_mask'],
                 chunksize=BASIN['chunksize'], #option to process in chunks
                 output=output_file.format('lu_et_monthly'), 
                 #option to save output as csv                  
                 quantity='volume')
    E=cf.calc_flux_per_LU_class(BASIN['data_cube']['monthly']['e'], 
                         BASIN['data_cube']['monthly']['lu'], 
                         BASIN['gis_data']['basin_mask'],
                 chunksize=BASIN['chunksize'], #option to process in chunks
                 output=output_file.format('lu_e_monthly'), 
                 #option to save output as csv                     
                 quantity='volume')
    T=cf.calc_flux_per_LU_class(BASIN['data_cube']['monthly']['t'], 
                         BASIN['data_cube']['monthly']['lu'], 
                         BASIN['gis_data']['basin_mask'],
                 chunksize=BASIN['chunksize'],
                 output=output_file.format('lu_t_monthly'),                               
                 quantity='volume')
    I=cf.calc_flux_per_LU_class(BASIN['data_cube']['monthly']['i'], 
                         BASIN['data_cube']['monthly']['lu'], 
                         BASIN['gis_data']['basin_mask'],
                 chunksize=BASIN['chunksize'],
                 output=output_file.format('lu_i_monthly'),               
                 quantity='volume') 
    sheet_folder=os.path.join(BASIN['output_folder'],'csv','sheet2') 
    if not os.path.exists(sheet_folder):
        os.makedirs(sheet_folder) #create sheet1 folder 
        
    #Fill data in Sheet 2 csv
    monthly_csvs=[]
    for i in range(len(E)):
        t_index=E.index[i]
        year=t_index.year
        month=t_index.month
        results=dict()
        results['LULC']=np.array(
                [float(s) for s in ET.columns])
        results['ET']=np.array(
                (ET.loc[t_index].values/unit_conversion))                
        results['E']=np.array(
                E.loc[t_index].values/unit_conversion)
        results['T']=np.array(
                T.loc[t_index].values/unit_conversion)
        results['I']=np.array(
                I.loc[t_index].values/unit_conversion)
        #write sheet 2 csv
        output_fh=os.path.join(sheet_folder,'sheet2_{0}_{1}.csv'.format(year,month))
        create_sheet2_csv(results,output_fh) 
        monthly_csvs.append(output_fh)
    ##calculate yearly sheets
    yearly_folder=os.path.join(sheet_folder,'yearly') 
    if not os.path.exists(yearly_folder):
        os.makedirs(yearly_folder) #create sheet1 folder  
    yearly_csvs=hl.calc_yearly_sheet(monthly_csvs,
                                     yearly_folder,
                                     hydroyear=BASIN['hydroyear'])
    return yearly_csvs
        
def create_sheet2_csv(results,output_fh):
    """
    Create the csv-file needed to plot sheet 1.
    
    Parameters
    ----------
    results : dict
        Dictionary generated by calc_sheet1.
    output_fh : str
        Filehandle to store the csv-file.
    """
    classes_dict=gd.get_sheet2_classes()  
    lulc_dict=gd.get_lulcs(lulc_version = '4.0')
    first_row = ['LAND_USE', 'CLASS', 'TRANSPIRATION', 'WATER',
                 'SOIL', 'INTERCEPTION', 'AGRICULTURE',
                 'ENVIRONMENT', 'ECONOMY', 'ENERGY', 'LEISURE',
                 'NON_BENEFICIAL']
    
    if not os.path.exists(os.path.split(output_fh)[0]):
        os.makedirs(os.path.split(output_fh)[0])
    
    csv_file = open(output_fh, 'w')
    writer = csv.writer(csv_file, delimiter=';', lineterminator = '\n')
    writer.writerow(first_row)
    
    for LAND_USE in list(classes_dict.keys()):
            for CLASS in list(classes_dict[LAND_USE].keys()):
                print(LAND_USE)
                print(CLASS)
                write_sheet2_row(LAND_USE, CLASS, 
                                 lulc_dict, 
                                 classes_dict,
                                 results,
                                 writer)
    
    csv_file.close()
    
def write_sheet2_row(LAND_USE, CLASS, 
                     lulc_dict, 
                     classes_dict,
                     results, 
                     writer):
    """
    Write a row with spatial aggregates to a sheet2 csv-file.
    
    Parameters
    ----------

    """
    LULC=results['LULC']
    ET=results['ET']
    E=results['E']
    T=results['T']
    I=results['I']
    # Get a list of the different landuse classes to be aggregated.
    lulcs = classes_dict[LAND_USE][CLASS]
    
    # Create a mask to ignore non relevant pixels.
    mask=np.logical_or.reduce([LULC == value for value in lulcs])
    
    # Calculate the spatial sum of the different parameters.
    evapotranspiration =np.nansum(ET[mask])
    transpiration = np.nansum(T[mask])
    interception = np.nansum(E[mask])
    evaporation = np.nansum(I[mask])
    sumcheck=transpiration + interception + evaporation
    if evapotranspiration != sumcheck:
        print('ETI split difference: {0}'.format(
                evapotranspiration-sumcheck))
    # Set special cases.
    if np.any([CLASS == 'Natural water bodies', CLASS == 'Managed water bodies']):
        soil_evaporation = 0
        water_evaporation = evaporation
    else:
        soil_evaporation = evaporation
        water_evaporation = 0            
        
    # Create some necessary variables.
    agriculture = 0
    environment = 0
    economy = 0
    energy = 0
    leisure = 0
    non_beneficial = 0
    
    # Calculate several landuse type specific variables.
    for lu_type in lulcs:
        
        # Get some constants for the landuse type.
        beneficial_percentages = np.array(lulc_dict[lu_type][3:6]) / 100
        service_contributions = np.array(lulc_dict[lu_type][6:11]) / 100         
          
        # Calculate the beneficial ET.
        benef_et = np.nansum([np.nansum(T[LULC == lu_type]) * beneficial_percentages[0],
               np.nansum(E[LULC == lu_type]) * beneficial_percentages[1],
               np.nansum(I[LULC == lu_type]) * beneficial_percentages[2]])
               
        # Determine the service contributions.
        agriculture += benef_et * service_contributions[0] 
        environment += benef_et * service_contributions[1] 
        economy += benef_et * service_contributions[2]
        energy += benef_et * service_contributions[3]
        leisure += benef_et * service_contributions[4]
       
        # Determine non-beneficial ET.
        non_beneficial += (np.nansum([np.nansum(T[LULC == lu_type]) * (1 - beneficial_percentages[0]),
               np.nansum(E[LULC == lu_type]) * (1 - beneficial_percentages[1]),
               np.nansum(I[LULC == lu_type]) * (1 - beneficial_percentages[2])]))
    
    # Create the row to be written
    row = [LAND_USE, CLASS, 
           "{0}".format(np.nansum([0, transpiration])),
           "{0}".format(np.nansum([0, water_evaporation])),
           "{0}".format(np.nansum([0, soil_evaporation])),
           "{0}".format(np.nansum([0, interception])),
           "{0}".format(np.nansum([0, agriculture])),
           "{0}".format(np.nansum([0, environment])),
           "{0}".format(np.nansum([0, economy])),
           "{0}".format(np.nansum([0, energy])),
           "{0}".format(np.nansum([0, leisure])),
           "{0}".format(np.nansum([0, non_beneficial]))]
    
    # Write the row.
    writer.writerow(row)