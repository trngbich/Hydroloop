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

def main(BASIN,unit_conversion=1):
    '''
    unit_conversion: 1 for TCM, 1000 for MCM, 1e6 for BCM or km3)
    '''
    #check requirements
    requirements=gd.get_requirements_for_sheet(3)    
    if not cf.check_requirement_sheet(BASIN,requirements):
        sys.exit("ERROR: Data requirements for Sheet 3 are not fulfilled")
    #create folder to save intermetidate data
    folder=os.path.join(BASIN['output_folder'],
                        'csv','timeseries')        
    if not os.path.exists(folder):
        os.makedirs(folder)
    output_file=os.path.join(folder,'sheet3_{0}.csv')
    #create folder to save sheet 3 csv
    sheet_folder=os.path.join(BASIN['output_folder'],'csv','sheet3') 
    if not os.path.exists(folder):
        os.makedirs(sheet_folder) 
        
    #Calulate yearly data to fill in Sheet 3    
        
    #Fill data in Sheet 3 csv
    for i in range(len(ET)):
        year=ET.index[i]
        results=dict()
        results['ET']=np.array(
                (ET.iloc[i].values/unit_conversion))                

        #write sheet 2 csv
        output_fh=os.path.join(sheet_folder,'sheet3_{0}.csv'.format(year))
        create_sheet2_csv(year,results,output_fh)    
        
    return True
        
def create_sheet4_csv():
    return