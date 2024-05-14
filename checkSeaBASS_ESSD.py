#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

'Script to doube-check that netcdf and seabass files have same data'
    

@author: tomjordan - tjor@pml.ac.uk
"""

import os
import sys
import csv
import netCDF4 as nc
import pandas as pd

import numpy as np
import glob   
import pandas as pd
import ephem
import xarray as xr

import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.cm as cm


import matplotlib.patches as mpatches
import matplotlib.colors as cl

import math
import matplotlib.dates as mdates
import datetime as dt
from scipy.interpolate import interp1d
import scipy.io as io
import datetime


import pandas as pd
import glob

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import SB_support as sbs

from xml.dom.minidom import *

import scipy
from scipy import stats

from scipy import signal as sg


# general





def pigment_check(data_nc , data_hplc):
        
      plt.figure()
      plt.scatter(data_nc['hplc_Tot_Chl_a'], data_hplc['tot_chl_a'])
      plt.scatter(data_nc['hplc_DP'], data_hplc['dp'])
      # plt.scatter(data_nc['hplc_Phytin_a'], data_hplc['phytin_a'])                         
      plt.scatter(data_nc['hplc_Zea'], data_hplc['zea'])
      plt.scatter(data_nc['hplc_Allo'], data_hplc['allo'])
        
      return

def chl_check(data_nc,data_chl):
    
    plt.figure()
    plt.plot(data_chl['LAT'],data_chl['CHL'],label='old processing chain')
    # plt.plot(data_nc['uway_lat'],data_nc['acx_chl_debiased'],label='new processing chain')
    plt.plot(data_nc['uway_lat'],data_nc['acs_chl'],label='new processing chain',linestyle='dashed')
    plt.scatter(data_nc['hplc_lat'],data_nc['hplc_Tot_Chl_a'],label='hplc',color='red')
    plt.yscale('log')
    plt.legend()
   
    return


def N_ACS():
    print('AMT 19')
    print(len(~np.isnan(data_nc_19['acx_chl_debiased']))) 
    print(len(~np.isnan(data_nc_19['acs_ap'][:,20]))) 
    print('AMT 22')
    print(len(~np.isnan(data_nc_22['acs_chl_debiased']))) 
    print(len(~np.isnan(data_nc_22['acs_ap'][:,20]))) 
    print('AMT 23')
    print(len(~np.isnan(data_nc_23['acs_chl_debiased'])))  
    print(len(~np.isnan(data_nc_23['acs_ap'][:,20]))) 
    print('AMT 24')
    print(len(~np.isnan(data_nc_24['acs_chl_debiased'])))  
    print(len(~np.isnan(data_nc_24['acs_ap'][:,20]))) 
    print('AMT 25')
    print(len(~np.isnan(data_nc_25['acs_chl_debiased'])))  
    print(len(~np.isnan(data_nc_25['acs_ap'][:,20]))) 
    print('AMT 26')
    print(len(~np.isnan(data_nc_26['acs_chl_debiased'])))  
    print(len(~np.isnan(data_nc_26['acs_ap'][:,20]))) 
    print('AMT 27')
    print(len(~np.isnan(data_nc_27['acs_chl_debiased'])))  
    print(len(~np.isnan(data_nc_27['acs_ap'][:,20]))) 
    print('AMT 28')
    print(len(~np.isnan(data_nc_28['acx_chl_debiased'])))  
    print(len(~np.isnan(data_nc_28['acs_chl'])))  
   # print(sum(~np.isnan(data_nc_28['acs_ap'][:,20]))) 
    print('AMT 29')
    print(len(~np.isnan(data_nc_29['acs_chl_debiased'])))  
    print(len(~np.isnan(data_nc_29['acs_ap'][:,20]))) 
    return

def N_HPLC():
    print('AMT 22')
    print(len(data_nc_22['hplc_Tot_Chl_a'].values) - np.sum(data_nc_22['hplc_Indicate_if_filters_are_replicates'].values))
    print(np.max(data_nc_22['hplc_Tot_Chl_a'].values)) 
    print(np.min(data_nc_22['hplc_Tot_Chl_a'].values)) 
    print('AMT 23')
    print(len(data_nc_23['hplc_Tot_Chl_a'].values))
    print(np.max(data_nc_23['hplc_Tot_Chl_a'].values)) 
    print(np.min(data_nc_23['hplc_Tot_Chl_a'].values)) 
    print('AMT 24')
    print(len(data_nc_24['hplc_Tot_Chl_a'].values))
    print(np.max(data_nc_24['hplc_Tot_Chl_a'].values)) 
    print(np.min(data_nc_24['hplc_Tot_Chl_a'].values)) 
    print('AMT 25')
    print(len(data_nc_25['hplc_Tot_Chl_a'].values))
    print(np.max(data_nc_25['hplc_Tot_Chl_a'].values)) 
    print(np.min(data_nc_25['hplc_Tot_Chl_a'].values)) 
    print('AMT 26')
    print(len(data_nc_26['hplc_Tot_Chl_a'].values))
    print(np.max(data_nc_26['hplc_Tot_Chl_a'].values)) 
    print(np.min(data_nc_26['hplc_Tot_Chl_a'].values)) 
    print('AMT 27')
    print(len(data_nc_27['hplc_Tot_Chl_a'].values))
    print(np.max(data_nc_27['hplc_Tot_Chl_a'].values)) 
    print(np.min(data_nc_27['hplc_Tot_Chl_a'].values)) 
    print('AMT 28')
    print(len(data_nc_28['hplc_Tot_Chl_a'].values))
    print(np.max(data_nc_28['hplc_Tot_Chl_a'].values)) 
    print(np.min(data_nc_28['hplc_Tot_Chl_a'].values)) 
    return
    


def conc():
    ' + cruise + '
    print('AMT 19')
    print(np.percentile(data_nc_19['acx_chl_debiased'].values,0.05))
    print(np.percentile(data_nc_19['acx_chl_debiased'].values,99.5))
    print('AMT 22')
    print(np.percentile(data_nc_22['acs_chl_debiased'].values,0.05))
    print(np.percentile(data_nc_22['acs_chl_debiased'].values,99.5))
    print('AMT 23')
    print(np.percentile(data_nc_23['acs_chl_debiased'].values,0.05))
    print(np.percentile(data_nc_23['acs_chl_debiased'].values,99.5))
    print('AMT 24')
    print(np.nanpercentile(data_nc_24['acs_chl_debiased'].values,0.05))
    print(np.nanpercentile(data_nc_24['acs_chl_debiased'].values,99.5))
    print(np.nanpercentile(data_nc_24['acs2_chl_debiased'].values,0.05))
    print(np.nanpercentile(data_nc_24['acs2_chl_debiased'].values,99.5))
    print('AMT 25')
    print(np.nanpercentile(data_nc_25['acs_chl_debiased'].values,0.05))
    print(np.nanpercentile(data_nc_25['acs_chl_debiased'].values,99.5))
    print(np.nanpercentile(data_nc_25['acs2_chl_debiased'].values,0.05))
    print(np.nanpercentile(data_nc_25['acs2_chl_debiased'].values,99.5))
    print('AMT 26')
    print(np.percentile(data_nc_26['acs_chl_debiased'].values,0.05))
    print(np.percentile(data_nc_26['acs_chl_debiased'].values,99.5))
    print('AMT 27')
    print(np.percentile(data_nc_27['acs_chl_debiased'].values,0.05))
    print(np.percentile(data_nc_27['acs_chl_debiased'].values,99.5))
    print('AMT 28')
    print(np.percentile(data_nc_28['acx_chl_debiased'].values,0.05))
    print(np.percentile(data_nc_28['acx_chl_debiased'].values,99.5))
    print('AMT 29')
    print(np.percentile(data_nc_19['acx_chl_debiased'].values,0.05))
    print(np.percentile(data_nc_19['acx_chl_debiased'].values,99.5))

    return


def check_matchup_stats(data_nc):

      delta  =data_nc['acs_chl_debiased'].delta
      sigma = data_nc['acs_chl_debiased'].sigma 
      N = np.sum(~np.isnan(data_nc['acs_chl_debiased'].HPLC_Tot_chla - data_nc['acs_chl_debiased'].acs_chl))
      print(delta)
      print(sigma)
      print(N)                

      return
 
    
def check_matchup_stats_acx(data_nc):

      delta  =data_nc['acx_chl_debiased'].delta
      sigma = data_nc['acx_chl_debiased'].sigma 
      N = np.sum(~np.isnan(data_nc['acx_chl_debiased'].HPLC_Tot_chla - data_nc['acx_chl_debiased'].acx_chl))
      print(delta)
      print(sigma)
      print(N)                

      return


def check_matchup_stats_2ACS(data_nc):

      delta  =data_nc['acs_chl_debiased'].delta
      sigma = data_nc['acs_chl_debiased'].sigma 
      N = np.sum(~np.isnan(data_nc['acs_chl_debiased'].HPLC_Tot_chla - data_nc['acs_chl_debiased'].acs_chl))
      N2 = np.sum(~np.isnan(data_nc['acs2_chl_debiased'].HPLC_Tot_chla - data_nc['acs2_chl_debiased'].acs_chl))
      print(delta)
      print(sigma)
      print(N)  
      print(N2)              

      return


def chl_check_AMT25(data_nc,data_chl):
    
    plt.figure(figsize=(12,12))
    plt.plot(data_chl['LAT'],data_chl['CHL'],label='old processing chain')
    # plt.plot(data_nc['uway_lat'],data_nc['acx_chl_debiased'],label='new processing chain')
    plt.plot(data_nc['uway_lat'],data_nc['acs_chl'],label='new processing: ACS 1',linestyle='dashed')
    plt.plot(data_nc['uway_lat'],data_nc['acs2_chl'],label='new processing: ACS 2 ',linestyle='dashed')
    plt.scatter(data_nc['hplc_lat'],data_nc['hplc_Tot_Chl_a'],label='hplc',color='red')
    plt.yscale('log')
    plt.legend()
    plt.ylabel('Tot Chl-a: mg m$^{-3}$')
    plt.xlabel('Lat: deg')

    plt.figure(figsize=(12,12))
    plt.plot(data_chl['LAT'],data_chl['CHL'],label='old processing chain')
    # plt.plot(data_nc['uway_lat'],data_nc['acx_chl_debiased'],label='new processing chain')
    plt.plot(data_nc['uway_lat'],data_nc['acs_chl'],label='new processing: ACS 1',linestyle='dashed')
    plt.plot(data_nc['uway_lat'],data_nc['acs2_chl'],label='new processing: ACS 2 ',linestyle='dashed')
    plt.scatter(data_nc['hplc_lat'],data_nc['hplc_Tot_Chl_a'],label='hplc',color='red')
    plt.yscale('log')
    plt.legend()
    plt.ylabel('Tot Chl-a: mg m$^{-3}$')
    plt.xlabel('Lat: deg')
    plt.xlim(20,25)


    return


def check_hplc_fields(data_nc, sb_hplc):

    keys_nc = ['hplc_Tot_Chl_a','hplc_Tot_Chl_b', 'hplc_Tot_Chl_c' , 
                'hplc_Allo','hplc_Alpha-beta-Car', 'hplc_But-fuco', 
                'hplc_Diadino','hplc_Diato','hplc_Fuco',
                'hplc_Hex-fuco','hplc_Perid', 'hplc_Zea', 
                'hplc_PPC', 'hplc_PSC',
                'hplc_lat',  'hplc_lon', 'hplc_depth']
    
         
    keys_sb = ['tot_chl_a','tot_chl_b', 'tot_chl_c',   
                'allo','alpha-beta-car', 'but-fuco', 
                'diadino','diato','fuco',
                'hex-fuco','perid', 'zea', 
                'ppc', 'psc', 
                'lat',  'lon', 'depth']
    
        
    for i in range(len(keys_nc)):
        if np.array_equal(data_nc[keys_nc[i]], sb_hplc[keys_sb[i]]) == True:
            print('array match')                            
        elif np.array_equal(data_nc[keys_nc[i]], sb_hplc[keys_sb[i]]) == False:
            print('array non match')
            breakpoint()
            
    return


def check_acs_fields(data_nc, sb_acs):


    keys_nc = ['acs_chl_debiased_nomedfilt', 
                'uway_lat',
                'uway_lon',
                'acs_ap',
                'acs_ap_u',
                'acs_cp',
                'acs_cp_u',
                'acs_ap',
                'acs_ap_u',
                'acs_cp',
                'acs_cp_u',
            ]
    
        
    keys_sb = ['chl_lineheight', 
                'lat',
                'lon',
                'ap400.0',
                'ap400.0_unc',
                'cp400.0',
                'cp400.0_unc',
                'ap750.0',
                'ap750.0_unc',
                'cp750.0',
                'cp750.0_unc',
            ]
            
        
                    
    for i in range(len(keys_nc)):
        if i < 3:
            if np.array_equal(data_nc[keys_nc[i]], sb_acs[keys_sb[i]]) == True:
                print('array match')                            
            elif np.array_equal(data_nc[keys_nc[i]], sb_acs[keys_sb[i]]) == False:
                print('array non match')
                breakpoint()
        if i > 2 and i < 7:
            if np.array_equal(data_nc[keys_nc[i]].data[:,0], sb_acs[keys_sb[i]]) == True:
               print('array match')                            
            elif np.array_equal(data_nc[keys_nc[i]].data[:,0], sb_acs[keys_sb[i]]) == False:
             print('array non match')  
             breakpoint()
        if i > 6:
            if np.array_equal(data_nc[keys_nc[i]].data[:,-1], sb_acs[keys_sb[i]]) == True:
                print('array match')                            
            elif np.array_equal(data_nc[keys_nc[i]].data[:,-1], sb_acs[keys_sb[i]]) == False:
              print('array non match')  
              breakpoint()

    return


if __name__ == '__main__':
    

    ## load netcdf files ##
    # AMT 19 #
    fn_nc_19 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT19/Processed/Step3/amt19_final_with_debiased_chl.nc'
    data_nc_19 = xr.open_dataset(fn_nc_19)
    print(sorted(list(data_nc_19.keys())))

    # AMT 22 #
    fn_nc_22 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT22/Processed/Step3/amt22_final_with_debiased_chl.nc'
    data_nc_22 = xr.open_dataset(fn_nc_22)
    print(sorted(list(data_nc_22.keys())))
    
    # AMT 23 #
    fn_nc_23 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT23/Processed/underway/Step3/amt23_final_with_debiased_chl.nc'
    data_nc_23 = xr.open_dataset(fn_nc_23)
    print(sorted(list(data_nc_23.keys()))) 
    
    # AMT 24 #
    fn_nc_24 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT24/Processed/Uway/Step3/amt24_final_with_debiased_chl.nc'
    data_nc_24 = xr.open_dataset(fn_nc_24)
    print(sorted(list(data_nc_24.keys()))) 
      
    # AMT 25 #
    fn_nc_25 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT25/Processed/UWay/Step3/amt25_final_with_debiased_chl.nc'
    data_nc_25 = xr.open_dataset(fn_nc_25)
    print(sorted(list(data_nc_25.keys()))) 
    
    # AMT 26 #
    fn_nc_26 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT26/Processed/Underway/Step3/amt26_final_with_debiased_chl.nc'
    data_nc_26 = xr.open_dataset(fn_nc_26)
    print(sorted(list(data_nc_26.keys()))) 

    # AMT 27 #
    fn_nc_27 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT27/Processed/Underway/Step3/amt27_final_with_debiased_chl.nc'
    data_nc_27 = xr.open_dataset(fn_nc_27)
    print(sorted(list(data_nc_27.keys()))) 

    # AMT 28 #  
    fn_nc_28 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT28/Processed/Underway/Step3/amt28_final_with_debiased_chl.nc'
    data_nc_28 = xr.open_dataset(fn_nc_28)
    print(sorted(list(data_nc_28.keys())))   
     
    # AMT 29 #   
    fn_nc_29 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT29/Processed/Step3/amt29_final_with_debiased_chl.nc'
    data_nc_29 = xr.open_dataset(fn_nc_29)
    print(sorted(list(data_nc_29.keys())))   


    
    # AMT 27 - seabass and nc data (for cross-check)
    dir_sb_27 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT27/Source/SeaBASS_submit/sb_processed/'
    fn_sb_acs_27 = dir_sb_27 + 'AMT27_InLine_ACS_20170924_20171101_Particulate_v20240305.sb'
    fn_sb_hplc_27 = dir_sb_27 + 'AMT27_HPLC_20170924_20171101_v20240305.sb'
   # fn_chl_27 = '/data/datasets/cruise_data/active/ACS_Chl/amtacs/AMT27_ACS_CHL-A_MEDFILT_BIAS_CORRECTED.csv'

    
    # load acs sb
    sb_acs_27 = sbs.readSB(fn_sb_acs_27)
    print(sb_acs_27.data.keys())
    sb_acs_27 = sb_acs_27.data
    
     # load hplc sb
    sb_hplc_27 = sbs.readSB(fn_sb_hplc_27) 
    print(sb_hplc_27.data.keys())
    sb_hplc_27 = sb_hplc_27.data

    check_hplc_fields(data_nc_27, sb_hplc_27)
    check_acs_fields(data_nc_27, sb_acs_27)


    ###########################################################################
    # AMT 26 - seabass and nc data (for cross-check)
    # dir_sb_26 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT26/Source/SeaBASS_submit/sb_processed/'
    #  fn_sb_acs_26 = dir_sb_26 + 'AMT26_InLine0_ACS_20160923_20161102_Particulate_v20230622.sb'
    #  fn_sb_hplc_26 = dir_sb_26 + 'AMT26_HPLC_20160923_20161102_v20230622.sb'
      
    ###########################################################################
    # AMT 28 - seabass and nc data (for cross-check)     
    # dir_sb_28 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT28/Source/SeaBASS_submit/sb_processed/'
    # fn_sb_acs_28 =  dir_sb_28 + 'AMT28_InLine0_ACS_20180925_20181028_Particulate_v20230622.sb'
    # fn_sb_hplc_28 =  dir_sb_28 + 'AMT28_HPLC_20180925_20181027_v20230622.sb'
    #  fn_nc_28 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT28/Processed/Underway/Step3/amt28_final_with_debiased_chl.nc'
    # fn_chl_28 = '/data/datasets/cruise_data/active/ACS_Chl/amtacs/AMT28_ACS_CHL-A_MEDFILT_BIAS_CORRECTED.csv'
   
    
    ###########################################################################
    # AMT 29 - seabass and nc data (for cross-check)     
    #fn_nc_29 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT29/Processed/Step3/amt29_final_with_debiased_chl.nc'
   # data_nc_29 = xr.open_dataset(fn_nc_29)
   # data_nc_29.keys()
                               
