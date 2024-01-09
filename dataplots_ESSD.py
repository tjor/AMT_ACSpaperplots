#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 11:36:02 2023

@author: tjor
"""

# Load net-cdf formatted data

# Load Seabass formatted data

# Check data sets for consistency
### %run write_sb_file.py /Users/gdal/Documents/AMT_processed/AMT29/Step3/amt29_final_with_debiased_chl.nc acs122.dev checklist_acs_particulate_inline_AMT29.rtf,checklist_acs_ag_cg_AMT29.rtf,AMT29_ACS_inline_ProcessingReport_v20220810.docx

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
def nearest(items, pivot): 
    'function to locate nearest values and indicies in list x relative to a pivot (used to locate nearest timestamps)'
    
    nearest_value = min(items, key=lambda x: abs(x - pivot)) # finds nearest value    
    nearest_index= items.index(min(items, key=lambda x: abs(x - pivot))) # finds nearest index
   
    return nearest_value, nearest_index


## Longhust functions

def longhurst_mask_HPLC(lat, lon, cruise): 
    
    province_mask = []
    for i in range(len(lat)):
      print(i)
      province_mask.append(longhurst(lat[i], lon[i]))
      
    df_province = pd.DataFrame(province_mask, columns=['LH_Province'])
    df_province.to_csv('LH_mask_HPLC_AMT.csv')

    return  df_province
        

def longhurst_mask(lat, lon, cruise): 
    
    province_mask = []
    for i in range(len(lat)):
      print(i)
      province_mask.append(longhurst(lat[i], lon[i]))
      
    df_province = pd.DataFrame(province_mask, columns=['LH_Province'])
    df_province.to_csv('LH_mask_AMT' + cruise + '.csv')

    return  df_province
        


def longhurst(myLat, myLon):

    ## Get lat and lon from command line argument list
    ###--------------------------------------------------------------------------

    #ppFileName = string(sys.argv[1])
    #imgFileName = string(sys.argv[2])	
    	
    	
    ### Parse GML data from longhurst.xml
    ###--------------------------------------------------------------------------

        
        provinces = {}
        tree = parse('longhurst.xml')
        
        for node in tree.getElementsByTagName('MarineRegions:longhurst'):
        
        	# 1. Get province code, name and bounding box from file
        	provCode = node.getElementsByTagName('MarineRegions:provcode')[0].firstChild.data
        	provName = node.getElementsByTagName('MarineRegions:provdescr')[0].firstChild.data
        	fid = node.getAttribute("fid")
        	b = node.getElementsByTagName('gml:coordinates')[0].firstChild.data
        
        	# 2. Parse bounding box coordinates
        	b = b.split(' ')
        	x1,y1 = b[0].split(',')
        	x2,y2 = b[1].split(',')
        	x1 = float(x1)
        	y1 = float(y1)
        	x2 = float(x2)
        	y2 = float(y2)
        
        	provinces[fid] = {'provName': provName, 'provCode': provCode, 'x1': x1, 'y1': y1, 'x2': x2, 'y2': y2}
            
       #print(provinces)
        
        
        ### Find which candidate provinces our coordinates come from.
        ###--------------------------------------------------------------------------
        
        inProvince = {}
        for p in provinces:
        	inLat = 0
        	inLon = 0
        	
        	if (myLat>=provinces[p]['y1'] and myLat<=provinces[p]['y2']):
        		inLat = 1
        		
        	if (myLon>=provinces[p]['x1'] and myLon<=provinces[p]['x2']):
        		inLon = 1
        	
        	if inLat and inLon:
        		inProvince[p] = True	
        		
        		
        ### Perform Crossings Test on each candidate province.
        ###--------------------------------------------------------------------------
        
        for node in tree.getElementsByTagName('MarineRegions:longhurst'):
        	fid = node.getAttribute("fid")
        	
        	if inProvince.get(fid):
        		crossings = 0
        		
        		## 1. Get all coordinate pairs for this province.
        		geom = node.getElementsByTagName('MarineRegions:the_geom')
        		
        		for g in geom:
        			c = g.getElementsByTagName('gml:coordinates')
        			
        			for i in c:
        				ii = i.childNodes
        				coordStr = ii[0].data		#<--- contains coordinate strings
        				P = coordStr.split(' ')
        				
        				pairs = []
        				for p in P:
        					[lon,lat] = p.split(',')
        					pairs.append([float(lon),float(lat)])	
        					
        				## 2. Use pair p and p+1 to perform Crossings Test.
        				for i in range(len(pairs)-1):
        					# test latitude
        					passLat = (pairs[i][1]>=myLat and pairs[i+1][1]<=myLat) or (pairs[i][1]<=myLat and pairs[i+1][1]>=myLat)
        
        					# test longitude
        					passLon = (myLon <= pairs[i+1][0])
        				
        					if passLon and passLat:
        						crossings += 1
        		 					
        		if crossings%2==1: 
        			inProvince[fid] = True
        		else:
        			inProvince[fid] = False
        

        solution = []
        for i in inProvince:
        	if inProvince[i] == True:
        		solution.append(provinces[i]['provCode'])# provinces[i]['provName']])
              #  solution.append(provinces[i]['provCode'])
              
        return solution
    
# coverage maps

def plot_coverage():

    plt.figure(figsize=(15,10))
    extent = [-65,20,-55, 60]
    request = cimgt.GoogleTiles(style='satellite')
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.add_image(request, 6)
    ax.set_extent(extent, ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels_top = gl.right_labels = False
    gl.xformatter =  LONGITUDE_FORMATTER
    gl.yformatter =  LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 14,  'rotation':45}
    gl.ylabel_style = {'size': 14,  'rotation': 0}
    
    
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.tick_params(labelsize=10)
    
    from cartopy.io.shapereader import Reader
    from cartopy.feature import ShapelyFeature
    
    fname = '/users/rsg/tjor/scratch_network/tjor/AMT_underway/Tests_and_plots/longhurst_v4_2010/Longhurst_world_v4_2010.shp'
    
    shape_feature = ShapelyFeature(Reader(fname).geometries(),
                                    ccrs.PlateCarree(),facecolor='none')
    ax.add_feature(shape_feature)
    
    
    ax.scatter(data_nc_19['uway_lon'], data_nc_19['uway_lat'], s=0.5, color='red',label='AMT 19')
    ax.scatter(data_nc_22['uway_lon'], data_nc_22['uway_lat'], s=0.5, color='orange',label='AMT 22')
    ax.scatter(data_nc_23['uway_lon'], data_nc_23['uway_lat'], s=0.5, color='yellow',label='AMT 23')
    ax.scatter(data_nc_24['uway_lon'], data_nc_24['uway_lat'], s=0.5, color='yellowgreen',label='AMT 24')
    ax.scatter(data_nc_25['uway_lon'], data_nc_25['uway_lat'], s=0.5, color='lime',label='AMT 25')
    ax.scatter(data_nc_26['uway_lon'], data_nc_26['uway_lat'], s=0.5, color='cyan',label='AMT 26')
    ax.scatter(data_nc_27['uway_lon'], data_nc_27['uway_lat'], s=0.5, color='magenta',label='AMT 27')
    ax.scatter(data_nc_28['uway_lon'], data_nc_28['uway_lat'], s=0.5, color='grey',label='AMT 28')
    ax.scatter(data_nc_29['uway_lon'], data_nc_29['uway_lat'], s=0.5, color='white',label='AMT 29')
    plt.legend(markerscale=8)            
    
    return


# IOP plot functions

def _9cruise_IOPsum(field, field2):
    
    combined_array = np.concatenate((data_nc_19[field], data_nc_22[field], data_nc_23[field], 
                                     data_nc_24[field], data_nc_24[field2], 
                                     data_nc_25[field], data_nc_25[field2], 
                                     data_nc_26[field], data_nc_27[field], 
                                     data_nc_28[field], data_nc_29[field]), axis=0)
                    
    return combined_array           
            
def int_norm(data): 
    'integral normalization, (follows Boss 2013)' 
    data_in = np.nan*np.ones([len(data),len(data.T)])
    for i in range(len(data)):
        data_in[i,:] = data[i,:]/np.mean(data[i,:])
        
    return data_in

    
def plot_uncertainties():
        
    apu_mat = _9cruise_IOPsum('acs_ap_u', 'acs2_ap_u')
    bpu_mat = _9cruise_IOPsum('acs_bp_u', 'acs2_bp_u')
    cpu_mat = _9cruise_IOPsum('acs_cp_u', 'acs2_cp_u')
    
    ap_mat = _9cruise_IOPsum('acs_ap', 'acs2_ap')
    bp_mat = _9cruise_IOPsum('acs_bp', 'acs2_bp')
    cp_mat = _9cruise_IOPsum('acs_cp', 'acs2_cp')
        
    plt.figure(figsize=(15,9))
    plt.subplot(2,3,1)
    plt.rcParams.update({'font.size': 16})
    plt.plot(data_nc_29['acs_wv'],np.nanmedian(apu_mat,axis=0),label='Median', color='blue')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(apu_mat,25,axis=0), color='blue', linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(apu_mat,75,axis=0), label='Quartiles', color='blue', linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(apu_mat,90,axis=0), label='10$^{th}$ & 90$^{th}$ percentiles', color='blue',linestyle='dotted')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(apu_mat,10,axis=0), color='blue',linestyle='dotted')
    plt.xlim(400,720)
    plt.ylim(0,0.015)
    #plt.xlabel('Wavelength [nm]')
    plt.ylabel('Uncertainty [m$^{-1}$]')
    plt.legend(fontsize=12)
    ax = plt.gca()
    plt.text(.05, .95,  'A', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
    
    plt.subplot(2,3,2)
    plt.rcParams.update({'font.size': 16})
    plt.plot(data_nc_29['acs_wv'],np.nanmedian(bpu_mat,axis=0),label='Median', color='black')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(bpu_mat,25,axis=0), color='black', linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(bpu_mat,75,axis=0), label='Quartiles', color='black', linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(bpu_mat,90,axis=0), label='10$^{th}$ & 90$^{th}$ percentiles', color='black',linestyle='dotted')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(bpu_mat,10,axis=0), color='black', linestyle='dotted')
    plt.xlim(400,720)
    plt.ylim(0,0.015)
    #plt.xlabel('Wavelength [nm]')
    #plt.ylabel('Uncertainty [m$^{-1}$]')
    plt.legend(fontsize=12)
    ax = plt.gca()
    plt.text(.05, .95,  'B', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
    
    
    plt.subplot(2,3,3)
    plt.rcParams.update({'font.size': 16})
    plt.plot(data_nc_29['acs_wv'],np.nanmedian(cpu_mat,axis=0),label='Median', color='red')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(cpu_mat,25,axis=0), color='red', linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(cpu_mat,75,axis=0), label='Quartiles', color='red', linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(cpu_mat,90,axis=0), label='10$^{th}$ & 90$^{th}$ percentiles', color='red',linestyle='dotted')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(cpu_mat,10,axis=0), color='red',linestyle='dotted')
    plt.xlim(400,720)
    plt.ylim(0,0.015)
    #plt.xlabel('Wavelength [nm]')
    #plt.ylabel('Uncertainty [m$^{-1}$]')
    plt.legend(fontsize=12)
    ax = plt.gca()
    plt.text(.05, .95,  'C', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
    
    
    
    plt.subplot(2,3,4)
    plt.rcParams.update({'font.size': 16})
    plt.plot(data_nc_29['acs_wv'],np.nanmedian(100*apu_mat/np.abs(ap_mat),axis=0),label='Median', color='blue')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*apu_mat/np.abs(ap_mat),25,axis=0), color='blue', linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*apu_mat/np.abs(ap_mat),75,axis=0), label='Quartiles', color='blue', linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*apu_mat/np.abs(ap_mat),90,axis=0), label='10$^{th}$ & 90$^{th}$ percentiles', color='blue',linestyle='dotted')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*apu_mat/np.abs(ap_mat),10,axis=0), color='blue',linestyle='dotted')
    plt.xlim(400,750)
    plt.ylim(0,100)
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Percentage uncertainty [%]')
    ax = plt.gca()
    plt.text(.05, .95,  'D', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
    # plt.legend(fontsize=12)
    
    
    plt.subplot(2,3,5)
    plt.rcParams.update({'font.size': 16})
    plt.plot(data_nc_29['acs_wv'],np.nanmedian(100*bpu_mat/np.abs(bp_mat),axis=0),label='Median', color='black')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*bpu_mat/np.abs(bp_mat),25,axis=0), color='black', linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*bpu_mat/np.abs(bp_mat),75,axis=0), label='Quartiles', color='black', linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*bpu_mat/np.abs(bp_mat),90,axis=0), label='10$^{th}$ & 90$^{th}$ percentiles', color='black',linestyle='dotted')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*bpu_mat/np.abs(bp_mat),10,axis=0), color='black',linestyle='dotted')
    plt.xlim(400,750)
    plt.ylim(0,30)
    plt.xlabel('Wavelength [nm]')
    ax = plt.gca()
    plt.text(.05, .95,  'E', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
    #plt.ylabel('Percentage uncertainty [%]')
    #plt.legend(fontsize=12)
    
    
    plt.subplot(2,3,6)
    plt.rcParams.update({'font.size': 16})
    plt.plot(data_nc_29['acs_wv'],np.nanmedian(100*cpu_mat/np.abs(cp_mat),axis=0),label='Median', color='red')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*cpu_mat/np.abs(cp_mat),25,axis=0), color='red', linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*cpu_mat/np.abs(cp_mat),75,axis=0), label='Quartiles', color='red', linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*cpu_mat/np.abs(cp_mat),90,axis=0), label='10$^{th}$ & 90$^{th}$ percentiles', color='red',linestyle='dotted')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*cpu_mat/np.abs(cp_mat),10,axis=0), color='red',linestyle='dotted')
    plt.xlim(400,750)
    plt.ylim(0,30)
    plt.xlabel('Wavelength [nm]')
    ax = plt.gca()
    plt.text(.05, .95,  'F', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
    #plt.ylabel('Percentage uncertainty [%]')
    #plt.legend(fontsize=12)
    
    plt.tight_layout(pad=1.6)
    
    return


def plot_median_ap_province(data_nc, mask):
        ' plots median ap, bp, cp by province'
    
        ap = data_nc['acs_ap'].data
        bp = data_nc['acs_bp'].data
        cp = data_nc['acs_cp'].data
            
        ap_int = int_norm(ap)
        bp_int = int_norm(bp)
        cp_int = int_norm(cp)
        
        plt.figure(figsize=(15,9))
        plt.subplot(2,3,1)
        plt.plot(data_nc['acs_wv'],np.median(ap[mask=='NECS',:],axis=0),label='NECS')
        plt.plot(data_nc['acs_wv'],np.median(ap[mask=='NADR',:],axis=0),label='NADR')
        plt.plot(data_nc['acs_wv'],np.median(ap[mask=='NASE',:],axis=0),label='NASE')
        plt.plot(data_nc['acs_wv'],np.median(ap[mask=='NATR',:],axis=0),label='NATR')
        plt.plot(data_nc['acs_wv'],np.median(ap[mask=='WTRA',:],axis=0),label='WTRA')
        plt.plot(data_nc['acs_wv'],np.median(ap[mask=='SATL',:],axis=0),label='SATL')
        plt.plot(data_nc['acs_wv'],np.median(ap[mask=='SANT',:],axis=0),label='SANT')
        plt.plot(data_nc['acs_wv'],np.median(ap[mask=='SSTC',:],axis=0),label='SSTC')
        plt.plot(data_nc['acs_wv'],np.median(ap[mask=='FKLD',:],axis=0),label='FKLD')
        plt.plot(data_nc['acs_wv'],np.median(ap[mask=='ANTA',:],axis=0),label='ANTA')
        #plt.legend(loc=(1.04, 0))
        plt.xlim(400,720)
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$a_{p}$ [m$^{-1}$]')
        ax = plt.gca()
        plt.text(.86, .95,  'A', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
                
        
        plt.subplot(2,3,2)
        plt.plot(data_nc['acs_wv'],np.median(bp[mask=='NECS',:],axis=0),label='NECS')
        plt.plot(data_nc['acs_wv'],np.median(bp[mask=='NADR',:],axis=0),label='NADR')
        plt.plot(data_nc['acs_wv'],np.median(bp[mask=='NASE',:],axis=0),label='NASE')
        plt.plot(data_nc['acs_wv'],np.median(bp[mask=='NATR',:],axis=0),label='NATR')
        plt.plot(data_nc['acs_wv'],np.median(bp[mask=='WTRA',:],axis=0),label='WTRA')
        plt.plot(data_nc['acs_wv'],np.median(bp[mask=='SATL',:],axis=0),label='SATL')
        plt.plot(data_nc['acs_wv'],np.median(bp[mask=='SANT',:],axis=0),label='SANT')
        plt.plot(data_nc['acs_wv'],np.median(bp[mask=='SSTC',:],axis=0),label='SSTC')
        plt.plot(data_nc['acs_wv'],np.median(bp[mask=='FKLD',:],axis=0),label='FKLD')
        plt.plot(data_nc['acs_wv'],np.median(bp[mask=='ANTA',:],axis=0),label='ANTA')
        #plt.legend(loc=(1.04, 0))
        plt.xlim(400,720)
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$b_{p}$ [m$^{-1}$]')
        ax = plt.gca()
        plt.text(.86, .95,   'B', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
        
         
        plt.subplot(2,3,3)   
        plt.plot(data_nc['acs_wv'],np.median(cp[mask=='NECS',:],axis=0),label='NECS')
        plt.plot(data_nc['acs_wv'],np.median(cp[mask=='NADR',:],axis=0),label='NADR')
        plt.plot(data_nc['acs_wv'],np.median(cp[mask=='NASE',:],axis=0),label='NASE')
        plt.plot(data_nc['acs_wv'],np.median(cp[mask=='NATR',:],axis=0),label='NATR')
        plt.plot(data_nc['acs_wv'],np.median(cp[mask=='WTRA',:],axis=0),label='WTRA')
        plt.plot(data_nc['acs_wv'],np.median(cp[mask=='SATL',:],axis=0),label='SATL')
        plt.plot(data_nc['acs_wv'],np.median(cp[mask=='SANT',:],axis=0),label='SANT')
        plt.plot(data_nc['acs_wv'],np.median(cp[mask=='SSTC',:],axis=0),label='SSTC')
        plt.plot(data_nc['acs_wv'],np.median(cp[mask=='FKLD',:],axis=0),label='FKLD')
        plt.plot(data_nc['acs_wv'],np.median(cp[mask=='ANTA',:],axis=0),label='ANTA')
        plt.legend(loc=(1.04, 0))
        plt.xlim(400,720)
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$c_{p}$ [m$^{-1}$]')
        ax = plt.gca()
        plt.text(.86, .95,   'C', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
        
        
        plt.subplot(2,3,4)
        plt.plot(data_nc['acs_wv'],np.median(ap_int[mask=='NECS',:],axis=0),label='NECS')
        plt.plot(data_nc['acs_wv'],np.median(ap_int[mask=='NADR',:],axis=0),label='NADR')
        plt.plot(data_nc['acs_wv'],np.median(ap_int[mask=='NASE',:],axis=0),label='NASE')
        plt.plot(data_nc['acs_wv'],np.median(ap_int[mask=='NATR',:],axis=0),label='NATR')
        plt.plot(data_nc['acs_wv'],np.median(ap_int[mask=='WTRA',:],axis=0),label='WTRA')
        plt.plot(data_nc['acs_wv'],np.median(ap_int[mask=='SATL',:],axis=0),label='SATL')
        plt.plot(data_nc['acs_wv'],np.median(ap_int[mask=='SANT',:],axis=0),label='SANT')
        plt.plot(data_nc['acs_wv'],np.median(ap_int[mask=='SSTC',:],axis=0),label='SSTC')
        plt.plot(data_nc['acs_wv'],np.median(ap_int[mask=='FKLD',:],axis=0),label='FKLD')
        plt.plot(data_nc['acs_wv'],np.median(ap_int[mask=='ANTA',:],axis=0),label='ANTA')
        #plt.legend(loc=(1.04, 0))
        plt.xlim(400,720)
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$<a_{p}>$')
        ax = plt.gca()
        plt.text(.86, .95, 'D', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
        
        
        plt.subplot(2,3,5)
        plt.plot(data_nc['acs_wv'],np.median(bp_int[mask=='NECS',:],axis=0),label='NECS')
        plt.plot(data_nc['acs_wv'],np.median(bp_int[mask=='NADR',:],axis=0),label='NADR')
        plt.plot(data_nc['acs_wv'],np.median(bp_int[mask=='NASE',:],axis=0),label='NASE')
        plt.plot(data_nc['acs_wv'],np.median(bp_int[mask=='NATR',:],axis=0),label='NATR')
        plt.plot(data_nc['acs_wv'],np.median(bp_int[mask=='WTRA',:],axis=0),label='WTRA')
        plt.plot(data_nc['acs_wv'],np.median(bp_int[mask=='SATL',:],axis=0),label='SATL')
        plt.plot(data_nc['acs_wv'],np.median(bp_int[mask=='SANT',:],axis=0),label='SANT')
        plt.plot(data_nc['acs_wv'],np.median(bp_int[mask=='SSTC',:],axis=0),label='SSTC')
        plt.plot(data_nc['acs_wv'],np.median(bp_int[mask=='FKLD',:],axis=0),label='FKLD')
        plt.plot(data_nc['acs_wv'],np.median(bp_int[mask=='ANTA',:],axis=0),label='ANTA')
        #plt.legend(loc=(1.04, 0))
        plt.xlim(400,720)
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$<b_{p}>$')
        ax = plt.gca()
        plt.text(.86, .95,  'E', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
        
         
        plt.subplot(2,3,6)   
        plt.plot(data_nc['acs_wv'],np.median(cp_int[mask=='NECS',:],axis=0),label='NECS')
        plt.plot(data_nc['acs_wv'],np.median(cp_int[mask=='NADR',:],axis=0),label='NADR')
        plt.plot(data_nc['acs_wv'],np.median(cp_int[mask=='NASE',:],axis=0),label='NASE')
        plt.plot(data_nc['acs_wv'],np.median(cp_int[mask=='NATR',:],axis=0),label='NATR')
        plt.plot(data_nc['acs_wv'],np.median(cp_int[mask=='WTRA',:],axis=0),label='WTRA')
        plt.plot(data_nc['acs_wv'],np.median(cp_int[mask=='SATL',:],axis=0),label='SATL')
        plt.plot(data_nc['acs_wv'],np.median(cp_int[mask=='SANT',:],axis=0),label='SANT')
        plt.plot(data_nc['acs_wv'],np.median(cp_int[mask=='SSTC',:],axis=0),label='SSTC')
        plt.plot(data_nc['acs_wv'],np.median(cp_int[mask=='FKLD',:],axis=0),label='FKLD')
        plt.plot(data_nc['acs_wv'],np.median(cp_int[mask=='ANTA',:],axis=0),label='ANTA')
        #plt.legend(loc=(1.04, 0))
        plt.xlim(400,720)
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$<c_{p}>$')
        ax = plt.gca()
        plt.text(.86, .95,   'F', ha='left', va='top', transform=ax.transAxes,fontsize=24)


        plt.tight_layout()

        return


# legacy

#def plot_median_ap_cp(data_nc):
        
 #   plt.figure()
  #  plt.rcParams.update({'font.size': 16})
   # plt.plot(data_nc['acs_wv'],np.nanmedian(data_nc['acs_ap'],axis=0),label='$a_{p}$')
   # plt.plot(data_nc['acs_wv'],np.nanmedian(data_nc['acs_cp'],axis=0),label='$c_{p}$')
   # plt.xlim(400,750)
   # plt.xlabel('Wavelength [nm]')
   # plt.ylabel('[m$^{-1}$]')
    #plt.legend()
    
   # return

#def hplc_ap_match_up(data_nc, two_systems=False):
    
 #   ' extracts ap spectra that match hplc timestamps - median average over 30 min bin.'''
  #  ' this is designed to replicate what was done in the tot-chl-a match up procedure'

   # tol = np.timedelta64(60*30,'s') # +/- 30 minutes of timestamp - window size
    #
    #if two_systems == False:
        
      #  ap = data_nc['acs_ap'].data 
     #   ap_matches = []
        
       # for i in range(len(data_nc['hplc_time'])):
        #    delta_t = abs(data_nc['time'].data  - data_nc['hplc_time'].data[i])
         #   ap_i = ap[delta_t < tol,:] 
          #  if len(ap_i) > 1: 
           #     print(len(ap_i))
            #    ap_matches.append(np.nanmedian(ap_i,axis=0))
           # else:
            #    print(0)
             #   ap_matches.append(np.nan*np.ones(176))
            
        #ap_matches = np.asarray(ap_matches)
        #ap_matches_2D = ap_matches.reshape(len(data_nc['hplc_time']),len(data_nc['acs_wv']))
    
    
    #elif two_systems == True:
              
       # ap1 = data_nc['acs_ap'].data 
        #ap2 = data_nc['acs2_ap'].data 
        
       # ap_matches = []
       # for i in range(len(data_nc['hplc_time'])):
         #   delta_t = abs(data_nc['time'].data  - data_nc['hplc_time'].data[i])
         #   ap1_i = ap1[delta_t < tol,:] 
         #   ap2_i = ap2[delta_t < tol,:] 
         #   if len(ap1_i) > 1: 
        #        ap_matches.append(np.nanmedian(ap1_i,axis=0))
       #     elif len(ap2_i) > 1: 
      #          ap_matches.append(np.nanmedian(a2_i,axis=0))
      #      else:
     #           ap_matches.append(np.nan*np.ones(176))
     #          
    #    ap_matches = np.asarray(ap_matches)
    #    ap_matches_2D = ap_matches.reshape(len(data_nc['hplc_time']),len(data_nc['acs2_wv']))
            
   # return ap_matches_2D

    
 #  pigment plots  
#
def _filter_combine_pigs(data_nc):
    
    'returns a filtered pigment data frame with deeper samples and replicates removed'
    
    df_hplc = pd.DataFrame(index = data_nc.hplc_time)
    df_hplc['hplc_lat'] = data_nc['hplc_lat']
    df_hplc['hplc_lon'] = data_nc['hplc_lon']
    
    pig_keys = ['hplc_Tot_Chl_a','hplc_Tot_Chl_b', 'hplc_Tot_Chl_c' , 
                'hplc_Allo','hplc_Alpha-beta-Car', 'hplc_But-fuco', 
                'hplc_Diadino','hplc_Diato','hplc_Fuco','hplc_Hex-fuco',
                'hplc_Perid', 'hplc_Zea']
       
    for i in range(len(pig_keys)):
        df_hplc[pig_keys[i]] = data_nc[pig_keys[i]]
    
    df_hplc = df_hplc[data_nc['hplc_depth'].values < 10]  # no deep 
    print(len(df_hplc))
    df_hplc = df_hplc.groupby(df_hplc.index).mean()    # 1 reps
    print(len(df_hplc))

    return df_hplc


def _pig_9cruises():

    df_hplc_29 = _filter_combine_pigs(data_nc_29)
    df_hplc_28 = _filter_combine_pigs(data_nc_28)
    df_hplc_27 = _filter_combine_pigs(data_nc_27)
    df_hplc_26 = _filter_combine_pigs(data_nc_26)
    df_hplc_25 = _filter_combine_pigs(data_nc_25)
    df_hplc_24 = _filter_combine_pigs(data_nc_24)
    df_hplc_23 = _filter_combine_pigs(data_nc_23)
    df_hplc_22 = _filter_combine_pigs(data_nc_22)
    df_hplc_19 = _filter_combine_pigs(data_nc_19)
    
    df_hplc_combined = pd.concat([df_hplc_29, df_hplc_28, df_hplc_27, df_hplc_26, df_hplc_25, df_hplc_24, df_hplc_23, df_hplc_22, df_hplc_19 ])
    
   # df_hplc_combined = pd.concat([df_hplc_29, df_hplc_28, df_hplc_27, df_hplc_22, df_hplc_19 ]) - removing PML
    
    
    return    df_hplc_combined 


def pig_cov(data):
    
    pig_keys = ['hplc_Tot_Chl_a','hplc_Tot_Chl_b', 'hplc_Tot_Chl_c' , 
                'hplc_Allo','hplc_Alpha-beta-Car', 'hplc_But-fuco', 
                'hplc_Diadino','hplc_Diato','hplc_Fuco',
                'hplc_Hex-fuco','hplc_Perid', 'hplc_Zea']
    
    pig_labels = ['Tot_Chl_a','Tot_Chl_b', 'Tot_Chl_c' , 
                  'Allo','$\alpha$-$\beta$-Car', 'But-fuco', 
                  'Diadino','Diato','Hex-fuco',
                  'Perid', 'Zea']
    
    plt.figure(figsize=(15,20))
    plt.subplot(2,1,1)
    plt.rcParams.update({'font.size': 22})
    C = np.nan*np.ones([len(pig_keys), len(pig_keys)])
    for i in range(len(pig_keys)):
        for j in range(len(pig_keys)):
            if i > j:
                C[i, j] = scipy.stats.pearsonr(data[pig_keys[i]], data[pig_keys[j]])[0]
            if i < j: 
                C[i, j] = scipy.stats.pearsonr(data[pig_keys[i]]/data['hplc_Tot_Chl_a'], data[pig_keys[j]]/data['hplc_Tot_Chl_a'])[0]
    
    
    plt.pcolor(np.flipud(C.T), cmap='bwr', vmin=-1, vmax=1)
    cbar = plt.colorbar()
    cbar.set_label('Correlation coefficient, $r$', rotation=90 ,labelpad=30)
    
    plt.xlim(1,12)
    plt.ylim(0,12)
    
    for i in range(C.shape[0]):
        for j in range(C.shape[1]):
            if np.isnan(C[i, j]) == 0:
                plt.text(i + 0.5, 11-j + 0.5, '%.2f' % C[i, j],
                         horizontalalignment='center',
                         verticalalignment='center',
                         )
    ax = plt.gca()
    ax.set_xticks([1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5], 
                  labels =  ['Tot_Chl_b', 'Tot_Chl_c' , 
                                'Allo','alpha-beta-Car', 'But-fuco', 
                                'Diadino','Diato','Fuco',
                                'Hex-fuco','Perid', 'Zea'])
           
    plt.xticks(rotation=45)                 
    ax.set_yticks([0.5, 1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5],
           labels = ['Zea','Perid', 'Hex-fuco', 'Fuco', 'Diato',  'Diadino', 'But-fuco', 'alpha-beta-Car', 'Allo',  'Tot_Chl_c' ,'Tot_Chl_b', 'Tot_Chl_a'])
          
    
    plt.subplot(2,1,2)
    plt.rcParams.update({'font.size': 22})
    C = np.nan*np.ones([len(pig_keys), len(pig_keys)])
    for i in range(len(pig_keys)):
        for j in range(len(pig_keys)):
            if i > j:
                C[i, j] = (scipy.stats.pearsonr(data[pig_keys[i]], data[pig_keys[j]])[0])**2
            if i < j: 
               C[i, j] = (scipy.stats.pearsonr(data[pig_keys[i]]/data['hplc_Tot_Chl_a'], data[pig_keys[j]]/data['hplc_Tot_Chl_a'])[0])**2
    
    plt.pcolor(np.flipud(C.T), cmap='Oranges', vmin=0, vmax=1)
    cbar = plt.colorbar()
    cbar.set_label('Coefficient of determination, $r^2$', rotation=90, labelpad=30)
    
    plt.xlim(1,12)
    plt.ylim(0,12)
    
    for i in range(C.shape[0]):
        for j in range(C.shape[1]):
            if np.isnan(C[i, j]) == 0:
                plt.text(i + 0.5, 11-j + 0.5, '%.2f' % C[i, j],
                         horizontalalignment='center',
                         verticalalignment='center',
                         )
    ax = plt.gca()
    ax.set_xticks([1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5], 
                    labels =  ['Tot_Chl_b', 'Tot_Chl_c' , 
                                  'Allo','alpha-beta-Car', 'But-fuco', 
                                  'Diadino','Diato','Fuco',
                                  'Hex-fuco','Perid', 'Zea'])
             
    plt.xticks(rotation=45)                 
    ax.set_yticks([0.5, 1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5],
    labels = ['Zea','Perid', 'Hex-fuco', 'Fuco', 'Diato',  'Diadino', 'But-fuco', 'alpha-beta-Car', 'Allo',  'Tot_Chl_c' ,'Tot_Chl_b', 'Tot_Chl_a'])
    plt.tight_layout(pad=1.2)
    
    return


def pigs_byprov(): 
       
    prov= ['NADR', 'NASE', 'NASW', 'NATR', 'WTRA', 'SATL', 'SANT', 'SSTC'] # FKLD==1, ANTA -4, NECS = 8. do not include
                               
    pig_keys = ['hplc_Tot_Chl_a','hplc_Tot_Chl_b', 'hplc_Tot_Chl_c' , 
                'hplc_Allo','hplc_Alpha-beta-Car', 'hplc_But-fuco', 
                'hplc_Diadino','hplc_Diato','hplc_Fuco',
                'hplc_Hex-fuco','hplc_Perid', 'hplc_Zea']
     
    labels =  ['Tot_Chl_a', 'Tot_Chl_b', 'Tot_Chl_c' , 
                  'Allo','alpha-beta-Car', 'But-fuco', 
                  'Diadino','Diato','Fuco',
                  'Hex-fuco','Perid', 'Zea']
    
    letters = ['A', 'B', 'C' ,'D' ,'E' ,'F','G' , 'H']
    
    
    plt.figure(figsize=(13,13))    
    for i in range(len(prov)):
        plt.subplot(3, 3, i+1)  
        plt.title(prov[i] + ': N = ' + str(np.sum(LH_HPLC['LH_Province'] ==str(prov[i]))))
    
        pig_vec = []
        for j in range(len(pig_keys)):
            pig_j = np.array(df_hplc[pig_keys[j]])[LH_HPLC['LH_Province'] ==str(prov[i])]
            pig_j[pig_j==0] = np.nan
            pig_vec.append(pig_j)                                                                   
    
        bp = plt.boxplot(pig_vec ,showfliers=True,patch_artist=True, medianprops=dict(color='black'), whis=[10,90]) 
        plt.yscale('log')
    
        
        ax = plt.gca()
        ax.set_xticks([])  
        ax.set_ylim([0.0005, 5])
        if i==0:
            plt.ylabel('Concentration mg m$^{-3}$')
        if i==3:
            plt.ylabel('Concentration mg m$^{-3}$')
        if i==5:
            plt.ylabel('Concentration mg m$^{-3}$')
            
    
        
        plt.text(.90, .95,  letters[i], ha='left', va='top', transform=ax.transAxes,fontsize=22)  
    
    
        colors = cm.Paired(np.linspace(0,1,len(pig_keys))) 
        patches = []
        for k in range(len(pig_keys)):
            bp['boxes'][k].set_facecolor(colors[k])
            color_string = cl.rgb2hex(colors[k])
            patches.append(mpatches.Patch(color=color_string, label=labels[k]))
        
            plt.tight_layout(pad=1.2)
       
    plt.subplot(3,3,9)
    plt.legend(handles=patches,fontsize=14,loc=2)
    ax = plt.gca()
    ax.set_xticks([])  
    ax.set_yticks([])  
    ax.set_axis_off()
    
    return


# ap match-up function. Note: desirble to do th 
def hplc_ap_match_up_V2(data_nc, pigment, wavelength, two_acs_systems=False):
    
    ' extracts ap spectra that match hplc timestamps - median average over 30 min bin.'
    ' this is designed to replicate what was done in the tot-chl-a match up procedure'
    ' important to first remove: (i) replicates, (ii) deeper waters'
    ' user specifies wavelength and pigment'
    
    #  wl = data_nc_28['acs_wv'].values #        
    wlindex = int(round((wavelength - 400)/2))

    if two_acs_systems == False:  
        
        # hplc data frame
        df_hplc = pd.DataFrame({pigment: data_nc[pigment].values}, index = data_nc.hplc_time)
        df_hplc = df_hplc[data_nc['hplc_depth'].values < 10] # no deep 
        df_hplc = df_hplc.groupby(df_hplc.index).mean() # 1 reps
             
        # ap data frame
        ap = data_nc['acs_ap'].data 
        ap_f = sg.medfilt(ap[:,wlindex],kernel_size=31) # med filt in time dimension; window size = 30 min        
        df_ap = pd.Series(ap_f, index=data_nc['time'])
        
        #
        df_hplc_ap = pd.DataFrame({pigment: df_hplc[pigment], 'ap_' +str(wavelength): df_ap})
        df_hplc_ap = df_hplc_ap.interpolate('index',limit=1).reindex(index=df_hplc.index,method='nearest',tolerance='30min')

    if str(data_nc['time'].values[0])[0:4] == '2013': # trim AMT23
        df_hplc_ap = df_hplc_ap.iloc[0:25] # clip off nans for back-half of cruise (no data)
      
    if str(data_nc['time'].values[0])[0:4] == '2012': # trim AMT22
        df_hplc_ap.iloc[192,1] = np.nan # for consistency with  Tchl- match-up
      
    elif two_acs_systems == True:  
        
        df_hplc = pd.DataFrame({pigment: data_nc[pigment].values}, index = data_nc.hplc_time)
        df_hplc =  df_hplc[data_nc['hplc_depth'].values < 10] # no deep 
        df_hplc = df_hplc.groupby(df_hplc.index).mean() # 1 reps
             
        # ap data frame
        ap1 = data_nc['acs_ap']
        ap_f1 = sg.medfilt(ap1[:,wlindex],kernel_size=31) # med filt in time dimension; window size = 30 min        
        
        # ap2 data frame
        ap2 = data_nc['acs2_ap'].data 
        ap_f2 = sg.medfilt(ap2[:,wlindex],kernel_size=31) # med filt in time dimension; window size = 30 min        
            
        # matrix
        mat_ap = np.zeros([2,len(ap1)])
        mat_ap[0,:] = ap_f1
        mat_ap[1,:] = ap_f2 
        combined_ap = np.nanmean(mat_ap, axis=0)
               
        #
        df_ap = pd.Series(combined_ap,  index=data_nc['time'])            
        df_hplc_ap = pd.DataFrame({pigment: df_hplc[pigment], 'ap_' +str(wavelength): df_ap})
        df_hplc_ap = df_hplc_ap.interpolate('index',limit=1).reindex(df_hplc.index,method='nearest',tolerance='30min')
                                      
    return  df_hplc_ap 


# legacy
#def med_unc_plot(data_nc):
   
 #   plt.figure()
 #   plt.title('Median percentage uncertanity at each wl')
  #  plt.rcParams.update({'font.size': 16})
 #   plt.plot(data_nc['acs_wv'],np.nanmedian(100*abs(data_nc['acs_ap_u']/data_nc['acs_ap']),axis=0),label='$a_{p}$ unc')
 #   plt.plot(data_nc['acs_wv'],np.nanmedian(100*abs(data_nc['acs_bp_u']/data_nc['acs_bp']),axis=0),label='$b_{p}$ unc')
 #   plt.plot(data_nc['acs_wv'],np.nanmedian(100*abs(data_nc['acs_cp_u']/data_nc['acs_cp']),axis=0),label='$c_{p}$ unc')
 #   plt.ylim(0,80)
  #  plt.xlim(400,750)
   # plt.xlabel('Wavelength [nm]')
   # plt.ylabel('Percentage uncertainty [%]')
   # plt.legend()
    
        #    
  #  plt.figure()
   # plt.title('Median uncertanity at each wl')
   # plt.rcParams.update({'font.size': 16})
   # plt.plot(data_nc['acs_wv'],np.median(data_nc['acs_ap_u'],axis=0),label='$a_{p}$ unc')
   # plt.plot(data_nc['acs_wv'],np.median(data_nc['acs_bp_u'],axis=0),label='$b_{p}$ unc')
   # plt.plot(data_nc['acs_wv'],np.median(data_nc['acs_cp_u'],axis=0),label='$c_{p}$ unc')
   # plt.xlim(400,750)
   # plt.xlabel('Wavelength [nm]')
   # plt.ylabel('Uncertainty [m$^{-1}$]')
   # plt.legend()
    
 #   plt.figure()
   # plt.title('Median spectra')
  #  plt.rcParams.update({'font.size': 16})
  #  plt.plot(data_nc['acs_wv'],np.median(data_nc['acs_ap'],axis=0),label='$a_{p}$')
  #  plt.plot(data_nc['acs_wv'],np.median(data_nc['acs_ap'],axis=0),label='$b_{p}$')
   # plt.plot(data_nc['acs_wv'],np.median(data_nc['acs_cp'],axis=0),label='$c_{p}$')
  #  plt.xlim(400,750)
   # plt.xlabel('Wavelength [nm]')
   # plt.ylabel('Uncertainty [m$^{-1}$]')
    #plt.legend()

  #  return


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


def conc_plot():

    plt.figure(figsize=(16,16))
    plt.suptitle('Median spectra by concentration range')
     
    ap_22 = data_nc_22['acs_ap'].data
    ap_23 = data_nc_23['acs_ap'].data
    ap_24 = data_nc_24['acs_ap'].data
    ap_24_2 = data_nc_24['acs2_ap'].data
    ap_25 = data_nc_25['acs_ap'].data
    ap_25_2 = data_nc_25['acs2_ap'].data
    ap_26 = data_nc_26['acs_ap'].data
    ap_27 = data_nc_27['acs_ap'].data
    ap_28 = data_nc_28['acs_ap'].data
   # ap_dy = data_nc_dy151['acs_ap'].data

    
    chl_22 = data_nc_22['acs_chl_debiased'].data
    chl_23 = data_nc_23['acs_chl_debiased'].data
    chl_24 = data_nc_24['acs_chl_debiased'].data
    chl_24_2 = data_nc_24['acs2_chl_debiased'].data  
    chl_25 = data_nc_25['acs_chl_debiased'].values         
    chl_25_2 = data_nc_25['acs2_chl_debiased'].values           
    chl_26 = data_nc_26['acs_chl_debiased'].values
    chl_27 = data_nc_27['acs_chl_debiased'].values
    chl_28 = data_nc_28['acx_chl_debiased'].values
   # chl_dy =  data_nc_dy151['acs_chl'].values
    

    plt.subplot(2,2,1)
    ap_22_f1 = ap_22[(chl_22<0.5) & (chl_22>0.49)]
    ap_23_f1 = ap_23[(chl_23<0.5) & (chl_23>0.49)]
    ap_24_f1 = ap_24[(chl_24<0.5) & (chl_24>0.49)]
    ap_24_2_f1 = ap_24_2[(chl_24_2<0.5) & (chl_24_2>0.49)]
    ap_25_f1 = ap_25[(chl_25<0.5) & (chl_25>0.49)]
    ap_25_2_f1 = ap_25_2[(chl_25_2<0.5) & (chl_25_2>0.49)]
    ap_26_f1 = ap_26[(chl_26<0.5) & (chl_26>0.49)]
    ap_27_f1 = ap_27[(chl_27<0.5) & (chl_27>0.49)] 
    ap_28_f1 = ap_28[(chl_28<0.5) & (chl_28>0.49)]
 #   ap_dy_f1 = ap_dy[(chl_dy<0.5) & (chl_dy>0.49)]
    
    
    plt.title('0.49 < ACS Chl < 0.50 mg m$^{-3}$')
    plt.plot(data_nc_22['acs_wv'],np.nanmedian(ap_22_f1,axis=0),label='AMT 22: n =' + str(len(ap_22_f1)))
    plt.plot(data_nc_23['acs_wv'],np.nanmedian(ap_23_f1,axis=0),label='AMT 23: n =' + str(len(ap_23_f1)))
    plt.plot(data_nc_24['acs_wv'],np.nanmedian(ap_24_f1,axis=0),label='AMT 24: ACS: 1: n =' + str(len(ap_24_f1)))
    plt.plot(data_nc_24['acs2_wv'],np.nanmedian(ap_24_2_f1,axis=0),label='AMT 24: ACS 2: n =' + str(len(ap_24_2_f1)))
    plt.plot(data_nc_25['acs_wv'],np.nanmedian(ap_25_f1,axis=0),label='AMT 25: ACS: 1: n =' + str(len(ap_25_f1)))
    plt.plot(data_nc_25['acs2_wv'],np.nanmedian(ap_25_2_f1,axis=0),label='AMT 25: ACS 2: n =' + str(len(ap_25_2_f1)))
    plt.plot(data_nc_26['acs_wv'],np.nanmedian(ap_26_f1,axis=0),label='AMT 26: n =' + str(len(ap_26_f1)))
    plt.plot(data_nc_27['acs_wv'],np.nanmedian(ap_27_f1,axis=0),label='AMT 27: n =' + str(len(ap_27_f1)))
    plt.plot(data_nc_28['acs_wv'],np.nanmedian(ap_28_f1,axis=0),label='AMT 28: n =' + str(len(ap_28_f1)))
 #   plt.plot(data_nc_dy151['acs_wv'],np.nanmedian(ap_dy_f1,axis=0),label='DY151: n =' + str(len(ap_dy_f1)))
    plt.legend()
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('$a_{p}$ [m$^{-1}$]')
    

    plt.subplot(2,2,2)
    plt.title('0.09 < ACS Chl < 0.1 mg m$^{-3}$')   
    ap_22_f1 = ap_22[(chl_22<0.1) & (chl_22>0.09)]
    ap_23_f1 = ap_23[(chl_23<0.1) & (chl_23>0.09)]
    ap_24_f1 = ap_24[(chl_24<0.1) & (chl_24>0.09)]
    ap_24_2_f1 = ap_24_2[(chl_24_2<0.1) & (chl_24_2>0.09)]
    ap_25_f1 = ap_25[(chl_25<0.1) & (chl_25>0.09)]
    ap_25_2_f1 = ap_25_2[(chl_25_2<0.1) & (chl_25_2>0.09)]
    ap_26_f1 = ap_26[(chl_26<0.1) & (chl_26 > 0.09)]
    ap_27_f1 = ap_27[(chl_27<0.1) & (chl_27 > 0.09)] 
    ap_28_f1 = ap_28[(chl_28<0.1) & (chl_28 > 0.09)] 
 #   ap_dy_f1 = ap_dy[(chl_dy<0.1) & (chl_dy >0.09)]
    
    
    plt.plot(data_nc_22['acs_wv'],np.nanmedian(ap_22_f1,axis=0),label='AMT 22: n =' + str(len(ap_22_f1)))
    plt.plot(data_nc_23['acs_wv'],np.nanmedian(ap_23_f1,axis=0),label='AMT 23: n =' + str(len(ap_23_f1)))
    plt.plot(data_nc_24['acs_wv'],np.nanmedian(ap_24_f1,axis=0),label='AMT 24: ACS: 1: n =' + str(len(ap_24_f1)))
    plt.plot(data_nc_24['acs2_wv'],np.nanmedian(ap_24_2_f1,axis=0),label='AMT 24: ACS 2: n =' + str(len(ap_24_2_f1)))
    plt.plot(data_nc_25['acs_wv'],np.nanmedian(ap_25_f1,axis=0),label='AMT 25: ACS: 1: n =' + str(len(ap_25_f1)))
    plt.plot(data_nc_25['acs2_wv'],np.nanmedian(ap_25_2_f1,axis=0),label='AMT 25: ACS 2: n =' + str(len(ap_25_2_f1)))
    plt.plot(data_nc_26['acs_wv'],np.nanmedian(ap_26_f1,axis=0),label='AMT 26: n =' + str(len(ap_26_f1)))
    plt.plot(data_nc_27['acs_wv'],np.nanmedian(ap_27_f1,axis=0),label='AMT 27: n =' + str(len(ap_27_f1)))
    plt.plot(data_nc_28['acs_wv'],np.nanmedian(ap_28_f1,axis=0),label='AMT 28: n =' + str(len(ap_28_f1)))
#    plt.plot(data_nc_dy151['acs_wv'],np.nanmedian(ap_dy_f1,axis=0),label='DY151: n =' + str(len(ap_dy_f1)))
    plt.legend()
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('$a_{p}$ [m$^{-1}$]')

    
    plt.subplot(2,2,3)
    ap_22_f1 = ap_22[(chl_22<1) & (chl_22>0.8)]
    ap_23_f1 = ap_23[(chl_23<1) & (chl_23>0.8)]
    ap_24_f1 = ap_24[(chl_24<1) & (chl_24>0.8)]
    ap_24_2_f1 = ap_24_2[(chl_24_2<1) & (chl_24_2>0.8)]
    ap_25_f1 = ap_25[(chl_25<1) & (chl_25>0.8)]
    ap_25_2_f1 = ap_25_2[(chl_25_2<1) & (chl_25_2>0.8)]
    ap_26_f1 = ap_26[(chl_26<1) & (chl_26 > 0.8)]
    ap_27_f1 = ap_27[(chl_27 <1) & (chl_27 > 0.8)] 
    ap_28_f1 = ap_28[(chl_28<1) & (chl_28 > 0.8)] 
 #   ap_dy_f1 = ap_dy[(chl_dy<1) & (chl_dy >0.8)]
    
    
    plt.title('0.8 < ACS Chl < 1 mg m$^{-3}$')
    plt.plot(data_nc_22['acs_wv'],np.nanmedian(ap_22_f1,axis=0),label='AMT 22: n =' + str(len(ap_22_f1)))
    plt.plot(data_nc_23['acs_wv'],np.nanmedian(ap_23_f1,axis=0),label='AMT 23: n =' + str(len(ap_23_f1)))
    plt.plot(data_nc_24['acs_wv'],np.nanmedian(ap_24_f1,axis=0),label='AMT 24: ACS: 1: n =' + str(len(ap_24_f1)))
    plt.plot(data_nc_24['acs2_wv'],np.nanmedian(ap_24_2_f1,axis=0),label='AMT 24: ACS 2: n =' + str(len(ap_24_2_f1)))
    plt.plot(data_nc_25['acs_wv'],np.nanmedian(ap_25_f1,axis=0),label='AMT 25: ACS: 1: n =' + str(len(ap_25_f1)))
    plt.plot(data_nc_25['acs2_wv'],np.nanmedian(ap_25_2_f1,axis=0),label='AMT 25: ACS 2: n =' + str(len(ap_25_2_f1)))
    plt.plot(data_nc_26['acs_wv'],np.nanmedian(ap_26_f1,axis=0),label='AMT 26: n =' + str(len(ap_26_f1)))
    plt.plot(data_nc_27['acs_wv'],np.nanmedian(ap_27_f1,axis=0),label='AMT 27: n =' + str(len(ap_27_f1)))
    plt.plot(data_nc_28['acs_wv'],np.nanmedian(ap_28_f1,axis=0),label='AMT 28: n =' + str(len(ap_28_f1)))
   # plt.plot(data_nc_dy151['acs_wv'],np.nanmedian(ap_dy_f1,axis=0),label='DY151: n =' + str(len(ap_dy_f1)))
    plt.legend()
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('$a_{p}$ [m$^{-1}$]')

    
   # plt.subplot(2,2,4)
   # plt.title('ACS Chl < 0.01 mg m$^{-3}$')
   # ap_26_f1 = ap_26[chl_26<0.01]
   # ap_27_f1 = ap_27[chl_27<0.01]
   # ap_28_f1 = ap_28[chl_28<0.01]
   # plt.plot(data_nc_26['acs_wv'],np.nanmedian(ap_26_f1,axis=0),label='AMT 26: n =' + str(len(ap_26_f1)))
    # plt.plot(data_nc_27['acs_wv'],np.nanmedian(ap_27_f1,axis=0),label='AMT 27: n =' + str(len(ap_27_f1)))
     # plt.plot(data_nc_28['acs_wv'],np.nanmedian(ap_28_f1,axis=0),label='AMT 28: n =' + str(len(ap_28_f1)))
    # plt.legend()
    # plt.xlabel('Wavelength [nm]')
    # plt.ylabel('$a_{p}$ [m$^{-1}$]')
    plt.tight_layout()  
     
    
    return

    
def ap_scatter(pigment, pigment_label, wl_i):

    
    wl = data_nc_28['acs_wv'].values        
    index = int(round((wl_i-400)/2))
    
    # total vectors - use in the regression stats
    total_pig = np.concatenate((data_nc_22[pigment], data_nc_23[pigment], data_nc_24[pigment], data_nc_25[pigment], data_nc_26[pigment], data_nc_27[pigment], data_nc_28[pigment]))
    total_ap_m = np.concatenate((ap_m_22[:,index], ap_m_23[:,index], ap_m_24[:,index], ap_m_25[:,index], ap_m_26[:,index], ap_m_27[:,index],  ap_m_28[:,index]))
    
    total_pig = total_pig[total_ap_m >0] # must be > 0 for log fit   
    total_ap_m = total_ap_m[total_ap_m >0]
    
    total_pig_log = np.log(total_pig[~np.isnan(total_ap_m)])
    total_ap_log = np.log(total_ap_m[~np.isnan(total_ap_m)])
    
    linear_mod = scipy.stats.linregress(total_pig_log,  total_ap_log) 
    
    A = np.round(1000*linear_mod.intercept)/1000
    B = np.round(1000*linear_mod.slope)/1000
    r_sq = np.round(1000*linear_mod.rvalue**2)/1000
   
    rho = np.round(1000*scipy.stats.spearmanr(total_pig_log,  total_ap_log).correlation)/1000
    
    pred = (total_ap_m[~np.isnan(total_ap_m)]/np.exp(A))**(1/B)
    eps = np.round(1000*np.median((np.abs(pred - total_pig[~np.isnan(total_ap_m)])/total_pig[~np.isnan(total_ap_m)])))/10
                                            

    if pigment == 'hplc_Tot_Chl_a' or  pigment == 'hplc_PPC' or  pigment == 'hplc_PSC':
        X = np.arange(0.01,10,0.01)
    else:
        X = np.arange(0.001,1,0.01)
    Y = np.exp(A+B*np.log(X))
    
    
    plt.figure(figsize=(8,8))
    plt.title('Type 1: log-log fit: A = ' + str(A) + '; B = '  + str(B) + '\n' +
               'r$^{2}$ = ' + str(r_sq) + '; ' + 'spearmans = ' + str(rho) +  '; ' + 'error = ' + str(eps) + '$\%$')
    plt.rcParams.update({'font.size': 18})
    plt.plot(X,Y,'k',linestyle='dashed')
    plt.scatter(data_nc_19[pigment], ap_m_19[:,index],label='AMT19' ,alpha=0.6)
    plt.scatter(data_nc_22[pigment], ap_m_22[:,index],label='AMT22' ,alpha=0.6)
    plt.scatter(data_nc_23[pigment], ap_m_23[:,index],label='AMT23' ,alpha=0.6)
    plt.scatter(data_nc_24[pigment], ap_m_24[:,index],label='AMT24' ,alpha=0.6)   
    plt.scatter(data_nc_25[pigment], ap_m_25[:,index],label='AMT25' ,alpha=0.6)
    plt.scatter(data_nc_26[pigment], ap_m_26[:,index],label='AMT26' ,alpha=0.6)
    plt.scatter(data_nc_27[pigment], ap_m_27[:,index],label='AMT27' ,alpha=0.6)   
    plt.scatter(data_nc_28[pigment], ap_m_28[:,index],label='AMT28' ,alpha=0.6)  
    plt.plot(X,Y,'k',linestyle='dashed')
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel('AC-S: $a_{p}$(' + str(wl_i) +') [m$^{-1}]$')
    plt.xlabel('HPLC: ' + pigment_label + '  [mg m$^{-3}]$')
    plt.legend()

    return


def AB_combined(pigment, pigment_label, wl_i):

    
    wl = data_nc_28['acs_wv'].values        
    index = int(round((wl_i-400)/2))
    
    # total vectors - use in the regression stats
    total_pig = np.concatenate((data_nc_22[pigment], data_nc_23[pigment], data_nc_24[pigment], data_nc_25[pigment], data_nc_26[pigment], data_nc_27[pigment], data_nc_28[pigment]))
    total_ap_m = np.concatenate((ap_m_22[:,index], ap_m_23[:,index], ap_m_24[:,index], ap_m_25[:,index], ap_m_26[:,index], ap_m_27[:,index],  ap_m_28[:,index]))
    
    total_pig = total_pig[total_ap_m >0] # must be > 0 for log fit   
    total_ap_m = total_ap_m[total_ap_m >0]
    
    total_pig_log = np.log(total_pig[~np.isnan(total_ap_m)])
    total_ap_log = np.log(total_ap_m[~np.isnan(total_ap_m)])
    
    linear_mod = scipy.stats.linregress(total_pig_log,  total_ap_log) 
    
    A = np.round(1000*linear_mod.intercept)/1000
    B = np.round(1000*linear_mod.slope)/1000
    r_sq = np.round(1000*linear_mod.rvalue**2)/1000
   
    rho = np.round(1000*scipy.stats.spearmanr(total_pig_log,  total_ap_log).correlation)/1000
    
    pred = (total_ap_m[~np.isnan(total_ap_m)]/np.exp(A))**(1/B)
    eps = np.round(1000*np.median((np.abs(pred - total_pig[~np.isnan(total_ap_m)])/total_pig[~np.isnan(total_ap_m)])))/10
                                            

    if pigment == 'hplc_Tot_Chl_a' or  pigment == 'hplc_PPC' or  pigment == 'hplc_PSC':
        X = np.arange(0.01,10,0.01)
    else:
        X = np.arange(0.001,1,0.01)
    Y = np.exp(A + B*np.log(X))
    
    
    plt.figure(figsize=(8,8))
    plt.title('Type 1: log-log fit: A = ' + str(A) + '; B = '  + str(B) + '\n' +
               'r$^{2}$ = ' + str(r_sq) + '; ' + 'spearmans = ' + str(rho) +  '; ' + 'error = ' + str(eps) + '$\%$')
    plt.rcParams.update({'font.size': 18})
    plt.plot(X,Y,'k',linestyle='dashed')
    plt.scatter(data_nc_19[pigment], ap_m_19[:,index],label='AMT19' ,alpha=0.6)
    plt.scatter(data_nc_22[pigment], ap_m_22[:,index],label='AMT22' ,alpha=0.6)
    plt.scatter(data_nc_23[pigment], ap_m_23[:,index],label='AMT23' ,alpha=0.6)
    plt.scatter(data_nc_24[pigment], ap_m_24[:,index],label='AMT24' ,alpha=0.6)   
    plt.scatter(data_nc_25[pigment], ap_m_25[:,index],label='AMT25' ,alpha=0.6)
    plt.scatter(data_nc_26[pigment], ap_m_26[:,index],label='AMT26' ,alpha=0.6)
    plt.scatter(data_nc_27[pigment], ap_m_27[:,index],label='AMT27' ,alpha=0.6)   
    plt.scatter(data_nc_28[pigment], ap_m_28[:,index],label='AMT28' ,alpha=0.6)  
    plt.plot(X,Y,'k',linestyle='dashed')
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel('AC-S: $a_{p}$(' + str(wl_i) +') [m$^{-1}]$')
    plt.xlabel('HPLC: ' + pigment_label + '  [mg m$^{-3}]$')
    plt.legend()

    return

   
def pig_scatter(pigment, pigment_label, pigment_ref, pigment_ref_label):

    'scatter plots for inter-pigment relationships. Important to filter out: (i) deeper water, (ii) replicates prior to statistics'

    # create series object - remove deeper waters and replicates
    ds_29 = pd.DataFrame({pigment: data_nc_29[pigment].values,pigment_ref:  data_nc_29[pigment_ref].values}, index = data_nc_29.hplc_time)
    ds_29 = ds_29[data_nc_29['hplc_depth'].values < 10]
    ds_29 = ds_29.groupby(ds_29.index).mean()

    ds_28 = pd.DataFrame({pigment: data_nc_28[pigment].values, pigment_ref:  data_nc_28[pigment_ref].values}, index = data_nc_28.hplc_time)
    ds_28 = ds_28[data_nc_28['hplc_depth'].values < 10]
    ds_28 = ds_28.groupby(ds_28.index).mean() # no reps!
    
    ds_27 = pd.DataFrame({pigment: data_nc_27[pigment].values, pigment_ref:  data_nc_27[pigment_ref].values}, index = data_nc_27.hplc_time)
    ds_27 = ds_27[data_nc_27['hplc_depth'].values < 6] # no deep 
    ds_27 = ds_27.groupby(ds_27.index).mean() # no reps
    
    ds_26 = pd.DataFrame({pigment: data_nc_26[pigment].values, pigment_ref:  data_nc_26[pigment_ref].values}, index = data_nc_26.hplc_time)
    ds_26 = ds_26[data_nc_26['hplc_depth'].values < 6] # no deep 
    ds_26 = ds_26.groupby(ds_26.index).mean() # 1 reps
    
    ds_25 = pd.DataFrame({pigment: data_nc_25[pigment].values, pigment_ref:  data_nc_25[pigment_ref].values}, index = data_nc_25.hplc_time)
    ds_25 = ds_25[data_nc_25['hplc_depth'].values < 10] # no deep 
    ds_25 = ds_25.groupby(ds_25.index).mean() # 1 reps
    
    ds_24 = pd.DataFrame({pigment: data_nc_24[pigment].values, pigment_ref:  data_nc_24[pigment_ref].values}, index = data_nc_24.hplc_time)
    ds_24 = ds_24[data_nc_24['hplc_depth'].values < 10] # no deep 
    ds_24 = ds_24.groupby(ds_24.index).mean() # 1 reps
        
    ds_23 = pd.DataFrame({pigment: data_nc_23[pigment].values, pigment_ref:  data_nc_23[pigment_ref].values}, index = data_nc_23.hplc_time)
    ds_23 = ds_23[data_nc_23['hplc_depth'].values < 10] # no deep 
    ds_23 = ds_23.groupby(ds_23.index).mean() # 1 reps
     
    ds_22 = pd.DataFrame({pigment: data_nc_22[pigment].values, pigment_ref:  data_nc_22[pigment_ref].values}, index = data_nc_22.hplc_time)
    ds_22 = ds_22[data_nc_22['hplc_depth'].values < 10] # no deep 
    ds_22 = ds_22.groupby(ds_22.index).mean() # 1 reps
    
    ds_19 = pd.DataFrame({pigment: data_nc_19[pigment].values, pigment_ref:  data_nc_19[pigment_ref].values}, index = data_nc_19.hplc_time)
    ds_19 = ds_19[data_nc_19['hplc_Depth_(meters)'].values < 10] # no deep 
    ds_19 = ds_19.groupby(ds_19.index).mean() # 1 reps

    # total vectors - use in the regression stats
    total_pig = np.concatenate((ds_19[pigment], ds_22[pigment], ds_23[pigment], ds_24[pigment], ds_25[pigment], ds_26[pigment], ds_27[pigment], ds_28[pigment], ds_29[pigment]))
    total_pig_ref = np.concatenate((ds_19[pigment_ref], ds_22[pigment_ref], ds_23[pigment_ref], ds_24[pigment_ref], ds_25[pigment_ref], ds_26[pigment_ref], ds_27[pigment_ref], ds_28[pigment_ref], ds_29[pigment_ref]))


    total_pig_f =  total_pig[~np.isnan(total_pig)] 
    total_pig_ref_f = total_pig_ref[~np.isnan(total_pig)] 

    total_pig_f2 =  total_pig_f[~np.isinf(np.log(total_pig_f))] 
    total_pig_ref_f2 = total_pig_ref_f[~np.isinf(np.log(total_pig_f))] 
    
    linear_mod = scipy.stats.linregress(np.log(total_pig_f2), np.log(total_pig_ref_f2)) 
    
    A = np.round(1000*linear_mod.intercept)/1000
    B = np.round(1000*linear_mod.slope)/1000
    N = len(total_pig_f2)
    
    r_sq = np.round(1000*linear_mod.rvalue**2)/1000
    rho = np.round(1000*scipy.stats.spearmanr(np.log(total_pig_f2), np.log(total_pig_ref_f2)).correlation)/1000
    
    pred = (total_pig_ref_f2/np.exp(A))**(1/B)
    eps = np.round(1000*np.nanmedian((np.abs(pred - total_pig_f2)/total_pig_f2)))/10

    if pigment == 'hplc_Tot_Chl_a' or  pigment == 'hplc_PPC' or  pigment == 'hplc_PSC':
        X = np.arange(0.01,10,0.01)
    else:
        X = np.arange(0.001,1,0.01)
    Y = np.exp(A+B*np.log(X))
    
    
    plt.figure(figsize=(8,8))
    plt.title('Type 1: log-log fit: A = ' + str(A) + '; B = '  + str(B) +'; N = ' + str(N)  + '\n'
               'r$^{2}$ = ' + str(r_sq) + '; ' + 'spearmans = ' + str(rho) +  '; ' + 'error = ' + str(eps) + '$\%$')
    plt.rcParams.update({'font.size': 18})
    plt.plot(X,Y,'k',linestyle='dashed')
    plt.scatter(data_nc_19[pigment], data_nc_19[pigment_ref],label='AMT19' ,alpha=0.6)
    plt.scatter(data_nc_22[pigment], data_nc_22[pigment_ref],label='AMT22' ,alpha=0.6)
    plt.scatter(data_nc_23[pigment], data_nc_23[pigment_ref],label='AMT23' ,alpha=0.6)
    plt.scatter(data_nc_24[pigment], data_nc_24[pigment_ref],label='AMT24' ,alpha=0.6)   
    plt.scatter(data_nc_25[pigment], data_nc_25[pigment_ref],label='AMT25' ,alpha=0.6)
    plt.scatter(data_nc_26[pigment], data_nc_26[pigment_ref],label='AMT26' ,alpha=0.6)
    plt.scatter(data_nc_27[pigment], data_nc_27[pigment_ref],label='AMT27' ,alpha=0.6)   
    plt.scatter(data_nc_28[pigment], data_nc_28[pigment_ref],label='AMT28' ,alpha=0.6)  
    plt.scatter(data_nc_29[pigment], data_nc_29[pigment_ref],label='AMT29' ,alpha=0.6)  
    plt.plot(X,Y,'k',linestyle='dashed')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('HPLC: ' + pigment_label + '   [mg m$^{-3}]$')
    plt.ylabel('HPLC: ' + pigment_ref_label + '   [mg m$^{-3}]$')
    plt.legend()

    return


  
def total_ABfit():
        '''routine to do power-law fit for combined dataset'''
        
        # concatenate total chl
        tot_chl = np.concatenate([data_nc_19.acx_chl_debiased.HPLC_Tot_chla, 
        data_nc_22.acs_chl_debiased.HPLC_Tot_chla, 
        data_nc_23.acs_chl_debiased.HPLC_Tot_chla, 
        data_nc_24.acs_chl_debiased.HPLC_Tot_chla,
        data_nc_24.acs2_chl_debiased.HPLC_Tot_chla, 
        data_nc_25.acs_chl_debiased.HPLC_Tot_chla, 
       # data_nc_25.acs2_chl_debiased.HPLC_Tot_chla, 
        data_nc_26.acs_chl_debiased.HPLC_Tot_chla, 
        data_nc_27.acs_chl_debiased.HPLC_Tot_chla, 
        data_nc_28.acx_chl_debiased.HPLC_Tot_chla, 
        data_nc_29.acs_chl_debiased.HPLC_Tot_chla]) 
        
        # concatenate AC-S
        tot_acs= np.concatenate([data_nc_19.acx_chl_debiased.acx_chl, 
        data_nc_22.acs_chl_debiased.acs_chl , 
        data_nc_23.acs_chl_debiased.acs_chl , 
        data_nc_24.acs_chl_debiased.acs_chl , 
        data_nc_24.acs2_chl_debiased.acs_chl,
        data_nc_25.acs_chl_debiased.acs_chl ,
      #  data_nc_25.acs2_chl_debiased.acs_chl , 
        data_nc_26.acs_chl_debiased.acs_chl , 
        data_nc_27.acs_chl_debiased.acs_chl , 
        data_nc_28.acx_chl_debiased.acx_chl , 
        data_nc_29.acs_chl_debiased.acs_chl ]) 
            
        print(np.sum(~np.isnan(data_nc_19.acx_chl_debiased.HPLC_Tot_chla + data_nc_19.acx_chl_debiased.acx_chl)))
        print(np.sum(~np.isnan(data_nc_22.acs_chl_debiased.HPLC_Tot_chla + data_nc_22.acs_chl_debiased.acs_chl)))
        print(np.sum(~np.isnan(data_nc_23.acs_chl_debiased.HPLC_Tot_chla + data_nc_23.acs_chl_debiased.acs_chl)))
        print(np.sum(~np.isnan(data_nc_24.acs_chl_debiased.HPLC_Tot_chla + data_nc_24.acs_chl_debiased.acs_chl)))
        print(np.sum(~np.isnan(data_nc_24.acs2_chl_debiased.HPLC_Tot_chla + data_nc_24.acs2_chl_debiased.acs_chl)))
        print(np.sum(~np.isnan(data_nc_25.acs_chl_debiased.HPLC_Tot_chla + data_nc_25.acs_chl_debiased.acs_chl)))
        print(np.sum(~np.isnan(data_nc_26.acs_chl_debiased.HPLC_Tot_chla + data_nc_26.acs_chl_debiased.acs_chl)))
        print(np.sum(~np.isnan(data_nc_27.acs_chl_debiased.HPLC_Tot_chla + data_nc_27.acs_chl_debiased.acs_chl)))
        print(np.sum(~np.isnan(data_nc_28.acx_chl_debiased.HPLC_Tot_chla + data_nc_28.acx_chl_debiased.acx_chl)))                                          
        print(np.sum(~np.isnan(data_nc_29.acs_chl_debiased.HPLC_Tot_chla + data_nc_29.acs_chl_debiased.acs_chl)))                                           
                    
        
        tot_acs = tot_acs[~np.isnan(tot_chl)]
        tot_chl = tot_chl[~np.isnan(tot_chl)]
        
        tot_chl = tot_chl[~np.isnan(tot_acs)]
        tot_acs = tot_acs[~np.isnan(tot_acs)]
        
        
        linear_mod = scipy.stats.linregress(np.log10(0.014*tot_acs),np.log10(tot_chl))
                                                
        
        A = np.round(100*linear_mod.intercept)/100
        B = np.round(1000*linear_mod.slope)/1000
        r_sq = np.round(1000*linear_mod.rvalue**2)/100
        
        r_sq = np.round(1000*linear_mod.rvalue**2)/1000
        stderr = np.round(1000*linear_mod.stderr)/1000
        interr = np.round(1000*linear_mod.intercept_stderr)/1000
        
        
        print('A = ' + str(10**A) + ' +/- ' + str(2*10**interr))
        print('B = ' + str(B) + ' +/- ' + str(2*stderr))
        
        print(r_sq)
            
        
        rres = tot_acs/tot_chl - 1
        
        rres_log = np.log10(tot_acs) / np.log10(tot_chl) - 1
        
        delta = np.nanmedian(rres)
        delta_log = np.nanmedian(rres_log)
        sigma = prcrng(rres)
        sigma_log = prcrng(rres_log)
        N = np.sum(~np.isnan(rres))
        
        print('delta, sigma, N')
        print(delta, sigma, N)
        
        print('delta_log, sigma_log, N')
        print(delta_log, sigma_log, N)
        return
        

def prcrng(x):
    return (np.nanpercentile(x,84) - np.nanpercentile(x,16))/2



if __name__ == '__main__':
    
    ###########################################################################
    # AMT 19 - seabass and nc data (for cross-check)
    fn_nc_19 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT19/Processed/Step3/amt19_final_with_debiased_chl.nc'
    fn_chl_19 = '/data/datasets/cruise_data/active/ACS_Chl/amtacs/AMT19_ACS_CHL-A_MEDFILT_BIAS_CORRECTED.csv'
    data_chl_19 = pd.read_csv(fn_chl_19)
    data_nc_19 = xr.open_dataset(fn_nc_19)
    print(sorted(list(data_nc_19.keys())))

    
    ###########################################################################
    # AMT 22 - seabass and nc data (for cross-check)
    fn_nc_22 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT22/Processed/Step3/amt22_final_with_debiased_chl.nc'
    fn_chl_22 = '/data/datasets/cruise_data/active/ACS_Chl/amtacs/AMT22_ACS_CHL-A_MEDFILT_BIAS_CORRECTED.csv'
    data_chl_22 = pd.read_csv(fn_chl_22)
    data_nc_22 = xr.open_dataset(fn_nc_22)
    print(sorted(list(data_nc_22.keys())))
    
    
    ############################################################################
    # AMT 23 - seabass and nc data (for cross-check)
    fn_nc_23 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT23/Processed/underway/Step3/amt23_final_with_debiased_chl.nc'
    fn_chl_23 = '/data/datasets/cruise_data/active/ACS_Chl/amtacs/AMT23_ACS_CHL-A_MEDFILT_BIAS_CORRECTED.csv'
    data_nc_23 = xr.open_dataset(fn_nc_23)
    print(sorted(list(data_nc_23.keys()))) 
    
    
    ###########################################################################
    # AMT 24 - seabass and nc data (for cross-check)
    fn_nc_24 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT24/Processed/Uway/Step3/amt24_final_with_debiased_chl.nc'
    fn_chl_24 = '/data/datasets/cruise_data/active/ACS_Chl/amtacs/AMT25_ACS_CHL-A_MEDFILT_BIAS_CORRECTED.csv'
    data_nc_24 = xr.open_dataset(fn_nc_24)
    print(sorted(list(data_nc_24.keys()))) 
    
      
    ###########################################################################
    # AMT 25 - seabass and nc data (for cross-check)
    fn_nc_25 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT25/Processed/UWay/Step3/amt25_final_with_debiased_chl.nc'
    fn_chl_25 = '/data/datasets/cruise_data/active/ACS_Chl/amtacs/AMT25_ACS_CHL-A_MEDFILT_BIAS_CORRECTED.csv'
    data_nc_25 = xr.open_dataset(fn_nc_25)
    data_chl_25 = pd.read_csv(fn_chl_25)
    
    # output initial plots 
    #  plot_median_ap_cp(data_nc_25)
    #  med_unc_plot(data_nc_25)
    chl_check(data_nc_25,data_chl_25)
    

    ###########################################################################
    # AMT 26 - seabass and nc data (for cross-check)
    fn_mask_26 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/LH_mask_AMT26.csv'
    mask_26 = pd.read_csv(fn_mask_26)
    plot_median_ap_cp(data_nc_25)

    #pigment_check(data_nc_26, data_hplc_25)
    med_unc_plot(data_nc_25)
    chl_check(data_nc_25,data_chl_25)
    
    
          
    ###########################################################################
    # AMT 26 - seabass and nc data (for cross-check)
    dir_sb_26 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT26/Source/SeaBASS_submit/sb_processed/'
    fn_sb_acs_26 = dir_sb_26 + 'AMT26_InLine0_ACS_20160923_20161102_Particulate_v20230622.sb'
    fn_sb_hplc_26 = dir_sb_26 + 'AMT26_HPLC_20160923_20161102_v20230622.sb'
       
    fn_nc_26 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT26/Processed/Underway/Step3/amt26_final_with_debiased_chl.nc'
    fn_chl_26 = '/data/datasets/cruise_data/active/ACS_Chl/amtacs/AMT26_ACS_CHL-A_MEDFILT_BIAS_CORRECTED.csv'
       
    # load acs sb
    data_acs_26 = sbs.readSB(fn_sb_acs_26)
    print(data_acs_26.data.keys())
    data_acs_26 = data_acs_26.data
    
     # load hplc sb
    data_hplc_26 = sbs.readSB(fn_sb_hplc_26) 
    print(data_hplc_26.data.keys())
    data_hplc_26 = data_hplc_26.data
    
    # load netcdf
    data_nc_26 = xr.open_dataset(fn_nc_26)
    print(sorted(list(data_nc_26.keys()))) 
    
    # load previously-processed chl
    data_chl_26 = pd.read_csv(fn_chl_26)
    
    # load longhurst mask
    fn_mask_26 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/LH_mask_AMT26.csv'
    mask_26 = pd.read_csv(fn_mask_26)
        
    # output initial plots 
    # plot_median_ap_cp(data_nc_26)
    # plot_median_ap_province(data_nc_26, mask_26['LH_Province'])
    # pigment_check(data_nc_26, data_hplc_26)
    # med_unc_plot(data_nc_26)
    chl_check(data_nc_26,data_chl_26)
       

      
    ###########################################################################
    # AMT 27 - seabass and nc data (for cross-check)
    dir_sb_27 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT27/SeaBASS_submit/sb_processed/'
    fn_sb_acs_27 = dir_sb_27 + 'AMT27_InLine_ACS_20170924_20171101_Particulate_v20230622.sb'
    fn_sb_hplc_27 = dir_sb_27 + 'AMT27_HPLC_20170924_20171101_v20230622.sb'
    fn_nc_27 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT27/Processed/Underway/Step3/amt27_final_with_debiased_chl.nc'
    fn_chl_27 = '/data/datasets/cruise_data/active/ACS_Chl/amtacs/AMT27_ACS_CHL-A_MEDFILT_BIAS_CORRECTED.csv'
      
    
    # load acs sb
    data_acs_27 = sbs.readSB(fn_sb_acs_27)
    print(data_acs_27.data.keys())
    data_acs_27 = data_acs_27.data
    
     # load hplc sb
    data_hplc_27 = sbs.readSB(fn_sb_hplc_27) 
    print(data_hplc_27.data.keys())
    data_hplc_27 = data_hplc_27.data
    
    # load netcdf
    data_nc_27 = xr.open_dataset(fn_nc_27)
    print(sorted(list(data_nc_27.keys())))   
          
    # load longhurst mask
    fn_mask_27 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/LH_mask_AMT27.csv'
    mask_27 = pd.read_csv(fn_mask_27)
          
    # load previously-processed chl
    data_chl_27 = pd.read_csv(fn_chl_27)
        
    # output initial plots 
    # plot_median_ap()
    # plot_median_ap_cp(data_nc_27)

    # pigment_check(data_nc_27, data_hplc_27)
    # med_unc_plot(data_nc_27)
    # chl_check(data_nc_27, data_chl_27)
    
      
    ###########################################################################
    # AMT 28 - seabass and nc data (for cross-check)     
    dir_sb_28 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT28/Source/SeaBASS_submit/sb_processed/'
    fn_sb_acs_28 =  dir_sb_28 + 'AMT28_InLine0_ACS_20180925_20181028_Particulate_v20230622.sb'
    fn_sb_hplc_28 =  dir_sb_28 + 'AMT28_HPLC_20180925_20181027_v20230622.sb'
    fn_nc_28 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT28/Processed/Underway/Step3/amt28_final_with_debiased_chl.nc'
    fn_chl_28 = '/data/datasets/cruise_data/active/ACS_Chl/amtacs/AMT28_ACS_CHL-A_MEDFILT_BIAS_CORRECTED.csv'
    
    # load acs sb
    data_acs_28 = sbs.readSB(fn_sb_acs_28)
    print(data_acs_28.data.keys())
    data_acs_28 = data_acs_28.data
    
     # load hplc sb
    data_hplc_28 = sbs.readSB(fn_sb_hplc_28) 
    print(data_hplc_28.data.keys())
    data_hplc_28 = data_hplc_28.data
    
    # load netcdf
    data_nc_28 = xr.open_dataset(fn_nc_28)
    print(sorted(list(data_nc_28.keys())))   
     
    # load longhurst mask
    fn_mask_28 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/LH_mask_AMT28.csv'
    mask_28 = pd.read_csv(fn_mask_28)
               
    # load previously-processed chl
    #data_chl_28 = pd.read_csv(fn_chl_28)
    
    
    ###########################################################################
    # AMT 29 - seabass and nc data (for cross-check)     
    fn_nc_29 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT29/Processed/Step3/amt29_final_with_debiased_chl.nc'
    data_nc_29 = xr.open_dataset(fn_nc_29)
    data_nc_29.keys()


    #LH29 = longhurst_mask(data_nc_29['uway_lat'], data_nc_29['uway_lon'], 'AMT29')
    #LH27 = longhurst_mask(data_nc_27['uway_lat'], data_nc_27['uway_lon'], 'AMT27')
    LH25 = longhurst_mask(data_nc_25['uway_lat'], data_nc_25['uway_lon'], 'AMT25')
    LH24 = longhurst_mask(data_nc_24['uway_lat'], data_nc_24['uway_lon'], 'AMT24')
    #LH23 = longhurst_mask(data_nc_23['uway_lat'], data_nc_24['uway_lon'], 'AMT23')
    #LH22 = longhurst_mask(data_nc_22['uway_lat'], data_nc_24['uway_lon'], 'AMT22')
    #LH19 = longhurst_mask(data_nc_19['uway_lat'], data_nc_24['uway_lon'], 'AMT19')
     
    
    LH29 = longhurst_mask_HPLC(data_nc_29['hplc_lat'], data_nc_29['hplc_lon'], 'AMT29')

    # keys = list(data_nc.data_vars)
    # for i in range(len(keys)):
    #     if "hplc" in keys()[i]:
    
    chase = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/Chase_AB.csv'
    chase_data = pd.read_csv(chase) 
    
    A = np.array(chase_data['Achl (m-1)'])
    B = np.array(chase_data['Bchl (unitless)'])
    wl =  np.array(chase_data['Wavelength (nm)'])
    
    # A*0.1**B  
    
   # s = pd.Series(data_nc[keys[i]].values, index = data_nc.hplc_time)
   # s =
            
     #       ds_28 = ds_28[data_nc_28['hplc_depth'].values < 10]
      #      ds_28 = ds_28.groupby(ds_28.index).mean() # no reps!
     #       
            
     #       data_nc_28['hplc_Tot_Chl_b'][data_nc_28['hplc_depth'].values < 10]
   # data_nc_28("time.dayofyear").mean("hplc_time")
    #
    #ds_29 = pd.Series(data_nc_29[pig_keys[i]].values, index = data_nc_29.hplc_time)
   # ds_29 = ds_29[data_nc_29['hplc_depth'].values < 10]
   # 
   
   
prov= ['NADR', 'NASE', 'NASW', 'NATR', 'WTRA', 'SATL', 'SANT', 'SSTC'] # FKLD==1, ANTA -4, NECS = 8. do not include
                           
pig_keys = ['hplc_Tot_Chl_a','hplc_Tot_Chl_b', 'hplc_Tot_Chl_c' , 
            'hplc_Allo','hplc_Alpha-beta-Car', 'hplc_But-fuco', 
            'hplc_Diadino','hplc_Diato','hplc_Fuco',
            'hplc_Hex-fuco','hplc_Perid', 'hplc_Zea']
 
labels =  ['Tot_Chl_a', 'Tot_Chl_b', 'Tot_Chl_c' , 
              'Allo','alpha-beta-Car', 'But-fuco', 
              'Diadino','Diato','Fuco',
              'Hex-fuco','Perid', 'Zea']

letters = ['A', 'B', 'C' ,'D' ,'E' ,'F','G' , 'H']


plt.figure(figsize=(13,13))    
for i in range(len(prov)):
    plt.subplot(3, 3, i+1)  
    plt.title(prov[i] + ': N = ' + str(np.sum(LH_HPLC['LH_Province'] ==str(prov[i]))))

    pig_vec = []
    for j in range(len(pig_keys)):
        pig_j = np.array(df_hplc[pig_keys[j]])[LH_HPLC['LH_Province'] ==str(prov[i])]
        pig_j[pig_j==0] = np.nan
        pig_vec.append(pig_j)                                                                   

    bp = plt.boxplot(pig_vec ,showfliers=True,patch_artist=True, medianprops=dict(color='black'), whis=[10,90]) 
    plt.yscale('log')


    ax = plt.gca()
    ax.set_xticks([])  
    ax.set_ylim([0.0005, 5])
    
    plt.text(.90, .95,  letters[i], ha='left', va='top', transform=ax.transAxes,fontsize=22)  


    colors = cm.Paired(np.linspace(0,1,len(pig_keys))) 
    patches = []
    for k in range(len(pig_keys)):
        bp['boxes'][k].set_facecolor(colors[k])
        color_string = cl.rgb2hex(colors[k])
        patches.append(mpatches.Patch(color=color_string, label=labels[k]))
    
        plt.tight_layout(pad=1.2)
   
plt.subplot(3,3,9)
plt.legend(handles=patches,fontsize=14,loc=2)
ax = plt.gca()
ax.set_xticks([])  
ax.set_yticks([])  
ax.set_axis_off()