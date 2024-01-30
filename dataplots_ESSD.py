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
    plt.legend(markerscale=12 ,loc=4)            
    
    
    plt.text(.33, .92,   'NADR', ha='left', va='top', color='white', transform=ax.transAxes,fontsize=16)
    plt.text(.12, .80,   'NASW', ha='left', va='top', color='white', transform=ax.transAxes,fontsize=16)
    plt.text(.53, .78,   'NASE', ha='left', va='top', color='white', transform=ax.transAxes,fontsize=16)
    plt.text(.16, .65,   'NATR', ha='left', va='top', color='white', transform=ax.transAxes,fontsize=16)
    plt.text(.26, .53,   'WTRA', ha='left', va='top', color='white', transform=ax.transAxes,fontsize=16)
    plt.text(.53, .28,   'SATL', ha='left', va='top', color='white', transform=ax.transAxes,fontsize=16)
    plt.text(.50, .13,   'SSTC', ha='left', va='top', color='white', transform=ax.transAxes,fontsize=16)
    plt.text(.42, .06,   'SANT', ha='left', va='top', color='white', transform=ax.transAxes,fontsize=16)    

    
    return


# IOP plot functions

def _9cruise_IOPsum(field, field2):
    
    combined_array = np.concatenate([np.array(data_nc_19[field].values), data_nc_22[field].values, data_nc_23[field].values, 
                                     data_nc_24[field].values,  data_nc_24[field2].values, 
                                     data_nc_25[field].values,  data_nc_25[field2].values, 
                                     data_nc_26[field].values, data_nc_27[field].values, 
                                     data_nc_28[field].values, data_nc_29[field].values], axis=0)
                    
    return combined_array           


def _9cruise_LH_mask(field = 'LH_Province'):
    'double-count for AMT 25 and 24 due to 2 AC-S systems'
    
    field = 'LH_Province'
    combined_LH = np.concatenate([mask_19[field], mask_22[field], mask_23[field], 
                                         mask_24[field],  mask_24[field], 
                                         mask_25[field], mask_25[field], mask_26[field], 
                                         mask_27[field], mask_28[field], mask_29[field]], axis=0)
    return combined_LH


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
        
    colors = cm.Paired(np.linspace(0,1,12)) 
    c1=colors[1]
    c2=colors[3]
    c3=colors[5]
    
    
    plt.figure(figsize=(15,9))
    plt.subplot(2,3,1)
    plt.rcParams.update({'font.size': 16})
    plt.plot(data_nc_29['acs_wv'],np.nanmedian(apu_mat,axis=0),label='Median', color=c1)
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(apu_mat,25,axis=0), color=c1, linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(apu_mat,75,axis=0), label='Quartiles', color=c1, linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(apu_mat,90,axis=0), label='10$^{th}$ & 90$^{th}$ percentiles', color=c1,linestyle='dotted')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(apu_mat,10,axis=0), color=c1,linestyle='dotted')
    plt.xlim(400,720)
    plt.ylim(0,0.015)
    #plt.xlabel('Wavelength [nm]')
    plt.ylabel('$\sigma_{a_p}$ [m$^{-1}$]')
    plt.legend(fontsize=12)
    ax = plt.gca()
    plt.text(.05, .95,  'A', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
    
    plt.subplot(2,3,2)
    plt.rcParams.update({'font.size': 16})
    plt.plot(data_nc_29['acs_wv'],np.nanmedian(bpu_mat,axis=0),label='Median', color=c2)
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(bpu_mat,25,axis=0), color=c2, linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(bpu_mat,75,axis=0), label='Quartiles', color=c2, linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(bpu_mat,90,axis=0), label='10$^{th}$ & 90$^{th}$ percentiles', color=c2,linestyle='dotted')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(bpu_mat,10,axis=0), color=c2, linestyle='dotted')
    plt.xlim(400,720)
    plt.ylim(0,0.015)
    #plt.xlabel('Wavelength [nm]')
    #plt.ylabel('Uncertainty [m$^{-1}$]')
    plt.ylabel('$\sigma_{b_p}$ [m$^{-1}$]')
    plt.legend(fontsize=12)
    ax = plt.gca()
    plt.text(.05, .95,  'B', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
    
    
    plt.subplot(2,3,3)
    plt.rcParams.update({'font.size': 16})
    plt.plot(data_nc_29['acs_wv'],np.nanmedian(cpu_mat,axis=0),label='Median', color=c3)
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(cpu_mat,25,axis=0), color=c3, linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(cpu_mat,75,axis=0), label='Quartiles', color=c3, linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(cpu_mat,90,axis=0), label='10$^{th}$ & 90$^{th}$ percentiles', color=c3,linestyle='dotted')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(cpu_mat,10,axis=0), color=c3,linestyle='dotted')
    plt.xlim(400,720)
    plt.ylim(0,0.015)
    plt.ylabel('$\sigma_{c_p}$ [m$^{-1}$]')
    #plt.xlabel('Wavelength [nm]')
    #plt.ylabel('Uncertainty [m$^{-1}$]')
    plt.legend(fontsize=12)
    ax = plt.gca()
    plt.text(.05, .95,  'C', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
    
    
    
    plt.subplot(2,3,4)
    plt.rcParams.update({'font.size': 16})
    plt.plot(data_nc_29['acs_wv'],np.nanmedian(100*apu_mat/np.abs(ap_mat),axis=0),label='Median', color=c1)
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*apu_mat/np.abs(ap_mat),25,axis=0), color=c1, linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*apu_mat/np.abs(ap_mat),75,axis=0), label='Quartiles', color=c1, linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*apu_mat/np.abs(ap_mat),90,axis=0), label='10$^{th}$ & 90$^{th}$ percentiles', color=c1,linestyle='dotted')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*apu_mat/np.abs(ap_mat),10,axis=0), color=c1,linestyle='dotted')
    plt.xlim(400,750)
    plt.ylim(0,100)
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('$\sigma_{a_{p}}$ [%]')
    ax = plt.gca()
    plt.text(.05, .95,  'D', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
    # plt.legend(fontsize=12)
    
    
    plt.subplot(2,3,5)
    plt.rcParams.update({'font.size': 16})
    plt.plot(data_nc_29['acs_wv'],np.nanmedian(100*bpu_mat/np.abs(bp_mat),axis=0),label='Median', color=c2)
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*bpu_mat/np.abs(bp_mat),25,axis=0), color=c2, linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*bpu_mat/np.abs(bp_mat),75,axis=0), label='Quartiles', color=c2, linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*bpu_mat/np.abs(bp_mat),90,axis=0), label='10$^{th}$ & 90$^{th}$ percentiles', color=c2,linestyle='dotted')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*bpu_mat/np.abs(bp_mat),10,axis=0), color=c2,linestyle='dotted')
    plt.xlim(400,750)
    plt.ylim(0,30)
    plt.xlabel('Wavelength [nm]')
    ax = plt.gca()
    plt.ylabel('$\sigma_{b_{p}}$ [%]')
    plt.text(.05, .95,  'E', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
    #plt.ylabel('Percentage uncertainty [%]')
    #plt.legend(fontsize=12)
    
    
    plt.subplot(2,3,6)
    plt.rcParams.update({'font.size': 16})
    plt.plot(data_nc_29['acs_wv'],np.nanmedian(100*cpu_mat/np.abs(cp_mat),axis=0),label='Median', color=c3)
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*cpu_mat/np.abs(cp_mat),25,axis=0), color=c3, linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*cpu_mat/np.abs(cp_mat),75,axis=0), label='Quartiles', color=c3, linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*cpu_mat/np.abs(cp_mat),90,axis=0), label='10$^{th}$ & 90$^{th}$ percentiles', color=c3,linestyle='dotted')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(100*cpu_mat/np.abs(cp_mat),10,axis=0), color=c3,linestyle='dotted')
    plt.xlim(400,750)
    plt.ylim(0,30)
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('$\sigma_{c_{p}}$ [%]')
    ax = plt.gca()
    plt.text(.05, .95,  'F', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
    #plt.ylabel('Percentage uncertainty [%]')
    #plt.legend(fontsize=12)
    
    plt.tight_layout(pad=1.6)
    
    return


def plot_median_ap_province():
        'plots median ap, bp, cp by province'
                
        mask = _9cruise_LH_mask(field = 'LH_Province')
  
        ap_mat = _9cruise_IOPsum('acs_ap', 'acs2_ap')
        bp_mat = _9cruise_IOPsum('acs_bp', 'acs2_bp')
        cp_mat = _9cruise_IOPsum('acs_cp', 'acs2_cp')
          
        ap_int = int_norm(ap_mat)
        bp_int = int_norm(bp_mat)
        cp_int = int_norm(cp_mat)
        
        N_NADR = np.sum(~np.isnan((ap_mat[mask=='NADR',50]))) # 50 is arbitrary
        N_NASE = np.sum(~np.isnan((ap_mat[mask=='NASE',50])))
        N_NASW = np.sum(~np.isnan((ap_mat[mask=='NASW',50])))
        N_NATR = np.sum(~np.isnan((ap_mat[mask=='NATR',50])))
        N_WTRA = np.sum(~np.isnan((ap_mat[mask=='WTRA',50])))
        N_SATL = np.sum(~np.isnan((ap_mat[mask=='SATL',50])))
        N_SSTC = np.sum(~np.isnan((ap_mat[mask=='SSTC',50])))
        N_SANT = np.sum(~np.isnan((ap_mat[mask=='SANT',50])))


        data_nc = data_nc_26 # also abitraty -use to get wavelength
        
        one_over_lamba =(1/data_nc['acs_wv'].values)/np.mean(1/data_nc['acs_wv'].values)


        plt.figure(figsize=(15,8))
        plt.subplot(2,3,1)
        plt.rcParams.update({'font.size': 16})
        
        
        colors = cm.Paired(np.linspace(0,1,10))
    
        
        data_nc['acs_wv']
        
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_mat[mask=='NADR',:],axis=0),label='NADR',color=colors[0])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_mat[mask=='NASW',:],axis=0),label='NASW',color=colors[1])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_mat[mask=='NASE',:],axis=0),label='NASE',color=colors[2])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_mat[mask=='NATR',:],axis=0),label='NATR',color=colors[3])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_mat[mask=='WTRA',:],axis=0),label='WTRA',color=colors[4])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_mat[mask=='SATL',:],axis=0),label='SATL',color=colors[5])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_mat[mask=='SSTC',:],axis=0),label='SSTC',color=colors[6])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_mat[mask=='SANT',:],axis=0),label='SANT',color=colors[7])
        
        #plt.legend(loc=(1.04, 0))
        plt.xlim(400,720)
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$a_{p}$ [m$^{-1}$]')
        ax = plt.gca()
        plt.text(.86, .95,  'A', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
                
        
        plt.subplot(2,3,2)
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_mat[mask=='NADR',:],axis=0),label='NADR',color=colors[0])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_mat[mask=='NASW',:],axis=0),label='NASW',color=colors[1])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_mat[mask=='NASE',:],axis=0),label='NASE',color=colors[2])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_mat[mask=='NATR',:],axis=0),label='NATR',color=colors[3])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_mat[mask=='WTRA',:],axis=0),label='WTRA',color=colors[4])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_mat[mask=='SATL',:],axis=0),label='SATL',color=colors[5])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_mat[mask=='SSTC',:],axis=0),label='SSTC',color=colors[6])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_mat[mask=='SANT',:],axis=0),label='SANT',color=colors[7])
        
        #plt.legend(loc=(1.04, 0))
        plt.xlim(400,720)
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$b_{p}$ [m$^{-1}$]')
        ax = plt.gca()
        plt.text(.86, .95,   'B', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
        
         
        plt.subplot(2,3,3)   
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_mat[mask=='NADR',:],axis=0),label='NADR: N = ' + str(N_NADR),color=colors[0])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_mat[mask=='NASW',:],axis=0),label='NASW: N = ' + str(N_NASW),color=colors[1])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_mat[mask=='NASE',:],axis=0),label='NASE: N = ' + str(N_NASE),color=colors[2])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_mat[mask=='NATR',:],axis=0),label='NATR: N = ' + str(N_NATR),color=colors[3])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_mat[mask=='WTRA',:],axis=0),label='WTRA: N = ' + str(N_WTRA),color=colors[4])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_mat[mask=='SATL',:],axis=0),label='SATL: N = ' + str(N_SATL),color=colors[5])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_mat[mask=='SSTC',:],axis=0),label='SSTC: N = ' + str(N_SSTC),color=colors[6])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_mat[mask=='SANT',:],axis=0),label='SANT: N = ' + str(N_SANT),color=colors[7])
        plt.plot(0, 0, color='black', linestyle = 'dashed', label ='$\lambda^{-1}$')
        plt.legend(loc=(1.04, 0))
        plt.xlim(400,720)
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$c_{p}$ [m$^{-1}$]')
        ax = plt.gca()
        plt.text(.86, .95,   'C', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
        
        
        plt.subplot(2,3,4)
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_int[mask=='NADR',:],axis=0),label='NADR',color=colors[0])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_int[mask=='NASW',:],axis=0),label='NASW',color=colors[1])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_int[mask=='NASE',:],axis=0),label='NASE',color=colors[2])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_int[mask=='NATR',:],axis=0),label='NATR',color=colors[3])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_int[mask=='WTRA',:],axis=0),label='WTRA',color=colors[4])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_int[mask=='SATL',:],axis=0),label='SATL',color=colors[5])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_int[mask=='SSTC',:],axis=0),label='SSTC',color=colors[6])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_int[mask=='SANT',:],axis=0),label='SANT',color=colors[7])
        #plt.legend(loc=(1.04, 0))
        plt.xlim(400,720)
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$<a_{p}>$')
        ax = plt.gca()
        plt.text(.86, .95, 'D', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
        
        
        plt.subplot(2,3,5)
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_int[mask=='NADR',:],axis=0),label='NADR',color=colors[0])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_int[mask=='NASW',:],axis=0),label='NASW',color=colors[1])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_int[mask=='NASE',:],axis=0),label='NASE',color=colors[2])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_int[mask=='NATR',:],axis=0),label='NATR',color=colors[3])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_int[mask=='WTRA',:],axis=0),label='WTRA',color=colors[4])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_int[mask=='SATL',:],axis=0),label='SATL',color=colors[5])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_int[mask=='SSTC',:],axis=0),label='SSTC',color=colors[6])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_int[mask=='SANT',:],axis=0),label='SANT',color=colors[7])
        plt.plot(data_nc['acs_wv'], one_over_lamba, color='black', linestyle = 'dashed')
        #plt.legend(loc=(1.04, 0))
        plt.xlim(400,720)
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$<b_{p}>$')
        ax = plt.gca()
        plt.text(.86, .95,  'E', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
        
         
        plt.subplot(2,3,6)   
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_int[mask=='NADR',:],axis=0),label='NADR',color=colors[0])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_int[mask=='NASW',:],axis=0),label='NASW',color=colors[1])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_int[mask=='NASE',:],axis=0),label='NASE',color=colors[2])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_int[mask=='NATR',:],axis=0),label='NATR',color=colors[3])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_int[mask=='WTRA',:],axis=0),label='WTRA',color=colors[4])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_int[mask=='SATL',:],axis=0),label='SATL',color=colors[5])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_int[mask=='SSTC',:],axis=0),label='SSTC',color=colors[6])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_int[mask=='SANT',:],axis=0),label='SANT',color=colors[7])
        plt.plot(data_nc['acs_wv'], one_over_lamba, color='black', linestyle = 'dashed')
        #plt.legend(loc=(1.04, 0))
        plt.xlim(400,720)
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('$<c_{p}>$')
        ax = plt.gca()
        plt.text(.86, .95,   'F', ha='left', va='top', transform=ax.transAxes,fontsize=24)


        plt.tight_layout()

        return

def plot_676_443():
 
    mask = _9cruise_LH_mask(field = 'LH_Province')
      
    ap_mat = _9cruise_IOPsum('acs_ap', 'acs2_ap')
    ap_int = int_norm(ap_mat)
    
    colors = cm.Paired(np.linspace(0,1,10))
    
    ap_443 = (ap_int[:, 22] + ap_int[:, 23])/2
    ap_676  = ap_int[:, 138]
    ap_rat = ap_676/ap_443
    
    prov= ['NADR', 'NASW',  'NASE' ,'NATR', 'WTRA', 'SATL', 'SSTC', 'SANT'] # FKLD==1, ANTA -4, NECS = 8. do not include
    
    rat_vec =[]
    for i in range(len(prov)):
        rat_i = np.array(ap_rat)[mask ==str(prov[i])]
        rat_i =  rat_i[~np.isnan(rat_i)]
        rat_vec.append(rat_i)                                                                   
    
    plt.figure()
    bp = plt.boxplot(rat_vec ,showfliers=True,patch_artist=True, medianprops=dict(color='black'), whis=[10,90],widths = 0.8) 
    plt.ylim(0,0.7)
    
    
    plt.ylabel('$a_{p}(676)/a_{p}(443)$')
    colors = cm.Paired(np.linspace(0,1,10))
    patches = []
    for k in range(len(prov)):
        bp['boxes'][k].set_facecolor(colors[k])
       
        color_string = cl.rgb2hex(colors[k])
        patches.append(mpatches.Patch(color=color_string, label=str(prov[k])))
        
        
   # plt.legend(handles=patches,fontsize=14,loc=2)
    ax = plt.gca()
    ax.set_xticklabels(labels= ['NADR', 'NASW',  'NASE' ,'NATR', 'WTRA', 'SATL', 'SSTC', 'SANT'],rotation=45)  
   # ax.set_xticks([])  
   # ax.set_yticks([])  
   # ax.set_axis_off()
    
    
    return



# legacy

def plot_median_ap_cp(data_nc):
        
    plt.figure()
    plt.rcParams.update({'font.size': 16})
    plt.plot(data_nc['acs_wv'],np.nanmedian(data_nc['acs_ap'],axis=0),label='$a_{p}$')
    plt.plot(data_nc['acs_wv'],np.nanmedian(data_nc['acs_cp'],axis=0),label='$c_{p}$')
    plt.xlim(400,750)
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('[m$^{-1}$]')
    plt.legend()
    
    return



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
                'hplc_Perid', 'hplc_Zea', 'hplc_PPC', 'hplc_PSC']
       
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
    
   # df_hplc_combined = pd.concat([df_hplc_29, df_hplc_28, df_hplc_27, df_hplc_22, df_hplc_19 ]) #- removing PML
    
    return  df_hplc_combined



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
       
    prov= ['NADR', 'NASW','NASE',  'NATR', 'WTRA', 'SATL', 'SSTC', 'SANT'] # FKLD==1, ANTA -4, NECS = 8. do not include
                               
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
    plt.rcParams.update({'font.size': 18})
    for i in range(len(prov)):
        plt.subplot(3, 3, i+1)  
        plt.title(prov[i] + ': N = ' + str(np.sum(LH_HPLC['LH_Province'] ==str(prov[i]))))
    
        pig_vec = []
        for j in range(len(pig_keys)):
            pig_j = np.array(df_hplc[pig_keys[j]])[LH_HPLC['LH_Province'] ==str(prov[i])]
            pig_j[pig_j==0] = np.nan
            pig_vec.append(pig_j)                                                                   
    
        bp = plt.boxplot(pig_vec ,showfliers=True,patch_artist=True, medianprops=dict(color='black'), whis=[10,90],widths = 0.8) 
        plt.yscale('log')
    
        
        ax = plt.gca()
        ax.set_xticks([])  
        ax.set_ylim([0.0005, 5])
        if i==0:
            plt.ylabel('Concentration [mg m$^{-3}]$')
        if i==3:
            plt.ylabel('Concentration [mg m$^{-3}]$')
        if i==6:
            plt.ylabel('Concentration [mg m$^{-3}]$')
            
    
        
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



def pigs_byprov_ratio(): 
       
    prov= ['NADR', 'NASW','NASE',  'NATR', 'WTRA', 'SATL', 'SSTC', 'SANT'] # FKLD==1, ANTA -4, NECS = 8. do not include
                               
    pig_keys = ['hplc_Tot_Chl_b', 'hplc_Tot_Chl_c' , 
                'hplc_Allo','hplc_Alpha-beta-Car', 'hplc_But-fuco', 
                'hplc_Diadino','hplc_Diato','hplc_Fuco',
                'hplc_Hex-fuco','hplc_Perid', 'hplc_Zea']
     
    labels =  ['Tot_Chl_b', 'Tot_Chl_c' , 
                  'Allo','alpha-beta-Car', 'But-fuco', 
                  'Diadino','Diato','Fuco',
                  'Hex-fuco','Perid', 'Zea']
    
    letters = ['A', 'B', 'C' ,'D' ,'E' ,'F','G' , 'H']
    
    
    plt.figure(figsize=(13,13))    
    plt.rcParams.update({'font.size': 18})
    for i in range(len(prov)):
        plt.subplot(3, 3, i+1)  
        plt.title(prov[i] + ': N = ' + str(np.sum(LH_HPLC['LH_Province'] ==str(prov[i]))))
    
        pig_vec = []
        for j in range(len(pig_keys)):
            pig_j = np.array(df_hplc[pig_keys[j]]/df_hplc['hplc_Tot_Chl_a'])[LH_HPLC['LH_Province'] ==str(prov[i])]
            pig_j[pig_j==0] = np.nan
            pig_vec.append(pig_j)                                                                   
    
        bp = plt.boxplot(pig_vec ,showfliers=True,patch_artist=True, medianprops=dict(color='black'), whis=[10,90],widths = 0.8) 
        plt.yscale('log')
    
        
        ax = plt.gca()
        ax.set_xticks([])  
        ax.set_ylim([0.003, 1])
        if i==0:
            plt.ylabel('Concentration ratio')
        if i==3:
            plt.ylabel('Concentration ratio')
        if i==6:
            plt.ylabel('Concentration ratio')
            
        plt.text(.05, .08,  letters[i], ha='left', va='top', transform=ax.transAxes,fontsize=22)  
    
        colors = cm.Paired(np.linspace(0,1,len(pig_keys)+1)) 
        patches = []
        for k in range(len(pig_keys)):
            bp['boxes'][k].set_facecolor(colors[k+1])
            color_string = cl.rgb2hex(colors[k+1])
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



def hplc_ap_match_up_V3(data_nc, two_acs_systems=False):
    
    ' extracts ap spectra that match hplc timestamps - median average over 30 min bin.'
    ' this is designed to replicate what was done in the tot-chl-a match up procedure'
    ' important to first remove: (i) replicates, (ii) deeper waters'
    ' user specifies wavelength and pigment'
    
    
    # hplc data frame
    df_hplc = _filter_combine_pigs(data_nc)
         
    if two_acs_systems == False:  
        
        # ap data frame
        ap = data_nc['acs_ap'].data       
        wv = data_nc['acs_wv'].values
        
        df_ap = pd.DataFrame(index = data_nc['time']) # pandas data frame format for down-sampled spectra
        for i in range(len(wv)):
            ap_f = sg.medfilt(ap[:,i],kernel_size=31) # 
            df_ap[str(wv[i])] = ap_f
        
            
        df_ap = df_ap.interpolate('index',limit=1).reindex(df_hplc.index,method='nearest',tolerance='30min')
    
         
        # repat for integral normalized
        ap_int = int_norm(ap)
        
        df_ap_int = pd.DataFrame(index = data_nc['time']) # pandas data frame format for down-sampled spectra
        for i in range(len(wv)):
            ap_int_f = sg.medfilt(ap_int[:,i],kernel_size=31) # 
            df_ap_int[str(wv[i])] = ap_int_f
        
            
        df_ap_int = df_ap_int.interpolate('index',limit=1).reindex(df_hplc.index,method='nearest',tolerance='30min')
        
    if two_acs_systems == True:  
        
        # ap data frame
        ap = data_nc['acs_ap'].data 
      #  ap = int_norm(ap)
        ap2 = data_nc['acs2_ap'].data 
      #  ap2 = int_norm(ap2)
        # integral normalization to go here
        wv = data_nc['acs_wv'].values
        
        mat_ap = np.zeros([2,len(ap), len(ap.T)])
        mat_ap[0,:,:] = ap
        mat_ap[1,:,:] = ap2 
        combined_ap = np.nanmean(mat_ap, axis=0)
                       
        
        df_ap = pd.DataFrame(index = data_nc['time']) # pandas data frame format for down-sampled spectra
        for i in range(len(wv)):
            ap_f = sg.medfilt(combined_ap[:,i],kernel_size=31) # 
            df_ap[str(wv[i])] = ap_f
        
        df_ap = df_ap.interpolate('index',limit=1).reindex(df_hplc.index,method='nearest',tolerance='30min')
            
        
        ap_int = int_norm(combined_ap)
        
        df_ap_int = pd.DataFrame(index = data_nc['time']) # pandas data frame format for down-sampled spectra
        for i in range(len(wv)):
            ap_int_f = sg.medfilt(ap_int[:,i],kernel_size=31) # 
            df_ap_int[str(wv[i])] = ap_int_f
        
            
        df_ap_int = df_ap_int.interpolate('index',limit=1).reindex(df_hplc.index,method='nearest',tolerance='30min')
        
        
        

    return  df_hplc, df_ap, df_ap_int


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


###############################################################################
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
        tot_acs= np.concatenate([0.014*data_nc_19.acx_chl_debiased.acx_chl, 
        0.014*data_nc_22.acs_chl_debiased.acs_chl , 
        0.014* data_nc_23.acs_chl_debiased.acs_chl , 
        0.014*data_nc_24.acs_chl_debiased.acs_chl , 
        0.014*data_nc_24.acs2_chl_debiased.acs_chl,
        0.014*data_nc_25.acs_chl_debiased.acs_chl ,
        #  data_nc_25.acs2_chl_debiased.acs_chl, 
        0.014*data_nc_26.acs_chl_debiased.acs_chl , 
        0.014*data_nc_27.acs_chl_debiased.acs_chl , 
        0.014* data_nc_28.acx_chl_debiased.acx_chl , 
        0.014*data_nc_29.acs_chl_debiased.acs_chl ]) 
            
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
        
        #linear_mod = scipy.stats.linregress(np.log10(0.014*tot_acs),np.log10(tot_chl))
                                                  
        linear_mod = scipy.stats.linregress(np.log10(tot_acs),np.log10(tot_chl)) #- no A
                                                  
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
           
        X = np.arange(0.001,10,0.010)
        Y = X*0.014
        Y2 = 10**(A+B*np.log10(Y))
        colors = cm.Paired(np.linspace(0,1,10))
           
           
        plt.figure(figsize=(6, 6))
        plt.rcParams.update({'font.size': 14})
        plt.scatter(0.014*data_nc_19.acx_chl_debiased.acx_chl, data_nc_19.acx_chl_debiased.HPLC_Tot_chla, color=colors[0],alpha=0.7, s=6, label = 'AMT 19')
        plt.scatter(0.014*data_nc_22.acs_chl_debiased.acs_chl, data_nc_22.acs_chl_debiased.HPLC_Tot_chla, color=colors[1], alpha=0.7, s=6, label = 'AMT 22')
        plt.scatter(0.014*data_nc_23.acs_chl_debiased.acs_chl, data_nc_23.acs_chl_debiased.HPLC_Tot_chla, color=colors[2], alpha=0.7, s=6, label = 'AMT 23')
        plt.scatter(0.014*data_nc_24.acs_chl_debiased.acs_chl, data_nc_24.acs_chl_debiased.HPLC_Tot_chla, color=colors[3], alpha=0.7, s=6)
        plt.scatter(0.014*data_nc_24.acs2_chl_debiased.acs_chl, data_nc_24.acs2_chl_debiased.HPLC_Tot_chla,  color=colors[3], alpha=0.7, s=6, label = 'AMT 24')
        plt.scatter(0.014*data_nc_25.acs_chl_debiased.acs_chl, data_nc_25.acs_chl_debiased.HPLC_Tot_chla, color=colors[4], alpha=0.7, s=6, label = 'AMT 25')
        plt.scatter(0.014*data_nc_26.acs_chl_debiased.acs_chl, data_nc_26.acs_chl_debiased.HPLC_Tot_chla, color=colors[5], alpha=0.7, s=6, label = 'AMT 26')
        plt.scatter(0.014*data_nc_27.acs_chl_debiased.acs_chl, data_nc_27.acs_chl_debiased.HPLC_Tot_chla, color=colors[6], alpha=0.7, s=6, label = 'AMT 27')
        plt.scatter(0.014*data_nc_28.acx_chl_debiased.acx_chl, data_nc_28.acx_chl_debiased.HPLC_Tot_chla,  color=colors[7], alpha=0.7, s=6, label = 'AMT 28')
        plt.scatter(0.014*data_nc_29.acs_chl_debiased.acs_chl, data_nc_29.acs_chl_debiased.HPLC_Tot_chla, color=colors[9], alpha=0.7, s=6, label = 'AMT 29')
        plt.yscale('log')
        plt.xscale('log')    
        plt.grid('on', ls='--')
        plt.plot(Y,X,'k',linestyle='dotted', label='0.014:1')
        plt.plot(Y,Y2,'magenta',linestyle='dashed', label='Best fit')
        plt.xlabel('$a_{ph}(676)$ [m$^{-1}]$')
        #plt.ylabel('$C_{a}$(HPLC)  [mg m$^{-3}]$')
        plt.ylabel('Tot_Chl_a (HPLC)  [mg m$^{-3}]$')
        plt.legend(fontsize=11)
  
        plt.axis('equal')
    
        plt.xlim(0.0003,0.01)
        plt.ylim(0.01,3)

        return
        

def prcrng(x):
    return (np.nanpercentile(x,84) - np.nanpercentile(x,16))/2



def Tchl_resid(data_nc): 
    'calculates residual ap, after removing Tchl baseline (follows  Chase 2013/Cael 2020)'
    
    data = data_nc.acs_ap.values
    Tchl = data_nc.acs_chl_debiased
    
    chase = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/Chase_AB.csv'
    chase_data = pd.read_csv(chase) 
    
    A = np.array(chase_data['Achl (m-1)'])
    B = np.array(chase_data['Bchl (unitless)'])
    wl_chase = np.array(chase_data['Wavelength (nm)'])
    
    data_resid = np.nan*np.ones([len(data),len(wl_chase)])
    for i in range(len(data)):
        Tchl_vec =Tchl[i].values*np.ones(len(wl_chase))
        data_resid[i,:] =  data[i,3:151] - (A*Tchl_vec**B)
  
    return data_resid, wv_chase



def Tchl_resid_with_intnorm(data_nc): 
    'calculates residual ap, after removing Tchl baseline (follows  Chase 2013/Cael 2020)'

    data = data_nc.acs_ap.values
    Tchl = data_nc.acs_chl_debiased
    
    chase = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/Chase_AB.csv'
    chase_data = pd.read_csv(chase) 
    
    A = np.array(chase_data['Achl (m-1)'])
    B = np.array(chase_data['Bchl (unitless)'])
    wl_chase = np.array(chase_data['Wavelength (nm)'])
    
    
    Tchl_vec = np.nan*np.ones([len(data),len(wl_chase)])
    for i in range(len(data)):
        Tchl_vec[i,:] = Tchl[i].values*np.ones(len(wl_chase))
    
    Tchl_vec_norm  = int_norm(Tchl_vec) 
    data_norm = int_norm(data) 
    
    for i in range(len(data)):
        data_resid[i,:] =  data_norm[i,3:151] - Tchl_vec_norm[i,:]
      
    return data_resid, wv_chase


    
def _relative_ACS_time(data_nc):

    timestamp = []
    for i in range(0,len(data_nc['time'])):
        str_i = '2023'+ str(data_nc['time'].values[i])[4:19] # create fake date all in same year
        timestamp.append(datetime.datetime.strptime(str_i, '%Y-%m-%dT%H:%M:%S'))
    timestamp = np.array(timestamp)

    return timestamp


    
def _relative_HPLC_time(data_nc):

    timestamp = []
    for i in range(0,len(data_nc['hplc_time'])):
        str_i = '2023'+ str(data_nc['hplc_time'].values[i])[4:19] # create fake date all in same yearep 
        timestamp.append(datetime.datetime.strptime(str_i, '%Y-%m-%dT%H:%M:%S'))
    timestamp = np.array(timestamp)
    timestamp = timestamp[data_nc['hplc_depth'].values < 10]  # no de

    return timestamp

    

def _color_by_prov(ts, tsh, mask, data_nc, plot_index, two_systems=False):
    
    prov_extra= ['NECS', 'FKLD', 'ANTA']
   
    if plot_index == 1 or plot_index == 8:
        for i in range(len(prov_extra)):
            plt.plot_date(plot_index*np.array(1*(~np.isnan(data_nc['acx_chl_debiased'])))[mask==prov_extra[i]],ts[mask==prov_extra[i]], marker="s",xdate=False, ms=10,color='black')
    else:
       for i in range(len(prov_extra)):
           plt.plot_date(plot_index*np.array(1*(~np.isnan(data_nc['acs_chl_debiased'])))[mask==prov_extra[i]],ts[mask==prov_extra[i]], marker="s", xdate=False, ms=10,color='gray')
           if two_systems ==True:
               combined = np.array(1*(~np.isnan(data_nc['acs2_chl_debiased']))+1*(~np.isnan(data_nc['acs_chl_debiased'])))
               combined[combined>0]=1
               plt.plot_date(plot_index*combined[mask==prov_extra[i]],ts[mask==prov_extra[i]], xdate=False, marker="s", ms=10,color='gray')

    prov= ['NADR', 'NASE', 'NASW', 'NATR', 'WTRA', 'SATL', 'SSTC', 'SANT'] 
   
    colors = cm.Paired(np.linspace(0,1,10))
    
    if plot_index == 1 or plot_index == 8:
        for i in range(len(prov)):
            plt.plot_date(plot_index*np.array(1*(~np.isnan(data_nc['acx_chl_debiased'])))[mask==prov[i]],ts[mask==prov[i]],xdate=False, marker="s",ms=10,color=colors[i])
        if plot_index==9:
                plt.plot_date(plot_index*np.array(1*(~np.isnan(data_nc['acx_chl_debiased'])))[mask==prov[i]],ts[mask==prov[i]], xdate=False, marker="s", ms=10,color=colors[i], label = 'ACS: ' + str(prov[i]))
    else:
       for i in range(len(prov)):
           plt.plot_date(plot_index*np.array(1*(~np.isnan(data_nc['acs_chl_debiased'])))[mask==prov[i]], ts[mask==prov[i]], xdate=False, marker="s",ms=10, color=colors[i])
           if plot_index ==9:
               plt.plot_date(plot_index*np.array(1*(~np.isnan(data_nc['acs_chl_debiased'])))[mask==prov[i]], ts[mask==prov[i]], xdate=False, marker="s", ms=10, color=colors[i], label = 'ACS: ' + str(prov[i]))
           if two_systems ==True:
               combined = np.array(1*(~np.isnan(data_nc['acs2_chl_debiased']))+1*(~np.isnan(data_nc['acs_chl_debiased'])))
               combined[combined>0]=1
               plt.plot_date(plot_index*combined[mask==prov[i]],ts[mask==prov[i]], xdate=False, ms=10, color=colors[i])
    

    plt.plot_date(plot_index*np.ones(len(tsh)), tsh, xdate=False, ms=3, color='black')

    if plot_index ==9:
        plt.plot_date(plot_index*np.array(1*(~np.isnan(data_nc['acs_chl_debiased'])))[mask==prov_extra[0]],ts[mask==prov_extra[0]], marker="s", xdate=False, ms=12,color='gray',label='ACS: OTHER PROVINCES')
        plt.plot_date(plot_index*np.ones(len(tsh)),tsh, xdate=False, ms=3, color='black', label='HPLC')

        plt.legend(fontsize=14)
        
    return 

def _AMT_timeline():
    
    plt.figure(figsize =(9,14))
    ts_19 = _relative_ACS_time(data_nc_19)
    tsh_19 = _relative_HPLC_time(data_nc_19)
    _color_by_prov(ts_19,tsh_19, mask_19['LH_Province'], data_nc_19, 1)

    ts_22 = _relative_ACS_time(data_nc_22)
    tsh_22 = _relative_HPLC_time(data_nc_22)
    _color_by_prov(ts_22,tsh_22, mask_22['LH_Province'], data_nc_22, 2)

    ts_23 = _relative_ACS_time(data_nc_23)
    tsh_23 = _relative_HPLC_time(data_nc_23)
    _color_by_prov(ts_23,tsh_23, mask_23['LH_Province'], data_nc_23, 3)

    ts_24 = _relative_ACS_time(data_nc_24)
    tsh_24 = _relative_HPLC_time(data_nc_24)
    _color_by_prov(ts_24, tsh_24,mask_24['LH_Province'], data_nc_24, 4, two_systems =True)

    ts_25 = _relative_ACS_time(data_nc_25)
    tsh_25 = _relative_HPLC_time(data_nc_25)
    _color_by_prov(ts_25,tsh_25, mask_25['LH_Province'], data_nc_25, 5, two_systems =True)

    ts_26 = _relative_ACS_time(data_nc_26)
    tsh_26 = _relative_HPLC_time(data_nc_26)
    _color_by_prov(ts_26, tsh_26, mask_26['LH_Province'], data_nc_26, 6)

    ts_27 = _relative_ACS_time(data_nc_27)
    tsh_27 = _relative_HPLC_time(data_nc_27)
    _color_by_prov(ts_27, tsh_27, mask_27['LH_Province'], data_nc_27, 7)

    ts_28 = _relative_ACS_time(data_nc_28)
    tsh_28 = _relative_HPLC_time(data_nc_28)
    _color_by_prov(ts_28,tsh_28, mask_28['LH_Province'], data_nc_28, 8)

    ts_29 = _relative_ACS_time(data_nc_29)
    tsh_29 = _relative_HPLC_time(data_nc_29)
    _color_by_prov(ts_29,tsh_29, mask_29['LH_Province'], data_nc_29, 9)

    plt.xlim(0.5,9.5)
    plt.gca().invert_yaxis()
    ax = plt.gca()
    ax.set_xticklabels(['','AMT 19','AMT 22','AMT 23','AMT 24','AMT 25','AMT 26','AMT 27','AMT 28','AMT 29'],rotation=45)  
    ax.yaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b-%d') )
    plt.ylabel('Date within each year')

    return


def _color_by_prov_chl(data_nc, mask, plot_index):
    
    plt.subplot(5,2,plot_index)
    letters = ['AMT 19', 'AMT 22' , 'AMT 23', 'AMT 24', 'AMT 25', 'AMT 26', 'AMT 27' ,'AMT 28', 'AMT 29']
        
    if plot_index == 1 or plot_index == 8:    
        acs_key ='acx_chl_debiased'
    else:
        acs_key ='acs_chl_debiased'
    if plot_index == 4 or plot_index == 5:    
        acs2_key = 'acs2_chl_debiased'
    
    
    
    prov_extra= ['NECS', 'FKLD', 'ANTA']
    for i in range(len(prov_extra)):
        plt.scatter(np.array(data_nc['time'])[mask==prov_extra[i]], np.array(data_nc[acs_key])[mask==prov_extra[i]], s=2, color='gray')
        if plot_index == 4 or plot_index == 5:    
            plt.scatter(np.array(data_nc['time'])[mask==prov_extra[i]], np.array(data_nc[acs2_key])[mask==prov_extra[i]], s=2, color='gray')
  
    prov= ['NADR', 'NASE', 'NASW', 'NATR', 'WTRA', 'SATL', 'SSTC', 'SANT'] 
    colors = cm.Paired(np.linspace(0,1,10))
    
    for i in range(len(prov)):
        plt.scatter(np.array(data_nc['time'])[mask==prov[i]], np.array(data_nc[acs_key])[mask==prov[i]],color=colors[i], s=2)
        if plot_index == 4 or plot_index == 5:  
            plt.scatter(np.array(data_nc['time'])[mask==prov[i]], np.array(data_nc[acs2_key])[mask==prov[i]],color=colors[i], s=2)
    
        
    if plot_index == 1 or plot_index == 8:
         plt.scatter(data_nc.acx_chl_debiased.match_up_dates, data_nc.acx_chl_debiased.HPLC_Tot_chla, s=10, color='black')
    elif plot_index == 4 or plot_index == 5: 
         plt.scatter(data_nc.acs2_chl_debiased.match_up_dates,data_nc.acs2_chl_debiased.HPLC_Tot_chla, s=10,color='black')    
    else:
         plt.scatter(data_nc.acs_chl_debiased.match_up_dates,data_nc.acs_chl_debiased.HPLC_Tot_chla, s=10,color='black')

    plt.yscale('log')
    plt.ylim(0.01,7)
    ax=plt.gca()
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%m-%d') )
    ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator(bymonthday=[1,5,10,15,20,25]))
    plt.xticks(rotation=45)

    plt.text(.02, .94,  letters[plot_index-1], ha='left', va='top', transform=ax.transAxes,fontsize=20) 

   # if plot_index == 8:
  

    return


def chl_time_series():
    
    fig = plt.figure(figsize=(15,20))
    plt.rcParams.update({'font.size': 18})
    _color_by_prov_chl(data_nc_19, mask_19['LH_Province'], 1)
    _color_by_prov_chl(data_nc_22, mask_22['LH_Province'], 2)
    _color_by_prov_chl(data_nc_23, mask_23['LH_Province'], 3)
    _color_by_prov_chl(data_nc_24, mask_24['LH_Province'], 4)
    _color_by_prov_chl(data_nc_25, mask_25['LH_Province'], 5)
    _color_by_prov_chl(data_nc_26, mask_26['LH_Province'], 6)
    _color_by_prov_chl(data_nc_27, mask_27['LH_Province'], 7)
    _color_by_prov_chl(data_nc_28, mask_28['LH_Province'], 8)
    _color_by_prov_chl(data_nc_29, mask_29['LH_Province'], 9)
    
    fig.supylabel('Tot_Chl_a concentration  [mg m$^{-3}]$')
    fig.supxlabel('Date (MM-DD) within each year')
    
    colors = cm.Paired(np.linspace(0,1,10)) 
    prov = ['ACS: NADR', 'ACS: NASE', 'ACS: NASW', 'ACS: NATR', 'ACS: WTRA', 'ACS: SATL', 'ACS: SSTC', 'ACS: SANT']     
    patches = []
    for k in range(8):
        color_string = cl.rgb2hex(colors[k])
        patches.append(mpatches.Patch(color=color_string, label=prov[k]))
    
    patches.append(mpatches.Patch(color='gray', label='ACS: OTHER'))
    #patches.append(mpatches.Patch(color='black', label='HPLC'))
    
    ax = plt.subplot(5,2,10)
    plt.legend(handles=patches,fontsize=14,loc=2)
    ax = plt.gca()
    ax.set_xticks([])  
    ax.set_yticks([])  
    ax.set_axis_off()
    ax = plt.gca()
    ax2 = ax.twinx()
    ax2.scatter(100,100, s=8,color='black', label='HPLC: ALL PROVINCES')    
    ax2.set_xlim(0,1)
    ax2.set_ylim(0,1)
    ax2.legend(loc=1,fontsize=14)
    plt.tight_layout()
    ax2.axis("off")

    return


def int_globalmed():
    
    ap_mat = _9cruise_IOPsum('acs_ap', 'acs2_ap')          
    ap_int = int_norm(ap_mat)
    global_med = np.nanmedian(ap_int,axis=0)

    return  global_med



def _ap_limitingcases():

    df_hplc_19, df_ap_19, df_ap_int_19 = hplc_ap_match_up_V3(data_nc_19, two_acs_systems=False)
    df_hplc_22, df_ap_22, df_ap_int_22 = hplc_ap_match_up_V3(data_nc_22, two_acs_systems=False)
    df_hplc_23, df_ap_23, df_ap_int_23 = hplc_ap_match_up_V3(data_nc_23, two_acs_systems=False)
    df_hplc_24, df_ap_24, df_ap_int_24 = hplc_ap_match_up_V3(data_nc_24, two_acs_systems=True)
    df_hplc_25, df_ap_25, df_ap_int_25 = hplc_ap_match_up_V3(data_nc_25, two_acs_systems=True)
    df_hplc_26, df_ap_26, df_ap_int_26 = hplc_ap_match_up_V3(data_nc_26, two_acs_systems=False)
    df_hplc_27, df_ap_27, df_ap_int_27 = hplc_ap_match_up_V3(data_nc_27, two_acs_systems=False)
    df_hplc_28, df_ap_28, df_ap_int_28 = hplc_ap_match_up_V3(data_nc_28, two_acs_systems=False)
    df_hplc_29, df_ap_29, df_ap_int_29 = hplc_ap_match_up_V3(data_nc_29, two_acs_systems=False)
    
    df_hplc_combined = pd.concat([df_hplc_29, df_hplc_28, df_hplc_27, df_hplc_26, df_hplc_25, df_hplc_24, df_hplc_23, df_hplc_22, df_hplc_19 ])
    df_ap_combined = pd.concat([df_ap_29, df_ap_28, df_ap_27, df_ap_26, df_ap_25, df_ap_24, df_ap_23, df_ap_22, df_ap_19 ])
    df_ap_int_combined = pd.concat([df_ap_int_29, df_ap_int_28, df_ap_int_27, df_ap_int_26, df_ap_int_25, df_ap_int_24, df_ap_int_23, df_ap_int_22, df_ap_int_19 ])
    
    
    df_hplc_combined = df_hplc_combined[np.isnan(df_ap_combined['400.0'])==0]
    df_ap_combined = df_ap_combined[np.isnan(df_ap_combined['400.0'])==0]
    df_ap_int_combined = df_ap_int_combined[np.isnan(df_ap_int_combined['400.0'])==0]
    
    global_med = int_globalmed()
    wv = data_nc_26['acs_wv'].values
    
   # pig_keys = ['hplc_Tot_Chl_b', 'hplc_Tot_Chl_c' , 
    #            'hplc_Allo','hplc_Alpha-beta-Car', 'hplc_But-fuco', 
     #           'hplc_Diadino','hplc_Diato','hplc_Fuco',
      #         'hplc_Hex-fuco','hplc_Perid', 'hplc_Zea']
    
   # labels =  [ 'Tot_Chl_b', 'Tot_Chl_c' , 
    #             'Allo','alpha-beta-Car', 'But-fuco', 
     #            'Diadino','Diato','Fuco',
      #           'Hex-fuco','Perid', 'Zea']
    
    pig_keys = ['hplc_Tot_Chl_b', 'hplc_Tot_Chl_c' , 
               'hplc_PPC', 'hplc_PSC']
              
    labels =  [ 'Tot_Chl_b', 'Tot_Chl_c' , 
                  'PPC', 'PSC']
    
        
    colors = cm.Paired(np.linspace(0,1,12))
    plt.rcParams.update({'font.size': 20})
    plt.figure(figsize=(18,18))
    plt
    plt.subplot(2,2,1)
    style_vec  = ['solid','dotted','dashed','dashdot','solid','dotted','dashed','dashdot','solid','dotted','dashed','dashdot']
    
    for i in range(len(pig_keys)):
    
        pig_name =  pig_keys[i]
        p_threshold = np.percentile(df_hplc_combined[pig_name]/df_hplc_combined['hplc_Tot_Chl_a'], 90)
        threshold_index = np.where(df_hplc_combined[pig_name]/df_hplc_combined['hplc_Tot_Chl_a'] > p_threshold)
        df_ap_threshold = df_ap_combined.iloc[threshold_index]
        ap_threshold = np.array(df_ap_threshold)
        
        plt.plot(wv,np.nanmedian(ap_threshold,axis=0),label=labels[i],linestyle=style_vec[i],linewidth=2)
        df_ap_threshold = []
        ap_threshold =[]
        
    #plt.legend()
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('$a_{p}$ m$^{-1}$')
    plt.xlim(400,720)
    
    plt.subplot(2,2,2)
    
    plt.plot(wv,global_med,label='Global median', color='black')    
    for i in range(len(pig_keys)):
    
        pig_name =  pig_keys[i]
        p_threshold = np.percentile(df_hplc_combined[pig_name]/df_hplc_combined['hplc_Tot_Chl_a'], 90)
        threshold_index = np.where(df_hplc_combined[pig_name]/df_hplc_combined['hplc_Tot_Chl_a'] > p_threshold)
        df_ap_threshold = df_ap_int_combined.iloc[threshold_index]
        ap_threshold = np.array(df_ap_threshold)
    
        plt.plot(wv,np.nanmedian(ap_threshold,axis=0),label=labels[i],linestyle=style_vec[i],linewidth=2)
        df_ap_threshold = []
        ap_threshold =[]
        
    plt.legend(fontsize=14)
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('$<a_{p}>$')
    plt.xlim(400,720)
    plt.grid()
    
    
    plt.subplot(2,2,3)
    for i in range(len(pig_keys)):
    
        pig_name =  pig_keys[i]
        p_threshold = np.percentile(df_hplc_combined[pig_name]/df_hplc_combined['hplc_Tot_Chl_a'], 90)
        threshold_index = np.where(df_hplc_combined[pig_name]/df_hplc_combined['hplc_Tot_Chl_a'] > p_threshold)
        df_ap_threshold = df_ap_int_combined.iloc[threshold_index]
        ap_threshold = np.array(df_ap_threshold)
        ap_resid = np.nan*np.ones([len(ap_threshold), len(ap_threshold.T)])
        for j in range(len(ap_threshold)):
              ap_resid[j,:] = (ap_threshold[j,:] - global_med)
        plt.plot(wv,np.nanmedian(ap_resid,axis=0),label=labels[i],linestyle=style_vec[i],linewidth=2)
        
        df_ap_threshold = []
        df_threshold =[]
        ap_resid[j,:] 
    
    #plt.legend()
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('$<a_{p}>$ residual')
    plt.xlim(400,720)
    plt.tight_layout()
    plt.grid()
    
    
    plt.subplot(2,2,4)
    for i in range(len(pig_keys)):
    
        pig_name =  pig_keys[i]
        p_threshold = np.percentile(df_hplc_combined[pig_name]/df_hplc_combined['hplc_Tot_Chl_a'], 90)
        threshold_index = np.where(df_hplc_combined[pig_name]/df_hplc_combined['hplc_Tot_Chl_a'] > p_threshold)
        df_ap_threshold = df_ap_int_combined.iloc[threshold_index]
        ap_threshold = np.array(df_ap_threshold)
        ap_resid = np.nan*np.ones([len(ap_threshold), len(ap_threshold.T)])
        for j in range(len(ap_threshold)):
              ap_resid[j,:] = 100*(ap_threshold[j,:] - global_med)/global_med
        plt.plot(wv,np.nanmedian(ap_resid,axis=0),label=labels[i],linestyle=style_vec[i],linewidth=2)
        
        df_ap_threshold = []
        df_threshold =[]
        ap_resid[j,:] 
    
    #plt.legend()
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('$<a_{p}>$ percentage residual [%]')
    plt.ylim(-40,40)
    plt.xlim(400,720)
    plt.tight_layout()
    plt.grid()

    

    return


if __name__ == '__main__':
    
    ###########################################################################
    # AMT 19 - seabass and nc data (for cross-check)
    fn_nc_19 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT19/Processed/Step3/amt19_final_with_debiased_chl.nc'
   # fn_chl_19 = '/data/datasets/cruise_data/active/ACS_Chl/amtacs/AMT19_ACS_CHL-A_MEDFILT_BIAS_CORRECTED.csv'
  #  data_chl_19 = pd.read_csv(fn_chl_19)
    data_nc_19 = xr.open_dataset(fn_nc_19)
    print(sorted(list(data_nc_19.keys())))

    
    ###########################################################################
    # AMT 22 - seabass and nc data (for cross-check)
    fn_nc_22 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT22/Processed/Step3/amt22_final_with_debiased_chl.nc'
    fn_chl_22 = '/data/datasets/cruise_data/active/ACS_Chl/amtacs/AMT22_ACS_CHL-A_MEDFILT_BIAS_CORRECTED.csv'
   # data_chl_22 = pd.read_csv(fn_chl_22)
    data_nc_22 = xr.open_dataset(fn_nc_22)
    print(sorted(list(data_nc_22.keys())))
    
    
    ############################################################################
    # AMT 23 - seabass and nc data (for cross-check)
    fn_nc_23 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT23/Processed/underway/Step3/amt23_final_with_debiased_chl.nc'
  #  fn_chl_23 = '/data/datasets/cruise_data/active/ACS_Chl/amtacs/AMT23_ACS_CHL-A_MEDFILT_BIAS_CORRECTED.csv'
    data_nc_23 = xr.open_dataset(fn_nc_23)
    print(sorted(list(data_nc_23.keys()))) 
    
    
    ###########################################################################
    # AMT 24 - seabass and nc data (for cross-check)
    fn_nc_24 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT24/Processed/Uway/Step3/amt24_final_with_debiased_chl.nc'
   # fn_chl_24 = '/data/datasets/cruise_data/active/ACS_Chl/amtacs/AMT25_ACS_CHL-A_MEDFILT_BIAS_CORRECTED.csv'
    data_nc_24 = xr.open_dataset(fn_nc_24)
    print(sorted(list(data_nc_24.keys()))) 
    
      
    ###########################################################################
    # AMT 25 - seabass and nc data (for cross-check)
    fn_nc_25 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT25/Processed/UWay/Step3/amt25_final_with_debiased_chl.nc'
# fn_chl_25 = '/data/datasets/cruise_data/active/ACS_Chl/amtacs/AMT25_ACS_CHL-A_MEDFILT_BIAS_CORRECTED.csv'
    data_nc_25 = xr.open_dataset(fn_nc_25)
   # data_chl_25 = pd.read_csv(fn_chl_25)
    
    # output initial plots 
    #  plot_median_ap_cp(data_nc_25)
    #  med_unc_plot(data_nc_25)
   # chl_check(data_nc_25,data_chl_25)
    

    ###########################################################################
    # AMT 26 - seabass and nc data (for cross-check)
   # fn_mask_26 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/LH_mask_AMT26.csv'
   # mask_26 = pd.read_csv(fn_mask_26)
   # plot_median_ap_cp(data_nc_25)

    #pigment_check(data_nc_26, data_hplc_25)
 #   med_unc_plot(data_nc_25)
  #  chl_check(data_nc_25,data_chl_25)
    
    
          
    ###########################################################################
    # AMT 26 - seabass and nc data (for cross-check)
   # dir_sb_26 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT26/Source/SeaBASS_submit/sb_processed/'
  #  fn_sb_acs_26 = dir_sb_26 + 'AMT26_InLine0_ACS_20160923_20161102_Particulate_v20230622.sb'
  #  fn_sb_hplc_26 = dir_sb_26 + 'AMT26_HPLC_20160923_20161102_v20230622.sb'
       
    fn_nc_26 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT26/Processed/Underway/Step3/amt26_final_with_debiased_chl.nc'
  #  fn_chl_26 = '/data/datasets/cruise_data/active/ACS_Chl/amtacs/AMT26_ACS_CHL-A_MEDFILT_BIAS_CORRECTED.csv'
       
    # load acs sb
   # data_acs_26 = sbs.readSB(fn_sb_acs_26)
   # print(data_acs_26.data.keys())
   # data_acs_26 = data_acs_26.data
    
     # load hplc sb
 #   data_hplc_26 = sbs.readSB(fn_sb_hplc_26) 
 #   print(data_hplc_26.data.keys())
 #   data_hplc_26 = data_hplc_26.data
    
    # load netcdf
    data_nc_26 = xr.open_dataset(fn_nc_26)
   # print(sorted(list(data_nc_26.keys()))) 
    
#    # load previously-processed chl
 #   data_chl_26 = pd.read_csv(fn_chl_26)
    
    # load longhurst mask
#    fn_mask_26 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/LH_mask_AMT26.csv'
 #   mask_26 = pd.read_csv(fn_mask_26)
        
    # output initial plots 
    # plot_median_ap_cp(data_nc_26)
    # plot_median_ap_province(data_nc_26, mask_26['LH_Province'])
    # pigment_check(data_nc_26, data_hplc_26)
    # med_unc_plot(data_nc_26)
    #chl_check(data_nc_26,data_chl_26)
       

      
    ###########################################################################
    # AMT 27 - seabass and nc data (for cross-check)
   # dir_sb_27 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT27/SeaBASS_submit/sb_processed/'
  #  fn_sb_acs_27 = dir_sb_27 + 'AMT27_InLine_ACS_20170924_20171101_Particulate_v20230622.sb'
  #  fn_sb_hplc_27 = dir_sb_27 + 'AMT27_HPLC_20170924_20171101_v20230622.sb'
    fn_nc_27 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT27/Processed/Underway/Step3/amt27_final_with_debiased_chl.nc'
   # fn_chl_27 = '/data/datasets/cruise_data/active/ACS_Chl/amtacs/AMT27_ACS_CHL-A_MEDFILT_BIAS_CORRECTED.csv'
      
    
    # load acs sb
   # data_acs_27 = sbs.readSB(fn_sb_acs_27)
  #  print(data_acs_27.data.keys())
  #  data_acs_27 = data_acs_27.data
    
     # load hplc sb
  #  data_hplc_27 = sbs.readSB(fn_sb_hplc_27) 
   # print(data_hplc_27.data.keys())
   # data_hplc_27 = data_hplc_27.data
    
    # load netcdf
    data_nc_27 = xr.open_dataset(fn_nc_27)
    print(sorted(list(data_nc_27.keys())))   
          
    # load longhurst mask
   #fn_mask_27 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/LH_mask_AMT27.csv'
   # mask_27 = pd.read_csv(fn_mask_27)
          
   # # load previously-processed chl
   # data_chl_27 = pd.read_csv(fn_chl_27)
        
    # output initial plots 
    # plot_median_ap()
    # plot_median_ap_cp(data_nc_27)

    # pigment_check(data_nc_27, data_hplc_27)
    # med_unc_plot(data_nc_27)
    # chl_check(data_nc_27, data_chl_27)
    
      
    ###########################################################################
    # AMT 28 - seabass and nc data (for cross-check)     
   # dir_sb_28 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT28/Source/SeaBASS_submit/sb_processed/'
    #fn_sb_acs_28 =  dir_sb_28 + 'AMT28_InLine0_ACS_20180925_20181028_Particulate_v20230622.sb'
   # fn_sb_hplc_28 =  dir_sb_28 + 'AMT28_HPLC_20180925_20181027_v20230622.sb'
    fn_nc_28 = '/users/rsg/tjor/scratch_network/AMT_underway/AMT28/Processed/Underway/Step3/amt28_final_with_debiased_chl.nc'
   # fn_chl_28 = '/data/datasets/cruise_data/active/ACS_Chl/amtacs/AMT28_ACS_CHL-A_MEDFILT_BIAS_CORRECTED.csv'
    
    # load acs sb
  #  data_acs_28 = sbs.readSB(fn_sb_acs_28)
  #  print(data_acs_28.data.keys())
  #  data_acs_28 = data_acs_28.data
    
     # load hplc sb
  #  data_hplc_28 = sbs.readSB(fn_sb_hplc_28) 
  #  print(data_hplc_28.data.keys())
  #  data_hplc_28 = data_hplc_28.data
    
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


    # LH29 = longhurst_mask(data_nc_29['uway_lat'], data_nc_29['uway_lon'], 'AMT29')
    # LH27 = longhurst_mask(data_nc_27['uway_lat'], data_nc_27['uway_lon'], 'AMT27')
    # LH25 = longhurst_mask(data_nc_25['uway_lat'], data_nc_25['uway_lon'], 'AMT25')
    # LH24 = longhurst_mask(data_nc_24['uway_lat'], data_nc_24['uway_lon'], 'AMT24')
    # LH23 = longhurst_mask(data_nc_23['uway_lat'], data_nc_24['uway_lon'], 'AMT23')
    # LH22 = longhurst_mask(data_nc_22['uway_lat'], data_nc_24['uway_lon'], 'AMT22')
    # LH19 = longhurst_mask(data_nc_19['uway_lat'], data_nc_24['uway_lon'], 'AMT19')
     
    
    #  LH29 = longhurst_mask_HPLC(data_nc_29['hplc_lat'], data_nc_29['hplc_lon'], 'AMT29')
 
    # keys = list(data_nc.data_vars)
    # for i in range(len(keys)):
    #     if "hplc" in keys()[i]:
    
  #  _9cruise_LH_mask(field = 'LH_Province')

    # Load Longhurst masks
    
    fn_mask_hplc = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/LH_mask_HPLC_AMTAMT.csv'
    LH_HPLC = pd.read_csv(fn_mask_hplc)
    

    
    
    fn_mask_19 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/LH_mask_AMT19.csv'
    mask_19 = pd.read_csv(fn_mask_19)
    
    fn_mask_22 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/LH_mask_AMT22.csv'
    mask_22 = pd.read_csv(fn_mask_22)
    
    fn_mask_23 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/LH_mask_AMT23.csv'
    mask_23 = pd.read_csv(fn_mask_23)
    
    fn_mask_24 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/LH_mask_AMT24.csv'
    mask_24 = pd.read_csv(fn_mask_24)
    
    fn_mask_25 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/LH_mask_AMT25.csv'
    mask_25 = pd.read_csv(fn_mask_25)
    
    fn_mask_26 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/LH_mask_AMT26.csv'
    mask_26 = pd.read_csv(fn_mask_26)
        
    fn_mask_27 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/LH_mask_AMT27.csv'
    mask_27 = pd.read_csv(fn_mask_27)
        
    fn_mask_28 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/LH_mask_AMT28.csv'
    mask_28 = pd.read_csv(fn_mask_28)
    
    fn_mask_29 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/LH_mask_AMT29.csv'
    mask_29 = pd.read_csv(fn_mask_29)


    pigment_ref = 'hplc_Tot_Chl_a'
    pigment_ref_label= 'Tot Chl-a'


