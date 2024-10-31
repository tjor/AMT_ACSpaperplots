#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script is for plotting functions for the manuscript:

Jordan et al. 2024,A compilation of surface inherent optical properties and 
phytoplankton pigment concentrations 
from the Atlantic Meridional Transect, submitted to ESSD in 2024.

It uses netcdf versions of the data files in the data release paper. The script 
checkSeaBASS_ESSD.py is used to check equivalance of data fields in netcdf and 
SeaBASS.

@author: tom jordan - tjor@pml.ac.uk
"""

# general imports
import pandas as pd
import numpy as np
import xarray as xr
import datetime
import scipy
from scipy import signal as sg


# matplotlib/cartopy
import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.colors as cl
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# use for longhurst
from xml.dom.minidom import *

# General functions
def prcrng(x):
    'robust std'
    return (np.nanpercentile(x,84) - np.nanpercentile(x,16))/2


# Longhust functions
def longhurst_mask_HPLC(lat, lon, cruise): 
    
    province_mask = []
    for i in range(len(lat)):
      print(i)
      province_mask.append(longhurst(lat[i], lon[i]))
      
    df_province = pd.DataFrame(province_mask, columns=['LH_Province'])
    df_province.to_csv('LH_mask_HPLC_AMT.csv')

    return  df_province
        

def longhurst_mask(lat, lon, cruise): 
    'used to generate LH masks. These are saved as csvs and called in the data analysis'
    
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

   
# Coverage map functions
def plot_coverage():
    'Fig 1A in MS'
    
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
                                    ccrs.PlateCarree(),facecolor='none',edgecolor='bisque',linestyle='--')
    ax.add_feature(shape_feature)
    
    colors = cm.tab10(np.linspace(0,1,10))    
       
    ax.scatter(data_nc_19['uway_lon'], data_nc_19['uway_lat'], s=0.5, color='white', label='AMT 19: 2009')
    ax.scatter(data_nc_22['uway_lon'], data_nc_22['uway_lat'], s=0.5, color=colors[1], label='AMT 22: 2012')
    ax.scatter(data_nc_23['uway_lon'], data_nc_23['uway_lat'], s=0.5, color=colors[2], label='AMT 23: 2013')
    ax.scatter(data_nc_24['uway_lon'], data_nc_24['uway_lat'], s=0.5, color=colors[3], label='AMT 24: 2014')
    ax.scatter(data_nc_25['uway_lon'], data_nc_25['uway_lat'], s=0.5, color=colors[4], label='AMT 25: 2015')
    ax.scatter(data_nc_26['uway_lon'], data_nc_26['uway_lat'], s=0.5, color=colors[5], label='AMT 26: 2016')
    ax.scatter(data_nc_27['uway_lon'], data_nc_27['uway_lat'], s=0.5, color=colors[6], label='AMT 27: 2017')
    ax.scatter(data_nc_28['uway_lon'], data_nc_28['uway_lat'], s=0.5, color=colors[7], label='AMT 28: 2018')
    ax.scatter(data_nc_29['uway_lon'], data_nc_29['uway_lat'], s=0.5, color=colors[8], label='AMT 29: 2019')
    plt.legend(markerscale=10 ,loc=4,fontsize=14)                
    
    plt.text(.33, .92,   'NADR', ha='left', va='top', color='white', transform=ax.transAxes,fontsize=16)
    plt.text(.12, .80,   'NASW', ha='left', va='top', color='white', transform=ax.transAxes,fontsize=16)
    plt.text(.53, .78,   'NASE', ha='left', va='top', color='white', transform=ax.transAxes,fontsize=16)
    plt.text(.16, .65,   'NATR', ha='left', va='top', color='white', transform=ax.transAxes,fontsize=16)
    plt.text(.26, .53,   'WTRA', ha='left', va='top', color='white', transform=ax.transAxes,fontsize=16)
    plt.text(.50, .28,   'SATL', ha='left', va='top', color='white', transform=ax.transAxes,fontsize=16)
    plt.text(.50, .13,   'SSTC', ha='left', va='top', color='white', transform=ax.transAxes,fontsize=16)
    plt.text(.42, .06,   'SANT', ha='left', va='top', color='white', transform=ax.transAxes,fontsize=16)    
  
    filename  =  fig_dir + '/'  + '_AMT_coverage.png'
    plt.savefig(filename,dpi=600)    

    return


    
def _relative_ACS_time(data_nc):
    'subroutine in Fig 1B'

    timestamp = []
    for i in range(0,len(data_nc['time'])):
        str_i = '2023'+ str(data_nc['time'].values[i])[4:19] # create fake date all in same year
        timestamp.append(datetime.datetime.strptime(str_i, '%Y-%m-%dT%H:%M:%S'))
    timestamp = np.array(timestamp)

    return timestamp


    
def _relative_HPLC_time(data_nc):
    'subroutine in Fig 1B'
    
    timestamp = []
    for i in range(0,len(data_nc['hplc_time'])):
        str_i = '2023'+ str(data_nc['hplc_time'].values[i])[4:19] # create fake date all in same year 
        timestamp.append(datetime.datetime.strptime(str_i, '%Y-%m-%dT%H:%M:%S'))
    timestamp = np.array(timestamp)
    timestamp = timestamp[data_nc['hplc_depth'].values < 10]  #

    return timestamp


def _color_by_prov(ts, tsh, mask, data_nc, plot_index, two_systems=False):
    'Sub-routine used in Figure 1B'
    
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
                plt.plot_date(plot_index*np.array(1*(~np.isnan(data_nc['acx_chl_debiased'])))[mask==prov[i]],ts[mask==prov[i]], xdate=False, marker="s", ms=10,color=colors[i], label = 'IOP: ' + str(prov[i]))
    else:
       for i in range(len(prov)):
           plt.plot_date(plot_index*np.array(1*(~np.isnan(data_nc['acs_chl_debiased'])))[mask==prov[i]], ts[mask==prov[i]], xdate=False, marker="s",ms=10, color=colors[i])
           if plot_index ==9:
               plt.plot_date(plot_index*np.array(1*(~np.isnan(data_nc['acs_chl_debiased'])))[mask==prov[i]], ts[mask==prov[i]], xdate=False, marker="s", ms=10, color=colors[i], label = 'IOP: ' + str(prov[i]))
           if two_systems ==True:
               combined = np.array(1*(~np.isnan(data_nc['acs2_chl_debiased'])) + 1*(~np.isnan(data_nc['acs_chl_debiased'])))
               combined[combined>0]=1
               plt.plot_date(plot_index*combined[mask==prov[i]],ts[mask==prov[i]], xdate=False, ms=10, color=colors[i])
    

    plt.plot_date(plot_index*np.ones(len(tsh)), tsh, xdate=False, ms=3, color='black')

    if plot_index ==9:
        plt.plot_date(plot_index*np.array(1*(~np.isnan(data_nc['acs_chl_debiased'])))[mask==prov_extra[0]],ts[mask==prov_extra[0]], marker="s", xdate=False, ms=12,color='gray',label='IOP: OTHER PROVINCES')
        plt.plot_date(plot_index*np.ones(len(tsh)),tsh, xdate=False, ms=3, color='black', label='HPLC')

        plt.legend(fontsize=15)
        
    return 


def plot_AMT_timeline():
    'Figure 1B/timeline'
    
    plt.figure(figsize =(9,14))
    plt.rcParams.update({'font.size': 16})
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
    
       
    filename  =  fig_dir + '/'  + '_AMT_timeline.png'
    plt.savefig(filename,dpi=600,bbox_inches='tight')

    return


# IOP plot functions
def _9cruise_IOPsum(field, field2):
    'sacks 2D arrays for each IOP variable for 9 cruises. field1 and field2 account for acs and acs2'
    
    combined_array = np.concatenate([np.array(data_nc_19[field].values), data_nc_22[field].values, data_nc_23[field].values, 
                                     data_nc_24[field].values,  data_nc_24[field2].values, 
                                     data_nc_25[field].values,  data_nc_25[field2].values, 
                                     data_nc_26[field].values, data_nc_27[field].values, 
                                     data_nc_28[field].values, data_nc_29[field].values], axis=0)
                    
    return combined_array           


def _9cruise_chlsum():
    'stacks chl transects for 9 cruises'
    
    field = 'acs_chl_debiased'
    field2= 'acs2_chl_debiased'
    field3 = 'acx_chl_debiased'
    
    combined_array = np.concatenate([np.array(data_nc_19[field3].values), data_nc_22[field].values, data_nc_23[field].values, 
                                     data_nc_24[field].values,  data_nc_24[field2].values, 
                                     data_nc_25[field].values,  data_nc_25[field2].values, 
                                     data_nc_26[field].values, data_nc_27[field].values, 
                                     data_nc_28[field3].values, data_nc_29[field].values], axis=0)
                    
    return combined_array      



def _9cruise_LH_mask(field = 'LH_Province'):
    'stacks LH mask'
    
    field = 'LH_Province'
    combined_LH = np.concatenate([mask_19[field], mask_22[field], mask_23[field], 
                                         mask_24[field],  mask_24[field],   # double-count for AMT 25 and 24 due to 2 AC-S systems'
                                         mask_25[field], mask_25[field], mask_26[field], 
                                         mask_27[field], mask_28[field], mask_29[field]], axis=0)
    return combined_LH


def int_norm(data): 
    'integral normalization of IOP spectra, (follows Boss 2013)' 
    data_in = np.nan*np.ones([len(data),len(data.T)])
    for i in range(len(data)):
        data_in[i,:] = data[i,:]/np.mean(data[i,:])
        
    return data_in


def int_globalmed():
    'Find global median integral (as used in Sect 3.5)'
    ap_mat = _9cruise_IOPsum('acs_ap', 'acs2_ap')       
    global_med = np.nanmedian(ap_mat,axis=0)
    ap_int = int_norm(ap_mat)
    global_med_int = np.nanmedian(ap_int,axis=0)

    return  global_med, global_med_int



def plot_median_ap_province():
        'plots median ap, bp, cp by province'
                
        mask = _9cruise_LH_mask(field = 'LH_Province')
  
        ap_mat = _9cruise_IOPsum('acs_ap', 'acs2_ap')
        bp_mat = _9cruise_IOPsum('acs_bp', 'acs2_bp')
        cp_mat = _9cruise_IOPsum('acs_cp', 'acs2_cp')
          
        ap_int = int_norm(ap_mat)
        bp_int = int_norm(bp_mat)
        cp_int = int_norm(cp_mat)
        
        chl_tot = _9cruise_chlsum()
        
        C_NADR = np.round(1000*np.nanmedian(chl_tot[mask=='NADR']))/1000
        C_NASE = np.round(1000*np.nanmedian(chl_tot[mask=='NASE']))/1000
        C_NASW = np.round(1000*np.nanmedian(chl_tot[mask=='NASW']))/1000
        C_NATR = np.round(1000*np.nanmedian(chl_tot[mask=='NATR']))/1000 
        C_WTRA = np.round(1000*np.nanmedian(chl_tot[mask=='WTRA']))/1000
        C_SATL = np.round(1000*np.nanmedian(chl_tot[mask=='SATL']))/1000
        C_SSTC = np.round(1000*np.nanmedian(chl_tot[mask=='SSTC']))/1000
        C_SANT = np.round(1000*np.nanmedian(chl_tot[mask=='SANT']))/1000
        
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
        
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_mat[mask=='NADR',:],axis=0),color=colors[0])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_mat[mask=='NASW',:],axis=0),color=colors[1])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_mat[mask=='NASE',:],axis=0),color=colors[2])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_mat[mask=='NATR',:],axis=0),color=colors[3])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_mat[mask=='WTRA',:],axis=0),color=colors[4])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_mat[mask=='SATL',:],axis=0),color=colors[5])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_mat[mask=='SSTC',:],axis=0),color=colors[6])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_mat[mask=='SANT',:],axis=0),color=colors[7])
        
        #plt.legend(loc=(1.04, 0))
        plt.xlim(400,720)
        plt.xlabel('Wavelength, $\lambda$ [nm]')
        plt.ylabel('$a_{p}(\lambda)$  [m$^{-1}$]')
        ax = plt.gca()
        ax.set_ylim(bottom=0)
        plt.text(.86, .95,  'A', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
                
        
        plt.subplot(2,3,2)
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_mat[mask=='NADR',:],axis=0),color=colors[0])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_mat[mask=='NASW',:],axis=0),color=colors[1])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_mat[mask=='NASE',:],axis=0),color=colors[2])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_mat[mask=='NATR',:],axis=0),color=colors[3])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_mat[mask=='WTRA',:],axis=0),color=colors[4])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_mat[mask=='SATL',:],axis=0),color=colors[5])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_mat[mask=='SSTC',:],axis=0),color=colors[6])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_mat[mask=='SANT',:],axis=0),color=colors[7])
        
        #plt.legend(loc=(1.04, 0))
        plt.xlim(400,720)
        plt.xlabel('Wavelength, $\lambda$ [nm]')
        plt.ylabel('$b^{\prime}_{p}(\lambda)$  [m$^{-1}$]')
        ax = plt.gca()
        ax.set_ylim(bottom=0)
        plt.text(.86, .95,   'B', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
        

        plt.subplot(2,3,3)   
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_mat[mask=='NADR',:],axis=0),label='NADR: N = ' + str(N_NADR) + ', Median $C_{a}(a_{p})$ = ' + str(C_NADR) + ' mg m$^{-3}$', color=colors[0])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_mat[mask=='NASW',:],axis=0),label='NASW: N = ' + str(N_NASW) + ', Median $C_{a}(a_{p})$ = ' + str(C_NASW) + ' mg m$^{-3}$', color=colors[1])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_mat[mask=='NASE',:],axis=0),label='NASE: N = ' + str(N_NASE) + ', Median $C_{a}(a_{p})$ = ' + str(C_NASE) + ' mg m$^{-3}$', color=colors[2])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_mat[mask=='NATR',:],axis=0),label='NATR: N = ' + str(N_NATR) + ', Median $C_{a}(a_{p})$ = ' + str(C_NATR) + ' mg m$^{-3}$', color=colors[3])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_mat[mask=='WTRA',:],axis=0),label='WTRA: N = ' + str(N_WTRA) + ', Median $C_{a}(a_{p})$ = ' + str(C_WTRA) + ' mg m$^{-3}$', color=colors[4])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_mat[mask=='SATL',:],axis=0),label='SATL: N = ' + str(N_SATL) + ', Median $C_{a}(a_{p})$ = ' + str(C_SATL) + ' mg m$^{-3}$', color=colors[5])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_mat[mask=='SSTC',:],axis=0),label='SSTC: N = ' + str(N_SSTC) + ', Median $C_{a}(a_{p})$ = ' + str(C_SSTC) + ' mg m$^{-3}$', color=colors[6])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_mat[mask=='SANT',:],axis=0),label='SANT: N = ' + str(N_SANT) + ', Median $C_{a}(a_{p})$ = ' + str(C_SANT) + ' mg m$^{-3}$', color=colors[7])
        plt.xlim(400,720)
        plt.xlabel('Wavelength, $\lambda$ [nm]')
        plt.ylabel('$c_{p}(\lambda)$  [m$^{-1}$]')
        ax = plt.gca()
        ax.set_ylim(bottom=0)
        plt.text(.86, .95,   'C', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
        
        plt.subplot(2,3,4)
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_int[mask=='NADR',:],axis=0),color=colors[0])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_int[mask=='NASW',:],axis=0),color=colors[1])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_int[mask=='NASE',:],axis=0),color=colors[2])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_int[mask=='NATR',:],axis=0),color=colors[3])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_int[mask=='WTRA',:],axis=0),color=colors[4])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_int[mask=='SATL',:],axis=0),color=colors[5])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_int[mask=='SSTC',:],axis=0),color=colors[6])
        plt.plot(data_nc['acs_wv'],np.nanmedian(ap_int[mask=='SANT',:],axis=0),color=colors[7])

        plt.xlim(400,720)
        plt.xlabel('Wavelength, $\lambda$ [nm]')
        plt.ylabel('$<a_{p}(\lambda)>$')
        ax = plt.gca()
        ax.set_ylim(bottom=0)
        plt.text(.86, .95, 'D', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
        
        
        plt.subplot(2,3,5)
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_int[mask=='NADR',:],axis=0),color=colors[0])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_int[mask=='NASW',:],axis=0),color=colors[1])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_int[mask=='NASE',:],axis=0),color=colors[2])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_int[mask=='NATR',:],axis=0),color=colors[3])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_int[mask=='WTRA',:],axis=0),color=colors[4])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_int[mask=='SATL',:],axis=0),color=colors[5])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_int[mask=='SSTC',:],axis=0),color=colors[6])
        plt.plot(data_nc['acs_wv'],np.nanmedian(bp_int[mask=='SANT',:],axis=0),color=colors[7])
        plt.plot(data_nc['acs_wv'], one_over_lamba, color='black', linestyle = 'dashed')

        #plt.legend(loc=(1.04, 0))
        plt.xlim(400,720)
        plt.xlabel('Wavelength, $\lambda$ [nm]')
        plt.ylabel('$<b^{\prime}_{p}(\lambda) >$')
        ax = plt.gca()
        plt.ylim(0.7,1.5)
        # ax.set_ylim(bottom=0)
        plt.text(.86, .95,  'E', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
        
        plt.subplot(2,3,6)   
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_int[mask=='NADR',:],axis=0),color=colors[0])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_int[mask=='NASW',:],axis=0),color=colors[1])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_int[mask=='NASE',:],axis=0),color=colors[2])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_int[mask=='NATR',:],axis=0),color=colors[3])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_int[mask=='WTRA',:],axis=0),color=colors[4])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_int[mask=='SATL',:],axis=0),color=colors[5])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_int[mask=='SSTC',:],axis=0),color=colors[6])
        plt.plot(data_nc['acs_wv'],np.nanmedian(cp_int[mask=='SANT',:],axis=0),color=colors[7])
        plt.plot(data_nc['acs_wv'], one_over_lamba, color='black', linestyle = 'dashed', label = '$\lambda^{-1}$')
        # plt.legend(loc=(1.04, 0))
        plt.xlim(400,720)
        plt.xlabel('Wavelength, $\lambda$ [nm]')
        plt.ylabel('$<c_{p}(\lambda)>$')
        ax = plt.gca()
        plt.ylim(0.7,1.5)
        # ax.set_ylim(bottom=0)
        plt.text(.86, .95,   'F', ha='left', va='top', transform=ax.transAxes,fontsize=24)
        
        plt.figlegend(loc='upper center', bbox_to_anchor=(0.5, -0.02),
          fancybox=True, ncol=2, fontsize=14)
        

        plt.tight_layout()
        
        filename  =  fig_dir + '/'  + '_provincemedians.png'
        plt.savefig(filename,dpi=600,bbox_inches='tight')

        return

    
def plot_uncertainties():
    'Uncertainty plot for paper'
        
    apu_mat = _9cruise_IOPsum('acs_ap_u', 'acs2_ap_u')
    bpu_mat = _9cruise_IOPsum('acs_bp_u', 'acs2_bp_u')
    cpu_mat = _9cruise_IOPsum('acs_cp_u', 'acs2_cp_u')
    
    ap_mat = _9cruise_IOPsum('acs_ap', 'acs2_ap')
    bp_mat = _9cruise_IOPsum('acs_bp', 'acs2_bp')
    cp_mat = _9cruise_IOPsum('acs_cp', 'acs2_cp')
        
    colors = cm.tab10(np.linspace(0,1,10)) 
    c1 = colors[0]
    c2 = colors[1]
    c3 = colors[3]
    
    # plt.figure() # ratio figure for reviewer comments
    # plt.rcParams.update({'font.size': 16})
    # plt.plot(data_nc_29['acs_wv'],np.nanmedian(cpu_mat,axis=0)/np.nanmedian(apu_mat,axis=0),label='Median', color=c1)
    # plt.ylim(0.8,1.6)
    # plt.xlabel('Wavelength, $\lambda$ [nm]')
    # plt.ylabel('Median $\sigma_{c_{p}(\lambda)}/\sigma_{a_{p}(\lambda)}$ ratio')
    # ax = plt.gca()
    # plt.text(.05, .95,  'D', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
    # plt.legend(fontsize=12)
    
    plt.figure(figsize=(15,9))
    plt.subplot(2,3,1)
    plt.rcParams.update({'font.size': 16})
    plt.plot(data_nc_29['acs_wv'],np.nanmedian(apu_mat,axis=0), label='Median', color=c1)
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(apu_mat,25,axis=0), color=c1, linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(apu_mat,75,axis=0), label='Quartiles', color=c1, linestyle='dashed')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(apu_mat,90,axis=0), label='10$^{th}$ & 90$^{th}$ percentiles', color=c1,linestyle='dotted')
    plt.plot(data_nc_29['acs_wv'],np.nanpercentile(apu_mat,10,axis=0), color=c1,linestyle='dotted')
    plt.rcParams.update({'font.size': 16})
    plt.xlim(400, 720)
    plt.ylim(0, 0.015)
    plt.xlabel('Wavelength, $\lambda$ [nm]')
    plt.ylabel('$\sigma_{a_p(\lambda)}$ [m$^{-1}$]')
    plt.rcParams.update({'font.size': 16})
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
    plt.xlabel('Wavelength, $\lambda$ [nm]')
    #plt.ylabel('Uncertainty [m$^{-1}$]')
    plt.ylabel('$\sigma_{b^{\prime}_{p}(\lambda)}$ [m$^{-1}$]')
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
    plt.ylabel('$\sigma_{c_p(\lambda)}$  [m$^{-1}$]')
    plt.xlabel('Wavelength, $\lambda$ [nm]')
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
    plt.xlabel('Wavelength, $\lambda$ [nm]')
    plt.ylabel('$\sigma_{a_{p}(\lambda)}$ [%]')
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
    plt.xlabel('Wavelength, $\lambda$ [nm]')
    ax = plt.gca()
    plt.ylabel('$\sigma_{b^{\prime}_{p}(\lambda)}$  [%]')
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
    plt.xlabel('Wavelength, $\lambda$ [nm]')
    plt.ylabel('$\sigma_{c_{p}(\lambda)}$  [%]')
    ax = plt.gca()
    plt.text(.05, .95,  'F', ha='left', va='top', transform=ax.transAxes,fontsize=24) 
    #plt.ylabel('Percentage uncertainty [%]')
    #plt.legend(fontsize=12)
    
    plt.tight_layout(pad=1.6)
    
    filename  =  fig_dir + '/'  + '_uncertainties.png'
    plt.savefig(filename,dpi=600)
    
    return


def plot_676_443():
      'Red-blue ratio plot'
 
      plt.figure(figsize=(13,5))
      plt.rcParams.update({'font.size': 16})
      
      plt.subplot(1,2,2)
      mask = _9cruise_LH_mask(field = 'LH_Province')
      #  
      ap_mat = _9cruise_IOPsum('acs_ap', 'acs2_ap')
      ap_int = int_norm(ap_mat)
      
      colors = cm.Paired(np.linspace(0,1,10))
      
      ap_443 = (ap_int[:, 22] + ap_int[:, 23])/2
      ap_676  = ap_int[:, 138]
      ap_rat = ap_676/ap_443
      
      prov= ['NADR', 'NASW', 'NASE' ,'NATR', 'WTRA', 'SATL', 'SSTC', 'SANT'] # FKLD==1, ANTA -4, NECS = 8. do not include
      
      rat_vec =[]
      for i in range(len(prov)):
          rat_i = np.array(ap_rat)[mask ==str(prov[i])]
          rat_i =  rat_i[~np.isnan(rat_i)]
          rat_vec.append(rat_i)                                                                   
      bp = plt.boxplot(rat_vec ,showfliers=True,patch_artist=True, medianprops=dict(color='black'), whis=[10,90],widths = 0.6) 
      plt.ylim(0,0.6)
      
      plt.ylabel('$a_{p}(676)/a_{p}(443)$')
      plt.xlabel('Longhurst Province')
      
      colors = cm.Paired(np.linspace(0,1,10))
      patches = []
      for k in range(len(prov)):
        bp['boxes'][k].set_facecolor(colors[k])
    
      color_string = cl.rgb2hex(colors[k])
      patches.append(mpatches.Patch(color=color_string, label=str(prov[k])))
    
      #plt.legend(handles=patches,fontsize=14,loc=2)  
      ax = plt.gca()
      ax.set_xticklabels(labels= ['NADR', 'NASW',  'NASE' ,'NATR', 'WTRA', 'SATL', 'SSTC', 'SANT'],rotation=45)  
      # ax.set_xticks([])  
      #ax.set_yticks([])  
      # ax.set_axis_off()
      ax=plt.gca()
      plt.text(0.90, .1,  'B', ha='left', va='top', transform=ax.transAxes,fontsize=22)
      
      plt.subplot(1,2,1)
      chl = _9cruise_chlsum()
         
      # concentrtaion ratio plots
      rat_vec =[]
      rat_i = ap_rat[np.where(chl < 0.05)[0]]  
      rat_i =  rat_i[~np.isnan(rat_i)]
      rat_vec.append(rat_i)                                                                   
      rat_i = ap_rat[np.where((chl < 0.1) & (chl >0.05))[0]]
      rat_i =  rat_i[~np.isnan(rat_i)]
      rat_vec.append(rat_i)   
      rat_i = ap_rat[np.where((chl < 0.5) & (chl >0.1))[0]]
      rat_i =  rat_i[~np.isnan(rat_i)]
      rat_vec.append(rat_i)   
      rat_i = ap_rat[np.where((chl < 1) & (chl >0.5))[0]]
      rat_i =  rat_i[~np.isnan(rat_i)]
      rat_vec.append(rat_i)   
      rat_i = ap_rat[np.where((chl > 1))[0]]
      rat_i =  rat_i[~np.isnan(rat_i)]
      rat_vec.append(rat_i)   
        
      bp = plt.boxplot(rat_vec ,showfliers=True,patch_artist=True, medianprops=dict(color='black'), whis=[10, 90], widths = 0.6) 
      plt.ylim(0,0.6)
      
      colors = cm.tab20(np.linspace(0,1,10))
      patches = []
      for k in range(5):
        bp['boxes'][k].set_facecolor(colors[k])
      plt.ylabel('$a_{p}(676)/a_{p}(443)$')
      #plt.legend(handles=patches,fontsize=14,loc=2)  
      ax = plt.gca()
      ax.set_xticklabels(labels= ['$C(a_{p})$ < 0.05 ', ' 0.05 < $C(a_{p})$ < 0.1',  ' 0.1 < $C(a_{p})$ < 0.5' ,' 0.5 < $C(a_{p})$ < 1','$C(a_{p})$ > 1', ],rotation=20,fontsize=14)  
      plt.xlabel('Longhurst Province')
      
      plt.xlabel('Tot_Chl_a concentration  [mg m$^{-3}]$')

      ax=plt.gca()
      plt.text(0.90, .1,  'A', ha='left', va='top', transform=ax.transAxes,fontsize=22)
    
      plt.tight_layout()
      filename  =  fig_dir + '/'  + '_675over473.png'
      plt.savefig(filename,dpi=600,bbox_inches='tight')
      
      return


# pigment plots
def _filter_combine_pigs(data_nc,removePML_Tchlb=False):
    
    'returns a combined/filtered pigment data frame with deeper samples and replicates removed'
    
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

    df_hplc['hplc_Tot_Chl_b_noPML'] =  df_hplc['hplc_Tot_Chl_b']    
    if removePML_Tchlb==True: # this is done for hplc_Tot_Chl_b sensitvity test 
        df_hplc['hplc_Tot_Chl_b_noPML'] = np.nan*np.ones(len(df_hplc['hplc_Tot_Chl_b_noPML']))
        
    return df_hplc


def _pig_9cruises():

    df_hplc_29 = _filter_combine_pigs(data_nc_29)
    df_hplc_28 = _filter_combine_pigs(data_nc_28)
    df_hplc_27 = _filter_combine_pigs(data_nc_27)
    df_hplc_26 = _filter_combine_pigs(data_nc_26, removePML_Tchlb = False)
    df_hplc_25 = _filter_combine_pigs(data_nc_25, removePML_Tchlb = False)
    df_hplc_24 = _filter_combine_pigs(data_nc_24, removePML_Tchlb = False)
    df_hplc_23 = _filter_combine_pigs(data_nc_23, removePML_Tchlb = False)
    df_hplc_22 = _filter_combine_pigs(data_nc_22)
    df_hplc_19 = _filter_combine_pigs(data_nc_19)
    
    df_hplc_combined = pd.concat([df_hplc_29, df_hplc_28, df_hplc_27, df_hplc_26, df_hplc_25, df_hplc_24, df_hplc_23, df_hplc_22, df_hplc_19 ])
   # df_hplc_combined = pd.concat([df_hplc_29, df_hplc_28, df_hplc_27, df_hplc_22, df_hplc_19 ]) #- removing PML 
   
    #df_hplc_combined = pd.concat([df_hplc_26, df_hplc_25, df_hplc_24, df_hplc_23]) # just PML
    
    
    return  df_hplc_combined


def plot_pigs_byprov(): 
    
    ' pigment box plot'

    prov=  ['NADR', 'NASW','NASE',  'NATR', 'WTRA', 'SATL', 'SSTC', 'SANT'] # FKLD==1, ANTA -4, NECS = 8. do not include

    pig_keys = ['hplc_Tot_Chl_a','hplc_Tot_Chl_b', 'hplc_Tot_Chl_c', 
                'hplc_PPC', 'hplc_Allo','hplc_Alpha-beta-Car', 'hplc_Diadino','hplc_Diato', 'hplc_Zea',
                'hplc_PSC', 'hplc_But-fuco', 'hplc_Hex-fuco', 'hplc_Fuco', 'hplc_Perid']

    labels =  ['Tot_Chl_a', 'Tot_Chl_b', 'Tot_Chl_c' , 
               'PPC', 'Allo','alpha-beta-Car', 'Diadino', 'Diato','Zea',
               'PSC', 'But-fuco', 'Hex-fuco', 'Fuco', 'Perid']


    letters = ['A', 'B', 'C' ,'D' ,'E' ,'F','G' , 'H']


    plt.figure(figsize=(15,15))    
    plt.rcParams.update({'font.size': 18})
    for i in range(len(prov)):
        plt.subplot(3, 3, i+1)  
        plt.title(prov[i] + ': N = ' + str(np.sum(LH_HPLC['LH_Province'] ==str(prov[i]))))

        pig_vec = []
        for j in range(len(pig_keys)):
            pig_j = np.array(df_hplc[pig_keys[j]])[LH_HPLC['LH_Province'] ==str(prov[i])]   
            pig_j[pig_j<0] = np.nan
            pig_j =pig_j[~np.isnan(pig_j)]
            pig_vec.append(pig_j)    

        bp = plt.boxplot(pig_vec ,showfliers=True,patch_artist=True, medianprops=dict(color='black'), whis=[10,90],widths = 0.8) 
        plt.yscale('log')

        ax = plt.gca()
        ax.set_xticks([])  
        ax.set_ylim([0.0003, 5])
        if i==0:
            plt.ylabel('Concentration [mg m$^{-3}]$')
        if i==3:
            plt.ylabel('Concentration [mg m$^{-3}]$')
        if i==6:
            plt.ylabel('Concentration [mg m$^{-3}]$')

        plt.text(.90, .95,  letters[i], ha='left', va='top', transform=ax.transAxes,fontsize=32)  


        colors = cm.tab20(np.linspace(0,1,len(pig_keys))) 
        patches = []
        for k in range(len(pig_keys)):
            bp['boxes'][k].set_facecolor(colors[k])
            color_string = cl.rgb2hex(colors[k])
            patches.append(mpatches.Patch(color=color_string, label=labels[k]))

            plt.tight_layout(pad=1.2)

    plt.subplot(3,3,9)
    plt.legend(handles=patches,fontsize=14,loc=2,ncol=2)
    ax = plt.gca()
    ax.set_xticks([])  
    ax.set_yticks([])  
    ax.set_axis_off()

    filename  =  fig_dir + '/'  + '_pigsbyprov.png'
    plt.savefig(filename,dpi=600)


    return


def plot_pigs_byprov_ratio(): 
    
    ' pigment ratio box plot'

    prov= ['NADR', 'NASW','NASE',  'NATR', 'WTRA', 'SATL', 'SSTC', 'SANT'] # FKLD==1, ANTA -4, NECS = 8. do not include

    pig_keys = ['hplc_Tot_Chl_b', 'hplc_Tot_Chl_c', 
                'hplc_PPC', 'hplc_Allo','hplc_Alpha-beta-Car', 'hplc_Diadino', 'hplc_Diato',
                'hplc_Zea', 'hplc_PSC', 'hplc_But-fuco', 'hplc_Hex-fuco', 'hplc_Fuco', 'hplc_Perid']

    labels =  ['Tot_Chl_b', 'Tot_Chl_c' , 
               'PPC', 'Allo','alpha-beta-Car', 'Diadino', 'Diato', 'Zea',
               'PSC', 'But-fuco', 'Hex-fuco', 'Fuco', 'Perid']


    letters = ['A', 'B', 'C' ,'D' ,'E' ,'F','G' , 'H']


    plt.figure(figsize=(15,15))    
    plt.rcParams.update({'font.size': 18})
    for i in range(len(prov)):
        plt.subplot(3, 3, i+1)  
        plt.title(prov[i] + ': N = ' + str(np.sum(LH_HPLC['LH_Province'] ==str(prov[i]))))

        pig_vec = []
        for j in range(len(pig_keys)):
            pig_j = np.array(df_hplc[pig_keys[j]]/df_hplc['hplc_Tot_Chl_a'])[LH_HPLC['LH_Province'] ==str(prov[i])]
            pig_j[pig_j<0] = np.nan
            pig_j =pig_j[~np.isnan(pig_j)]
            pig_vec.append(pig_j)    

        bp = plt.boxplot(pig_vec ,showfliers=True,patch_artist=True, medianprops=dict(color='black'), whis=[10,90],widths = 0.8) 
        plt.yscale('log')


        ax = plt.gca()
        ax.set_xticks([])  
        ax.set_ylim([0.003, 3])
        if i==0:
            plt.ylabel('Concentration ratio')
        if i==3:
            plt.ylabel('Concentration ratio')
        if i==6:
            plt.ylabel('Concentration ratio')

        plt.text(.90, .95,  letters[i], ha='left', va='top', transform=ax.transAxes,fontsize=22)  

        colors = cm.tab20(np.linspace(0,1,len(pig_keys)+1)) 
        patches = []
        for k in range(len(pig_keys)):
            bp['boxes'][k].set_facecolor(colors[k+1])
            color_string = cl.rgb2hex(colors[k+1])
            patches.append(mpatches.Patch(color=color_string, label=labels[k]))

            plt.tight_layout(pad=1.2)

    plt.subplot(3,3,9)
    plt.legend(handles=patches,fontsize=14,loc=2,ncols=2)
    ax = plt.gca()
    ax.set_xticks([])  
    ax.set_yticks([])  
    ax.set_axis_off()

    filename  =  fig_dir + '/'  + '_pigsbyprovratio.png'
    plt.savefig(filename,dpi=600)

    return




def plot_pig_cov(data):
    'plot to show pigment correlation matrices'

    pig_keys = ['hplc_Tot_Chl_a',' ','hplc_Tot_Chl_b', 'hplc_Tot_Chl_b_noPML' ,'hplc_Tot_Chl_c' , 
                'hplc_PPC', 'hplc_Allo','hplc_Alpha-beta-Car', 'hplc_Diadino', 'hplc_Diato', 'hplc_Zea',
                'hplc_PSC', 'hplc_But-fuco', 'hplc_Hex-fuco', 'hplc_Fuco', 'hplc_Perid']



     #   pig_labels =  ['Tot_Chl_a', 'Tot_Chl_b', 'Tot_Chl_b$^{*}$', 'Tot_Chl_c' , 
     #              'PPC', 'Allo','alpha-beta-Car', 'hplc_Diato', 'Diadino', 'Zea',
      #             'PSC', 'But-fuco', 'Hex-fuco', 'Fuco', 'Perid']

    
    plt.figure(figsize=(18,22))
    plt.subplot(2,1,1)
    plt.rcParams.update({'font.size': 22})
    C = np.nan*np.ones([len(pig_keys), len(pig_keys)])
    for i in range(len(pig_keys)):
        for j in range(len(pig_keys)):
            if i > j:
              #  breakpoint()
                if (pig_keys[i] != ' ') & (pig_keys[j] != ' '):
                    pig_i = data[pig_keys[i]][(data[pig_keys[i]] > 0) & (data[pig_keys[j]] > 0)]
                    pig_j = data[pig_keys[j]][(data[pig_keys[i]] > 0) & (data[pig_keys[j]] > 0)]
                    C[i, j] = scipy.stats.pearsonr(pig_i, pig_j)[0]
                    if pig_keys[i] =='hplc_Tot_Chl_b' and pig_keys[j] == 'hplc_Tot_Chl_b_noPML':
                       C[i, j] = np.nan
                    if pig_keys[j] =='hplc_Tot_Chl_b' and pig_keys[i] == 'hplc_Tot_Chl_b_noPML':
                       C[i, j] = np.nan
  
          
            if i < j:     
                if (pig_keys[i] != ' ') & (pig_keys[j] != ' '):
                    pig_i = data[pig_keys[i]][(data[pig_keys[i]] > 0) & (data[pig_keys[j]] > 0)]
                    pig_j = data[pig_keys[j]][(data[pig_keys[i]] > 0) & (data[pig_keys[j]] > 0)]
                    Tchla = data['hplc_Tot_Chl_a'][(data[pig_keys[i]] > 0) & (data[pig_keys[j]] > 0)]
                    if i > 0:
                        C[i, j] = scipy.stats.pearsonr(pig_i/Tchla, pig_j/Tchla)[0]
                    elif i == 0:
                        C[i, j] = scipy.stats.pearsonr(pig_i, pig_j/Tchla)[0]
                    if pig_keys[i] =='hplc_Tot_Chl_b' and pig_keys[j] == 'hplc_Tot_Chl_b_noPML':
                        C[i, j] = np.nan
                    if pig_keys[j] =='hplc_Tot_Chl_b' and pig_keys[i] == 'hplc_Tot_Chl_b_noPML':
                        C[i, j] = np.nan
   
    plt.rcParams.update({'font.size': 22})
    plt.pcolor(np.flipud(C.T), cmap='bwr', vmin=-1, vmax=1)
    cbar = plt.colorbar()
    cbar.set_label('Correlation coefficient, $r$', rotation=90 ,labelpad=30)

    plt.xlim(0,16)
    plt.ylim(0,16)

    ax = plt.gca()
    ax.text(-0.20, 0.95,  'A', ha='left', va='top', transform=ax.transAxes,fontsize=32)  
        

    for i in range(C.shape[0]):
        for j in range(C.shape[1]):
            if np.isnan(C[i, j]) == 0:
                plt.text(i + 0.5, 15 -j + 0.5, '%.2f' % C[i, j],
                         horizontalalignment='center',
                         verticalalignment='center',
                         )
    ax = plt.gca()
    ax.set_xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5], 
    labels =  ['Tot_Chl_a', ' ', 'Tot_Chl_b','Tot_Chl_b$^{*}$' ,'Tot_Chl_c' , 
               'PPC', 'Allo','alpha-beta-Car', 'Diadino', 'Diato', 'Zea',
               'PSC', 'But-fuco', 'Hex-fuco', 'Fuco', 'Perid'])



    plt.xticks(rotation=90)                 
    ax.set_yticks([0.5, 1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5],
                  labels =  ['Perid', 'Fuco', 'Hex-fuco', 'But-fuco', 'PSC', 'Zea', 'Diato','Diadino','alpha-beta-Car','Allo','PPC','Tot_Chl_c','Tot_Chl_b$^{*}$', 'Tot_Chl_b', '  ', 'Tot_Chl_a'])

    plt.subplot(2,1,2)
    plt.rcParams.update({'font.size': 22})
    C = np.nan*np.ones([len(pig_keys), len(pig_keys)])
    for i in range(len(pig_keys)):
        for j in range(len(pig_keys)):
           if i > j:
             #  breakpoint()
               if (pig_keys[i] != ' ') & (pig_keys[j] != ' '):
                   pig_i = data[pig_keys[i]][(data[pig_keys[i]] > 0) & (data[pig_keys[j]] > 0)]
                   pig_j = data[pig_keys[j]][(data[pig_keys[i]] > 0) & (data[pig_keys[j]] > 0)]
                   C[i, j] = scipy.stats.pearsonr(pig_i, pig_j)[0]**2
                   if pig_keys[i] =='hplc_Tot_Chl_b' and pig_keys[j] == 'hplc_Tot_Chl_b_noPML':
                      C[i, j] = np.nan
                   if pig_keys[j] =='hplc_Tot_Chl_b' and pig_keys[i] == 'hplc_Tot_Chl_b_noPML':
                      C[i, j] = np.nan
 
         
           if i < j:     
               if (pig_keys[i] != ' ') & (pig_keys[j] != ' '):
                   pig_i = data[pig_keys[i]][(data[pig_keys[i]] > 0) & (data[pig_keys[j]] > 0)]
                   pig_j = data[pig_keys[j]][(data[pig_keys[i]] > 0) & (data[pig_keys[j]] > 0)]
                   Tchla = data['hplc_Tot_Chl_a'][(data[pig_keys[i]] > 0) & (data[pig_keys[j]] > 0)]
                   if i > 0:
                       C[i, j] = scipy.stats.pearsonr(pig_i/Tchla, pig_j/Tchla)[0]**2
                   elif i == 0:
                       C[i, j] = scipy.stats.pearsonr(pig_i, pig_j/Tchla)[0]**2
                   if pig_keys[i] =='hplc_Tot_Chl_b' and pig_keys[j] == 'hplc_Tot_Chl_b_noPML':
                        C[i, j] = np.nan
                   if pig_keys[j] =='hplc_Tot_Chl_b' and pig_keys[i] == 'hplc_Tot_Chl_b_noPML':
                        C[i, j] = np.nan

    plt.pcolor(np.flipud(C.T), cmap='Oranges', vmin=0, vmax=1)
    cbar = plt.colorbar()
    cbar.set_label('Coefficient of determination, $r^2$', rotation=90, labelpad=30)

    plt.xlim(0,16)
    plt.ylim(0,16)

    for i in range(C.shape[0]):
         for j in range(C.shape[1]):
             if np.isnan(C[i, j]) == 0:
                 plt.text(i + 0.5, 15 -j + 0.5, '%.2f' % C[i, j],
                          horizontalalignment='center',
                          verticalalignment='center',
                          )
    ax = plt.gca()
    ax.set_xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5], 
    labels =  ['Tot_Chl_a', ' ', 'Tot_Chl_b', 'Tot_Chl_b$^{*}$', 'Tot_Chl_c' , 
                'PPC', 'Allo','alpha-beta-Car', 'Diadino', 'Diato', 'Zea',
                'PSC', 'But-fuco', 'Hex-fuco', 'Fuco', 'Perid'])


    plt.xticks(rotation=90)                 
    ax.set_yticks([0.5, 1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5],
                  labels =  ['Perid', 'Fuco', 'Hex-fuco', 'But-fuco', 'PSC', 'Zea', 'Diato','Diadino','alpha-beta-Car','Allo','PPC','Tot_Chl_c','Tot_Chl_b$^{*}$', 'Tot_Chl_b', ' ', 'Tot_Chl_a'])

    ax = plt.gca()
    ax.text(-0.20, 0.95,  'B', ha='left', va='top', transform=ax.transAxes,fontsize=32)  
        

    plt.tight_layout(pad=1.2)
    filename  =  fig_dir + '/'  + '_pigscovmat_Kramerformat.png'
    plt.savefig(filename,dpi=600)


    return


def plot_pig_cov_v2(data):

    pig_keys = ['hplc_Tot_Chl_a',' ','hplc_Tot_Chl_b', 'hplc_Tot_Chl_c' , 
                'hplc_PPC', 'hplc_Allo','hplc_Alpha-beta-Car', 'hplc_Diadino', 'hplc_Diato', 'hplc_Zea',
                'hplc_PSC', 'hplc_But-fuco', 'hplc_Hex-fuco', 'hplc_Fuco', 'hplc_Perid']


    pig_labels =  ['Tot_Chl_a','Tot_Chl_b', 'Tot_Chl_c' , 
               'PPC', 'Allo','alpha-beta-Car', 'hplc_Diato', 'Diadino', 'Zea',
               'PSC', 'But-fuco', 'Hex-fuco', 'Fuco', 'Perid']

    plt.figure(figsize=(18,22))
    plt.subplot(2,1,1)
    plt.rcParams.update({'font.size': 22})
    C = np.nan*np.ones([len(pig_keys), len(pig_keys)])
    for i in range(len(pig_keys)):
        for j in range(len(pig_keys)):
            if i > j:
              #  breakpoint()
                if (pig_keys[i] != ' ') & (pig_keys[j] != ' '):
                    pig_i = data[pig_keys[i]][(data[pig_keys[i]] > 0) & (data[pig_keys[j]] > 0)]
                    pig_j = data[pig_keys[j]][(data[pig_keys[i]] > 0) & (data[pig_keys[j]] > 0)]
                    P = scipy.stats.linregress(pig_i, pig_j, alternative='two-sided').pvalue
                    if P < 0.01:
                        C[i, j] = scipy.stats.pearsonr(pig_i, pig_j)[0]
          
            if i < j:     
                if (pig_keys[i] != ' ') & (pig_keys[j] != ' '):
                    pig_i = data[pig_keys[i]][(data[pig_keys[i]] > 0) & (data[pig_keys[j]] > 0)]
                    pig_j = data[pig_keys[j]][(data[pig_keys[i]] > 0) & (data[pig_keys[j]] > 0)]
                    Tchla = data['hplc_Tot_Chl_a'][(data[pig_keys[i]] > 0) & (data[pig_keys[j]] > 0)]
                    if i > 0:
                       # breakpoint()
                        P = scipy.stats.linregress(pig_i/Tchla, pig_j/Tchla, alternative='two-sided').pvalue
                        if P < 0.01:
                            C[i, j] = scipy.stats.pearsonr(pig_i/Tchla, pig_j/Tchla)[0]
                    elif i == 0:
                        P = scipy.stats.linregress(pig_i, pig_j/Tchla, alternative='two-sided').pvalue
                        if P < 0.01:
                            C[i, j] = scipy.stats.pearsonr(pig_i, pig_j/Tchla)[0]
   
    plt.rcParams.update({'font.size': 22})
    plt.pcolor(np.flipud(C.T), cmap='bwr', vmin=-1, vmax=1)
    cbar = plt.colorbar()
    cbar.set_label('Correlation coefficient, $r$', rotation=90 ,labelpad=30)

    plt.xlim(0,15)
    plt.ylim(0,15)

    ax = plt.gca()
    ax.text(-0.20, 0.95,  'A', ha='left', va='top', transform=ax.transAxes,fontsize=32)  
        

    for i in range(C.shape[0]):
        for j in range(C.shape[1]):
            if np.isnan(C[i, j]) == 0:
                plt.text(i + 0.5, 14 -j + 0.5, '%.2f' % C[i, j],
                         horizontalalignment='center',
                         verticalalignment='center',
                         )
    ax = plt.gca()
    ax.set_xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5], 
    labels =  ['Tot_Chl_a', ' ', 'Tot_Chl_b', 'Tot_Chl_c' , 
               'PPC', 'Allo','alpha-beta-Car', 'Diadino', 'Diato', 'Zea',
               'PSC', 'But-fuco', 'Hex-fuco', 'Fuco', 'Perid'])


    plt.xticks(rotation=90)                 
    ax.set_yticks([0.5, 1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5],
                  labels =  ['Perid', 'Fuco', 'Hex-fuco', 'But-fuco', 'PSC', 'Zea', 'Diato','Diadino','alpha-beta-Car','Allo','PPC','Tot_Chl_c','Tot_Chl_b', '  ', 'Tot_Chl_a'])

    plt.subplot(2,1,2)
    plt.rcParams.update({'font.size': 22})
    C = np.nan*np.ones([len(pig_keys), len(pig_keys)])
    for i in range(len(pig_keys)):
        for j in range(len(pig_keys)):
           if i > j:
             #  breakpoint()
               if (pig_keys[i] != ' ') & (pig_keys[j] != ' '):
                   pig_i = data[pig_keys[i]][(data[pig_keys[i]] > 0) & (data[pig_keys[j]] > 0)]
                   pig_j = data[pig_keys[j]][(data[pig_keys[i]] > 0) & (data[pig_keys[j]] > 0)]
                   P = scipy.stats.linregress(pig_i, pig_j, alternative='two-sided').pvalue
                   if P < 0.01:
                       C[i, j] = scipy.stats.pearsonr(pig_i, pig_j)[0]**2
         
           if i < j:     
               if (pig_keys[i] != ' ') & (pig_keys[j] != ' '):
                   pig_i = data[pig_keys[i]][(data[pig_keys[i]] > 0) & (data[pig_keys[j]] > 0)]
                   pig_j = data[pig_keys[j]][(data[pig_keys[i]] > 0) & (data[pig_keys[j]] > 0)]
                   Tchla = data['hplc_Tot_Chl_a'][(data[pig_keys[i]] > 0) & (data[pig_keys[j]] > 0)]
                   if i > 0:
                       P = scipy.stats.linregress(pig_i/Tchla, pig_j/Tchla, alternative='two-sided').pvalue
                       if P < 0.01:
                           C[i, j] = scipy.stats.pearsonr(pig_i/Tchla, pig_j/Tchla)[0]**2
                   elif i == 0:
                       P = scipy.stats.linregress(pig_i, pig_j/Tchla, alternative='two-sided').pvalue
                       if P < 0.01:
                           C[i, j] = scipy.stats.pearsonr(pig_i, pig_j/Tchla)[0]**2

    plt.pcolor(np.flipud(C.T), cmap='Oranges', vmin=0, vmax=1)
    cbar = plt.colorbar()
    cbar.set_label('Coefficient of determination, $r^2$', rotation=90, labelpad=30)

    plt.xlim(0,15)
    plt.ylim(0,15)

    for i in range(C.shape[0]):
         for j in range(C.shape[1]):
             if np.isnan(C[i, j]) == 0:
                 plt.text(i + 0.5, 14 -j + 0.5, '%.2f' % C[i, j],
                          horizontalalignment='center',
                          verticalalignment='center',
                          )
    ax = plt.gca()
    ax.set_xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5], 
    labels =  ['Tot_Chl_a', ' ', 'Tot_Chl_b', 'Tot_Chl_c' , 
                'PPC', 'Allo','alpha-beta-Car', 'Diadino', 'Diato', 'Zea',
                'PSC', 'But-fuco', 'Hex-fuco', 'Fuco', 'Perid'])


    plt.xticks(rotation=90)                 
    ax.set_yticks([0.5, 1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5],
                  labels =  ['Perid', 'Fuco', 'Hex-fuco', 'But-fuco', 'PSC', 'Zea', 'Diato','Diadino','alpha-beta-Car','Allo','PPC','Tot_Chl_c','Tot_Chl_b', ' ', 'Tot_Chl_a'])

    ax = plt.gca()
    ax.text(-0.20, 0.95,  'B', ha='left', va='top', transform=ax.transAxes,fontsize=32)  
        

    plt.tight_layout(pad=1.2)

    filename  =  fig_dir + '/'  + '_pigscovmat_Kramerformat.png'
    plt.savefig(filename,dpi=600)


    return



# Chl:HPLC plot, Chl transects/histograms
def plot_total_ABfit():
        '''routine to do power-law fit & plot for combined dataset (Fig 3 in paper)'''
                
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
        0.014*data_nc_23.acs_chl_debiased.acs_chl , 
        0.014*data_nc_24.acs_chl_debiased.acs_chl , 
        0.014*data_nc_24.acs2_chl_debiased.acs_chl,
        0.014*data_nc_25.acs_chl_debiased.acs_chl ,
        # data_nc_25.acs2_chl_debiased.acs_chl; 
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
                    
        # clip
        tot_acs = tot_acs[~np.isnan(tot_chl)]
        tot_chl = tot_chl[~np.isnan(tot_chl)]
        
        tot_chl = tot_chl[~np.isnan(tot_acs)]
        tot_acs = tot_acs[~np.isnan(tot_acs)]
        
        # linear_mod = scipy.stats.linregress(np.log10(0.014*tot_acs),np.log10(tot_chl))
                                                  
        linear_mod = scipy.stats.linregress(np.log10(tot_acs),np.log10(tot_chl)) #- no A
                                                  
        log10A = np.round(100*linear_mod.intercept)/100 # note - A_ 
        A = 10**log10A
        B = np.round(1000*linear_mod.slope)/1000
        r_sq = np.round(1000*linear_mod.rvalue**2)/100
        
        r_sq = np.round(1000*linear_mod.rvalue**2)/1000
        stderr = np.round(1000*linear_mod.stderr)/1000
        interr = np.round(1000*linear_mod.intercept_stderr)/1000
        
        
        print('A = ' + str(10**log10A) + ' +/- ' + str(2*10**interr)) #
        print('B = ' + str(B) + ' +/- ' + str(2*stderr))
        print(r_sq)
        
        
        rres = tot_acs/(0.014*tot_chl) - 1
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
        Y = X*0.014 # used in 1-1 plot
        Y2 = (X/A)**(1/B)
        colors = cm.Paired(np.linspace(0,1,10))
           
           
        plt.figure(figsize=(13, 6))
        plt.subplot(1,2,1)
        plt.rcParams.update({'font.size': 18})
        plt.scatter(0.014*data_nc_19.acx_chl_debiased.acx_chl, data_nc_19.acx_chl_debiased.HPLC_Tot_chla, color=colors[0],alpha=0.5, s=8, label = 'AMT 19')
        plt.scatter(0.014*data_nc_22.acs_chl_debiased.acs_chl, data_nc_22.acs_chl_debiased.HPLC_Tot_chla, color=colors[1], alpha=0.5, s=8, label = 'AMT 22')
        plt.scatter(0.014*data_nc_23.acs_chl_debiased.acs_chl, data_nc_23.acs_chl_debiased.HPLC_Tot_chla, color=colors[2], alpha=0.5, s=8, label = 'AMT 23')
        plt.scatter(0.014*data_nc_24.acs_chl_debiased.acs_chl, data_nc_24.acs_chl_debiased.HPLC_Tot_chla, color=colors[3], alpha=0.5, s=8)
        plt.scatter(0.014*data_nc_24.acs2_chl_debiased.acs_chl, data_nc_24.acs2_chl_debiased.HPLC_Tot_chla, color=colors[3], alpha=0.5, s=8, label = 'AMT 24')
        plt.scatter(0.014*data_nc_25.acs_chl_debiased.acs_chl, data_nc_25.acs_chl_debiased.HPLC_Tot_chla, color=colors[4], alpha=0.5, s=8, label = 'AMT 25')
        plt.scatter(0.014*data_nc_26.acs_chl_debiased.acs_chl, data_nc_26.acs_chl_debiased.HPLC_Tot_chla, color=colors[5], alpha=0.5, s=8, label = 'AMT 26')
        plt.scatter(0.014*data_nc_27.acs_chl_debiased.acs_chl, data_nc_27.acs_chl_debiased.HPLC_Tot_chla, color=colors[6], alpha=0.5, s=8, label = 'AMT 27')
        plt.scatter(0.014*data_nc_28.acx_chl_debiased.acx_chl, data_nc_28.acx_chl_debiased.HPLC_Tot_chla,  color=colors[7], alpha=0.5, s=8, label = 'AMT 28')
        plt.scatter(0.014*data_nc_29.acs_chl_debiased.acs_chl, data_nc_29.acs_chl_debiased.HPLC_Tot_chla, color=colors[9], alpha=0.5, s=8, label = 'AMT 29')
        plt.yscale('log')
        plt.xscale('log')    
        plt.grid('on', ls='--')
        plt.plot(Y,X,'k',linestyle='dotted', label='0.014:1')
        plt.plot(Y2,X,'magenta',linestyle='dashed', label='Best fit')
        plt.xlabel('$a_{ph}(676)$ [m$^{-1}]$')
        #plt.ylabel('$C_{a}$(HPLC)  [mg m$^{-3}]$')
        plt.ylabel('$C_{a}$(HPLC) [mg m$^{-3}]$')
        plt.legend(fontsize=14)
        plt.axis('equal')
    
        plt.xlim(0.001,0.01)
        plt.ylim(0.01,3)

        ax = plt.gca()
        plt.text(.05, .95,  'A', ha='left', va='top', transform=ax.transAxes,fontsize=22)  
       

        plt.subplot(1,2,2)
        resid_19 = 100*(data_nc_19.acx_chl_debiased.HPLC_Tot_chla - A*(0.014*data_nc_19.acx_chl_debiased.acx_chl)**B)/(A*(0.014*data_nc_19.acx_chl_debiased.acx_chl)**B)
        resid_22 = 100*(data_nc_22.acs_chl_debiased.HPLC_Tot_chla - A*(0.014*data_nc_22.acs_chl_debiased.acs_chl)**B)/(A*(0.014*data_nc_22.acs_chl_debiased.acs_chl)**B)
        resid_23 = 100*(data_nc_23.acs_chl_debiased.HPLC_Tot_chla - A*(0.014*data_nc_23.acs_chl_debiased.acs_chl)**B)/(A*(0.014*data_nc_23.acs_chl_debiased.acs_chl)**B)
        resid_24 = 100*(data_nc_24.acs_chl_debiased.HPLC_Tot_chla - A*(0.014*data_nc_24.acs_chl_debiased.acs_chl)**B)/(A*(0.014*data_nc_24.acs_chl_debiased.acs_chl)**B)
        resid_24_2 = 100*(data_nc_24.acs2_chl_debiased.HPLC_Tot_chla - A*(0.014*data_nc_24.acs2_chl_debiased.acs_chl)**B)/(A*(0.014*data_nc_24.acs2_chl_debiased.acs_chl)**B)
        resid_25 = 100*(data_nc_25.acs_chl_debiased.HPLC_Tot_chla - A*(0.014*data_nc_25.acs_chl_debiased.acs_chl)**B)/(A*(0.014*data_nc_25.acs_chl_debiased.acs_chl)**B)
        resid_26 = 100*(data_nc_26.acs_chl_debiased.HPLC_Tot_chla - A*(0.014*data_nc_26.acs_chl_debiased.acs_chl)**B)/(A*(0.014*data_nc_26.acs_chl_debiased.acs_chl)**B)
        resid_27 = 100*(data_nc_27.acs_chl_debiased.HPLC_Tot_chla - A*(0.014*data_nc_27.acs_chl_debiased.acs_chl)**B)/(A*(0.014*data_nc_27.acs_chl_debiased.acs_chl)**B)
        resid_28 = 100*(data_nc_28.acx_chl_debiased.HPLC_Tot_chla - A*(0.014*data_nc_28.acx_chl_debiased.acx_chl)**B)/(A*(0.014*data_nc_28.acx_chl_debiased.acx_chl)**B)
        resid_29 = 100*(data_nc_29.acs_chl_debiased.HPLC_Tot_chla - A*(0.014*data_nc_29.acs_chl_debiased.acs_chl)**B)/(A*(0.014*data_nc_29.acs_chl_debiased.acs_chl)**B)
        
        plt.rcParams.update({'font.size': 16})
        plt.scatter(0.014*data_nc_19.acx_chl_debiased.acx_chl, resid_19, color=colors[0],alpha=0.5, s=8, label = 'AMT 19')
        plt.scatter(0.014*data_nc_22.acs_chl_debiased.acs_chl, resid_22, color=colors[1],alpha=0.5, s=8, label = 'AMT 22')
        plt.scatter(0.014*data_nc_23.acs_chl_debiased.acs_chl, resid_23, color=colors[2],alpha=0.5, s=8, label = 'AMT 23')
        plt.scatter(0.014*data_nc_24.acs_chl_debiased.acs_chl, resid_24, color=colors[3],alpha=0.5, s=8, label = 'AMT 24')
        plt.scatter(0.014*data_nc_24.acs2_chl_debiased.acs_chl, resid_24_2, color=colors[3],alpha=0.5, s=8)
        plt.scatter(0.014*data_nc_25.acs_chl_debiased.acs_chl, resid_25, color=colors[4],alpha=0.5, s=8, label = 'AMT 25')
        plt.scatter(0.014*data_nc_26.acs_chl_debiased.acs_chl, resid_26, color=colors[5],alpha=0.5, s=8, label = 'AMT 26')
        plt.scatter(0.014*data_nc_27.acs_chl_debiased.acs_chl, resid_27, color=colors[6],alpha=0.5, s=8, label = 'AMT 27')
        plt.scatter(0.014*data_nc_28.acx_chl_debiased.acx_chl, resid_28, color=colors[7],alpha=0.5, s=8, label = 'AMT 28')
        plt.scatter(0.014*data_nc_29.acs_chl_debiased.acs_chl, resid_29, color=colors[9],alpha=0.5, s=8, label = 'AMT 29')
    

        plt.xscale('log')      
        plt.xlabel('$a_{ph}(676)$ [m$^{-1}]$')
        plt.ylabel('$C_{a}$(HPLC) percentage residual [$\%$]')
      #  plt.legend(fontsize=13)
        plt.xscale('log') 
        plt.tight_layout()
        plt.grid('on', ls='--')
        #plt.xlim(0.001,0.01)
    
        # plt.xlim(0.001,0.01)
        plt.ylim(-100,100)
        
        ax = plt.gca()
        plt.text(.05, .95,  'B', ha='left', va='top', transform=ax.transAxes,fontsize=22)  
 
        
        filename  =  fig_dir + '/'  + '_ABscatterplot.png'
        plt.savefig(filename,dpi=600)
       
        return
        

def plot_chl_hist():
      'Chl-histogram plot'
    
      plt.figure(figsize=(13,5))
      plt.rcParams.update({'font.size': 16})
      
      plt.subplot(1,2,1)
      chl = _9cruise_chlsum()
      bins = np.arange(0.01,10,0.25)
      hist, bins = np.histogram(chl, bins=bins)
      logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
      plt.hist(chl, bins=logbins,edgecolor='black',label='Total $N_{ACS/AC9}$ = ' + str(np.sum(~np.isnan(chl))))
      plt.xscale('log')
      plt.ylabel('$N_{ACS/AC9}$')
      plt.xlabel('Tot_Chl_a concentration  [mg m$^{-3}]$')
      plt.xlim(0.01,10)
      
      plt.legend()
      ax=plt.gca()
      plt.text(0.05, .95,  'A', ha='left', va='top', transform=ax.transAxes,fontsize=22)  
        
      plt.subplot(1,2,2)
      chl = df_hplc['hplc_Tot_Chl_a']
      bins = np.arange(0.01,10,0.25)
      hist, bins = np.histogram(chl, bins=bins)
      logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
      plt.hist(chl, bins=logbins,edgecolor='black',color='orange',label='Total $N_{HPLC}$ = ' + str(len(df_hplc)) )
      plt.xscale('log')
      plt.ylabel('$N_{HPLC}$')
      plt.xlabel('Tot_Chl_a concentration  [mg m$^{-3}]$')
      plt.xlim(0.01,10)
  
      plt.legend()
      ax=plt.gca()
      plt.text(0.05, 0.95,  'B', ha='left', va='top', transform=ax.transAxes,fontsize=22)  
        
      filename  =  fig_dir + '/'  + '_TChla_hist.png'
      plt.savefig(filename,dpi=600,bbox_inches='tight')

      return


def _color_by_prov_chl(data_nc, mask, plot_index):
    
    plt.subplot(5,2,plot_index)
    letters = ['A: AMT 19', 'B: AMT 22' , 'C: AMT 23', 'D: AMT 24', 'E: AMT 25', 'F: AMT 26', 'G: AMT 27' ,'H: AMT 28', 'I: AMT 29']
        
    if plot_index == 1 or plot_index == 8:    
        acs_key ='acx_chl_debiased'
    else:
        acs_key ='acs_chl_debiased'
    if plot_index == 4 or plot_index == 5:    
        acs2_key = 'acs2_chl_debiased'
    
    
    prov_extra= ['NECS', 'FKLD', 'ANTA']
    for i in range(len(prov_extra)):
        plt.scatter(np.array(data_nc['uway_lat'])[mask==prov_extra[i]], np.array(data_nc[acs_key])[mask==prov_extra[i]], s=2, color='gray')
        if plot_index == 4 or plot_index == 5:    
            plt.scatter(np.array(data_nc['uway_lat'])[mask==prov_extra[i]], np.array(data_nc[acs2_key])[mask==prov_extra[i]], s=2, color='gray')
  
    prov= ['NADR', 'NASE', 'NASW', 'NATR', 'WTRA', 'SATL', 'SSTC', 'SANT'] 
    colors = cm.Paired(np.linspace(0,1,10))
    
    for i in range(len(prov)):
        plt.scatter(np.array(data_nc['uway_lat'])[mask==prov[i]], np.array(data_nc[acs_key])[mask==prov[i]],color=colors[i], s=2)
        if plot_index == 4 or plot_index == 5:  
            plt.scatter(np.array(data_nc['uway_lat'])[mask==prov[i]], np.array(data_nc[acs2_key])[mask==prov[i]],color=colors[i], s=2)
    
        
    #if plot_index == 1 or plot_index == 8:
    #      plt.scatter(data_nc.acx_chl_debiased.match_up_dates, data_nc.acx_chl_debiased.HPLC_Tot_chla, s=10, color='black')
    #  elif plot_index == 4 or plot_index == 5: 
    #      plt.scatter(data_nc.acs2_chl_debiased.match_up_dates,data_nc.acs2_chl_debiased.HPLC_Tot_chla, s=10,color='black')    
    # else:
    #      plt.scatter(data_nc.acs_chl_debiased.match_up_dates,data_nc.acs_chl_debiased.HPLC_Tot_chla, s=10,color='black')

        
    plt.scatter(np.array(data_nc['hplc_lat'])[data_nc['hplc_depth']<10], np.array(data_nc['hplc_Tot_Chl_a'])[data_nc['hplc_depth']<10],s=10, color='black')

      # plt.scatter(data_nc['hplc_lat'], data_nc['hplc_Tot_Chl_a'], s=10, color='black')

    plt.yscale('log')
    plt.ylim(0.01,7)
    ax=plt.gca()
    plt.xlim(-56,52)
    plt.gca().invert_xaxis()
    # ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%m-%d') )
    #ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator(bymonthday=[1,5,10,15,20,25]))
    #plt.xticks(rotation=45)

    plt.text(.02, .94,  letters[plot_index-1], ha='left', va='top', transform=ax.transAxes,fontsize=20) 

   # if plot_index == 8:
  

    return


def plot_chl_time_series():
    ' plot for chl transects (note these are plotted against latitude NOT time)'
    
    fig = plt.figure(figsize=(16,16))
    plt.rcParams.update({'font.size': 16})
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
    fig.supxlabel('Latitude [$^{\circ}$]')
    
    colors = cm.Paired(np.linspace(0,1,10)) 
    prov = ['IOP: NADR', 'IOP: NASE', 'IOP: NASW', 'IOP: NATR', 'IOP: WTRA', 'IOP: SATL', 'IOP: SSTC', 'IOP: SANT']     
    patches = []
    for k in range(8):
        color_string = cl.rgb2hex(colors[k])
        patches.append(mpatches.Patch(color=color_string, label=prov[k]))
    
    patches.append(mpatches.Patch(color='gray', label='IOP: OTHER'))
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
    
    filename  =  fig_dir + '/'  + '_Chl_transects.png'
    plt.savefig(filename,dpi=900)

    return


# Plots for Sect 3.5 (variaton in ap spectra)
def hplc_ap_match_up_V3(data_nc, two_acs_systems=False):
    
    ' extracts ap spectra that match hplc timestamps - median average over 30 min bin.'
    ' this is designed to replicate what was done in the tot-chl-a match up procedure'    
    
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



def plot_ap_limitingcases():
    ' ap for pigment ratios plot (sect 3.5)'

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
    
    global_med, global_med_int = int_globalmed()
    wv = data_nc_26['acs_wv'].values
    
    pig_keys = ['hplc_Tot_Chl_b', 'hplc_Tot_Chl_c' , 
               'hplc_PPC', 'hplc_PSC']

    style_vec  = ['solid','dotted','dashed','dashdot','solid']
    
    
    plt.rcParams.update({'font.size': 18})    
    plt.figure(figsize=(15,5))
    plt.subplot(1,3,1)

    plt.plot(wv,global_med, color='black')    
    for i in range(len(pig_keys)):
        pig_name =  pig_keys[i]
        p_threshold = np.percentile(df_hplc_combined[pig_name]/df_hplc_combined['hplc_Tot_Chl_a'], 90)
        threshold_index = np.where(df_hplc_combined[pig_name]/df_hplc_combined['hplc_Tot_Chl_a'] > p_threshold)
        df_ap_threshold = df_ap_combined.iloc[threshold_index]
        ap_threshold = np.array(df_ap_threshold)
        
        plt.plot(wv,np.nanmedian(ap_threshold,axis=0),linestyle=style_vec[i],linewidth=2)
        df_ap_threshold = []
        ap_threshold =[]
        
    #plt.legend()
    ax = plt.gca()
    plt.text(.90, .95,  'A', ha='left', va='top', transform=ax.transAxes,fontsize=22) 
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('$a_{p}(\lambda)$ [m$^{-1}$]')
    plt.xlim(400,720)
 
    
    p_labels=[]
    plt.subplot(1,3,2)
    plt.plot(wv,global_med_int, color='black')    
    for i in range(len(pig_keys)):
     
        pig_name =  pig_keys[i]
        p_threshold = np.percentile(df_hplc_combined[pig_name]/df_hplc_combined['hplc_Tot_Chl_a'], 90)
        threshold_index = np.where(df_hplc_combined[pig_name]/df_hplc_combined['hplc_Tot_Chl_a'] > p_threshold)
        df_ap_threshold = df_ap_int_combined.iloc[threshold_index]
        ap_threshold = np.array(df_ap_threshold)
        p_labels.append((1/100)*np.round(p_threshold*100))
        plt.plot(wv,np.nanmedian(ap_threshold,axis=0),linestyle=style_vec[i],linewidth=2)
        df_ap_threshold = []
        ap_threshold =[]
    ax = plt.gca()
    plt.text(.90, .95,'B', ha='left', va='top', transform=ax.transAxes,fontsize=22)  
   #  plt.legend(fontsize=14,loc=7)
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('$<a_{p}(\lambda)>$')
    plt.xlim(400,720)
 
    # plt.grid()
    #  breakpoint()
    labels =  ['Median(Tot_Chl_b/Tot_Chl_a > $P_{90}$),  ($P_{90}$ = ' + str(p_labels[0]) +')', 'Median(Tot_Chl_c/Tot_Chl_a > $P_{90}$), ($P_{90}$ = ' + str(p_labels[1]) +')',
               'Median(PPC/Tot_Chl_a > $P_{90}$), ($P_{90}$ = ' + str(p_labels[2]) +')', 'Median(PSC/Tot_Chl_a > $P_{90}$), ($P_{90}$ = ' + str(p_labels[3]) +')']
               
               
    plt.subplot(1,3,3)
    plt.plot(wv,np.zeros(len(wv)),label='AMT median', color='black')    
    for i in range(len(pig_keys)):
    
        pig_name =  pig_keys[i]
        p_threshold = np.percentile(df_hplc_combined[pig_name]/df_hplc_combined['hplc_Tot_Chl_a'], 90)
        threshold_index = np.where(df_hplc_combined[pig_name]/df_hplc_combined['hplc_Tot_Chl_a'] > p_threshold)
        df_ap_threshold = df_ap_int_combined.iloc[threshold_index]
        ap_threshold = np.array(df_ap_threshold)
        ap_resid = np.nan*np.ones([len(ap_threshold), len(ap_threshold.T)])
        for j in range(len(ap_threshold)):
              ap_resid[j,:] = (ap_threshold[j,:] - global_med_int)
        plt.plot(wv,np.nanmedian(ap_resid,axis=0),label=labels[i],linestyle=style_vec[i],linewidth=2)
        
        df_ap_threshold = []
        df_threshold =[]
        ap_resid[j,:] 
    
    #plt.legend()
    ax = plt.gca()
    plt.text(.90, .95, 'C', ha='left', va='top', transform=ax.transAxes,fontsize=22)
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('$<a_{p}(\lambda)>$ residual')
    plt.xlim(400,720)
    plt.tight_layout()
    plt.grid()

    plt.figlegend(loc='upper center', bbox_to_anchor=(0.5, -0.02),
    fancybox=True, ncol=2, fontsize=13)

    
    filename  =  fig_dir + '/'  + '_limitingcases.png'
    plt.savefig(filename,dpi=600,bbox_inches='tight')

    return


if __name__ == '__main__':

    
    fig_dir = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/paperplots_essd/'

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
     

    # load Longhurst Province masks
    fn_mask_19 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/Longhurst_masks/LH_mask_AMT19.csv'
    mask_19 = pd.read_csv(fn_mask_19)
    
    fn_mask_22 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/Longhurst_masks/LH_mask_AMT22.csv'
    mask_22 = pd.read_csv(fn_mask_22)
    
    fn_mask_23 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/Longhurst_masks/LH_mask_AMT23.csv'
    mask_23 = pd.read_csv(fn_mask_23)
    
    fn_mask_24 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/Longhurst_masks/LH_mask_AMT24.csv'
    mask_24 = pd.read_csv(fn_mask_24)
    
    fn_mask_25 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/Longhurst_masks/LH_mask_AMT25.csv'
    mask_25 = pd.read_csv(fn_mask_25)
    
    fn_mask_26 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/Longhurst_masks/LH_mask_AMT26.csv'
    mask_26 = pd.read_csv(fn_mask_26)
        
    fn_mask_27 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/Longhurst_masks/LH_mask_AMT27.csv'
    mask_27 = pd.read_csv(fn_mask_27)
        
    fn_mask_28 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/Longhurst_masks/LH_mask_AMT28.csv'
    mask_28 = pd.read_csv(fn_mask_28)
    
    fn_mask_29 = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/Longhurst_masks/LH_mask_AMT29.csv'
    mask_29 = pd.read_csv(fn_mask_29)
    
    fn_mask_hplc = '/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/Longhurst_masks/LH_mask_HPLC_AMTAMT.csv'
    LH_HPLC = pd.read_csv(fn_mask_hplc)

    _9cruise_LH_mask(field = 'LH_Province') # stacks LH mask for 9 cruise dataset
    df_hplc = _pig_9cruises() # stacks pigments into dataframe


    # Generate figures as used in MS submission
    # Methods figs - Sect. 2
    plot_coverage() # Fig 1 A
    plot_AMT_timeline() # Fig 1 B
  #   # see step5/AMT 28 for Fig 2 plots                    
    plot_total_ABfit() # Fig 3
 #   
    # IOP results figs - Sect 3.1
    plot_median_ap_province() # Fig 4
    plot_676_443() # Fig 5
    plot_uncertainties() # Fig 6
#
    # Tchla results - Sect 3.2
    plot_chl_time_series() # Fig 7
    plot_chl_hist()  # Fig 8
 
    # Pigment results - Sect 3.3
    plot_pigs_byprov() # Fig 9 
    plot_pigs_byprov_ratio()  # Fig 10

    # Pigment correlationss - Sect 3.4
    plot_pig_cov(df_hplc)
    plot_pig_cov_v2(df_hplc)
    
    # ap limting case - Sect 3.5
    plot_ap_limitingcases()
    
    # code example to generate LH mask
   # mask19 = longhurst_mask(data_nc_19['uway_lat'], data_nc_19['uway_lon'], 'AMT19')
    #mask19.to_csv('/data/abitibi1/scratch/scratch_disk/tjor/AMT_underway/Tests_and_plots/Longhurst_masks/LH_mask_AMT19_v2.csv' ,index=False)

    