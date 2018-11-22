#Animation of ensemble simulations for ElEvoHI

# Author: C. Moestl, IWF Graz, Austria
# twitter @chrisoutofspace, https://github.com/cmoestl
# November 2018
# This work is published under the MIT LICENSE (see bottom)

import numpy as np
import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import cm
from scipy import stats
import scipy.io
import sunpy.time
import time
import pickle
import seaborn as sns
import math

####################################################### functions

#for reading catalogues  
def getcat(filename):
  print( 'reading CAT '+filename)
  cat=scipy.io.readsav(filename)#, verbose='false')  
  print( 'done reading CAT')
  return cat  
  
def decode_array(bytearrin):
 #for decoding the strings from the IDL .sav file to a list of python strings, not bytes 
 #make list of python lists with arbitrary length
 bytearrout= ['' for x in range(len(bytearrin))]
 for i in range(0,len(bytearrin)-1):
  bytearrout[i]=bytearrin[i].decode()
 #has to be np array so to be used with numpy "where"
 bytearrout=np.array(bytearrout)
 return bytearrout  
  
def time_to_num_cat(time_in):  
  #for time conversion from catalogue .sav to numerical time
  #this for 1-minute data or lower time resolution
  #for all catalogues
  #time_in is the time in format: 2007-11-17T07:20:00 or 2007-11-17T07:20Z
  #for times help see: 
  #http://docs.sunpy.org/en/latest/guide/time.html
  #http://matplotlib.org/examples/pylab_examples/date_demo2.html
  
  j=0
  #time_str=np.empty(np.size(time_in),dtype='S19')
  time_str= ['' for x in range(len(time_in))]
  #=np.chararray(np.size(time_in),itemsize=19)
  time_num=np.zeros(np.size(time_in))
  
  
  for i in time_in:
   #convert from bytes (output of scipy.readsav) to string
   time_str[j]=time_in[j][0:16].decode()+':00'
   year=int(time_str[j][0:4])
   time_str[j]
   #convert time to sunpy friendly time and to matplotlibdatetime
   #only for valid times so 9999 in year is not converted
   #pdb.set_trace()
   if year < 2100:
    	  time_num[j]=mdates.date2num(sunpy.time.parse_time(time_str[j]))
   j=j+1  
   #the date format in matplotlib is e.g. 735202.67569444
   #this is time in days since 0001-01-01 UTC, plus 1.
   
   #return time_num which is already an array and convert the list of strings to an array
  return time_num, np.array(time_str)



########################################################################################
################################# main program ##########################################
########################################################################################


###################################### CONTROLS

#directory of current event

current_event='Nov2010'

#IDL sav file with ensemble simulation results
ensemble_results='formovie_all_flag.sav'

#set 1 on the first run to produce .p save files for interpolated variables needed for the movie
read_data=0

#how much time is between frames in days
dayjump=np.double(1/24.0) 

#how long the movie takes in days
duration_days=2

movie_start_date_time='2010-Nov-3 18:00:00'
cme_start_date_time='2010-Nov-4 08:00:00'

#how long an in situ arrival stays visible in fade mode
fadedays=20

#font size on bottom labels
labelfontsize=8

bscale=4

#whether HEEQ or HEE positions are used for plot
HEEQ=0
HEE=1

#save file with elongation tracks
#tracksav='track_B_img_ccsds.sav'	

##########################################

plt.close('all')
current_event_dir='events/'+current_event+'/'
if os.path.isdir('movies') == False: os.mkdir('movies')
if os.path.isdir(current_event_dir+'/frames') == False: os.mkdir(current_event_dir+'/frames')

print()
print( 'Start ElEvoHI animation program.')
print()
print('Current event ', current_event)
print()



##########get ICMECAT
filename_icmecat='cats/HELCATS_ICMECAT_v10_SCEQ.sav'
i=getcat(filename_icmecat)
print()

#get parameters
bmean=i.icmecat['MO_BMEAN']*bscale #bscale makes circles larger in movie
long=i.icmecat['SC_LONG_HEEQ']*np.pi/180 #hee longitude converted to radians
rdist=i.icmecat['sc_heliodistance'] #AU
sc=i.icmecat['sc_insitu'] #string
sc=decode_array(sc)

#get indices of events in different spacecraft
vexind=np.where(sc == 'VEX')
staind=np.where(sc == 'STEREO-A')
stbind=np.where(sc == 'STEREO-B')
winind=np.where(sc == 'Wind')
mesind=np.where(sc == 'MESSENGER')
ulyind=np.where(sc == 'ULYSSES')

#make time conversion for all icme_start_time variables
#save it as string
icme_start_time_str=i.icmecat['icme_start_time']
#save it as matplotlib date number
[icme_start_time_num,icme_start_time_str]=time_to_num_cat(icme_start_time_str)

#for each spacecraft, make a zeros array 
active_icme_vex=np.zeros(np.size(icme_start_time_num))
active_icme_stb=np.zeros(np.size(icme_start_time_num))
active_icme_sta=np.zeros(np.size(icme_start_time_num))
active_icme_win=np.zeros(np.size(icme_start_time_num))
active_icme_mes=np.zeros(np.size(icme_start_time_num))
active_icme_uly=np.zeros(np.size(icme_start_time_num))


#####get spacecraft and planet positions
if HEEQ == 1: pos=getcat('cats/positions_2007_2023_HEEQ_6hours.sav')
if HEE == 1:  pos=getcat('cats/positions_2007_2023_HEE_6hours.sav')
pos_time_num=time_to_num_cat(pos.time)[0]
#positions are available as pos.mercury etc.

print()

#define times
CME_start_time=mdates.date2num(sunpy.time.parse_time(cme_start_date_time))
#define cme frame times
h_time_num=np.arange(CME_start_time,CME_start_time+duration_days,dayjump)
h_time_str=mdates.num2date(h_time_num)



# ########### read and interpolate e-t profile to movie frame times - used for making line from spacecraft to front
# 
# #et_time_num
# #h_time_num
# #et_elon
# #h_et_elon
# 
# 
# # get elongation-time profile from track
# et=getcat(current_event_dir+tracksav)
# et_time=et.track.track_date[0]
# et_time_num=time_to_num_cat(et_time)[0]
# et_elon= et.track['elon'][0]
# 
# #linearly interpolate to hourly values make automatic later
# et_start_time=mdates.date2num(sunpy.time.parse_time(cme_start_date_time))
# et_time_num_interp=np.arange(et_start_time,et_start_time+10,dayjump)
# et_elon_interp= np.interp(et_time_num_interp, et_time_num, et_elon)
# 

############### read file with ensemble results, dump as pickle to use later
if read_data ==1:

  h=getcat(current_event_dir+ensemble_results)
  all_apex_t=h.elevo_kin.all_apex_t[0]
  [all_apex_t_num_non_interp,all_apex_t_num_non_interp_str]=time_to_num_cat(all_apex_t)
  #get all parameters
  all_apex_r_non_interp=h.elevo_kin.all_apex_r[0]
  all_apex_lat_non_interp=h.elevo_kin.all_apex_lat[0]  #degree
  all_apex_lon_non_interp=h.elevo_kin.all_apex_lon[0]  #degree
  #f
  all_apex_f_non_interp=h.elevo_kin.all_apex_f[0]
  #width
  all_apex_w_non_interp=np.deg2rad(h.elevo_kin.all_apex_w[0])
  
  #constants
  all_apex_s_non_interp=decode_array(h.elevo_kin.all_apex_s[0])
  all_apex_run_non_interp=h.elevo_kin.runnumber[0]
  all_apex_flag_non_interp=h.elevo_kin.colorflag[0]

  #go through each run and interpolate data for each run
  #final array size -> time array of CME frames * run numbers
  finarrs=np.size(h_time_num)*np.max(all_apex_run_non_interp)
  
  eventsize=np.size(h_time_num)
  
  #initialise arrays 
  all_apex_t=np.zeros(finarrs)
  all_apex_r=np.zeros(finarrs)
  all_apex_lat=np.zeros(finarrs)
  all_apex_lon=np.zeros(finarrs)
  all_apex_f=np.zeros(finarrs)
  all_apex_w=np.zeros(finarrs)
  all_apex_s=['']*finarrs
  all_apex_run=np.zeros(finarrs)
  all_apex_flag=np.zeros(finarrs)
  
  print()
  print('start interpolation')
  for q in np.arange(0,np.max(all_apex_run_non_interp)):
    
    #print(q)
    #get indices of kinematic data for this run
    thisrunind=np.where(all_apex_run_non_interp == q)  
  
    #if there is data available for this run, interpolate to CME times
    if np.size(thisrunind) >0:
    
      #these variables change with time
      
      #this is time, fill with frame times
      all_apex_t[eventsize*q:eventsize*(q+1)]=h_time_num
      
      #fill with interpolation variables
      all_apex_r[eventsize*q:eventsize*(q+1)] = np.interp(h_time_num, all_apex_t_num_non_interp[thisrunind],all_apex_r_non_interp[thisrunind]) 
      all_apex_lon[eventsize*q:eventsize*(q+1)] = np.interp(h_time_num, all_apex_t_num_non_interp[thisrunind],all_apex_lon_non_interp[thisrunind]) 
      all_apex_lat[eventsize*q:eventsize*(q+1)] = np.interp(h_time_num, all_apex_t_num_non_interp[thisrunind],all_apex_lat_non_interp[thisrunind]) 
      all_apex_f[eventsize*q:eventsize*(q+1)] = np.interp(h_time_num, all_apex_t_num_non_interp[thisrunind],all_apex_f_non_interp[thisrunind]) 
      all_apex_w[eventsize*q:eventsize*(q+1)] = np.interp(h_time_num, all_apex_t_num_non_interp[thisrunind],all_apex_w_non_interp[thisrunind]) 
      
      #fill with run numbers
      all_apex_run[eventsize*q:eventsize*(q+1)] = all_apex_run_non_interp[thisrunind][0:eventsize]

      #fill with flag numbers
      all_apex_flag[eventsize*q:eventsize*(q+1)] = all_apex_flag_non_interp[thisrunind][0:eventsize]
      
      #fill with observatory string
      all_apex_s[eventsize*q:eventsize*(q+1)] = all_apex_s_non_interp[thisrunind][0:eventsize] 
      
    else: #set all to np.nan for this run
       all_apex_t[eventsize*q:eventsize*(q+1)]=np.nan
       all_apex_r[eventsize*q:eventsize*(q+1)] = np.nan 
       all_apex_lon[eventsize*q:eventsize*(q+1)] = np.nan 
       all_apex_lat[eventsize*q:eventsize*(q+1)] = np.nan 
       all_apex_f[eventsize*q:eventsize*(q+1)] = np.nan 
       all_apex_w[eventsize*q:eventsize*(q+1)] = np.nan
       all_apex_run[eventsize*q:eventsize*(q+1)] = np.nan      
       all_apex_s[eventsize*q:eventsize*(q+1)] = ''
       all_apex_flag[eventsize*q:eventsize*(q+1)] = np.nan
      
  print('end interpolation')
  pickle.dump((all_apex_t,all_apex_r, all_apex_lat, all_apex_lon,all_apex_f,all_apex_w,all_apex_s, all_apex_run,all_apex_flag), open( current_event_dir+"all_apex_variables.p", "wb" ) )

if read_data == 0:
  [all_apex_t,all_apex_r, all_apex_lat, all_apex_lon,all_apex_f,all_apex_w,all_apex_s, all_apex_run, all_apex_flag] = pickle.load( open(current_event_dir+'all_apex_variables.p', "rb" ) )




################################### MAKE MOVIE FRAMES

#initiate plot
plt.figure(1, figsize=(8, 6), dpi=100, facecolor='w', edgecolor='w')


sns.set_context('talk')  
sns.set_style('darkgrid')

#set start time of movie 
frame_time_num=mdates.date2num(sunpy.time.parse_time(movie_start_date_time))


###### loop over all movie frames
for k in np.arange(0,duration_days,dayjump):  
  
 #to current frame time, the days need to be added, so +k is done 
 #save frame time as string to write on plot
 
 framestr = '%04i' % np.round(k*1.0/dayjump) 
 frame_time_str=str(mdates.num2date(frame_time_num+k))

 print( 'frame ', framestr,'  ', frame_time_str)
 #difference array of current frame time frame_time_num+k to position time frame_time_num
 cmedt=frame_time_num+k-all_apex_t
 #get indices where difference is less than half the time resolution
 #use this to avoid nan in np.where
 cmedt[np.isnan(cmedt)]=10000
 cmeind=np.where(np.abs(cmedt) < dayjump/2)
 #print( 'cmeind', cmeind)
  

  
 ############################################ make plot 
 
 ax = plt.subplot(111,projection='polar')
 
 #difference array of current frame time frame_time_num+k to position time frame_time_num
 dct=frame_time_num+k-pos_time_num
 #get index of closest to 0, use this for position
 timeind=np.argmin(abs(dct))
 #print('index pos')
 #print(timeind) 


 ############################### plot all active CME ellipses
 if np.size(cmeind) >0:
  for p in range(0,np.size(cmeind)):
   
   #print('CME active ',p)
   
   #derive values for ellipse 
   theta=np.arctan((all_apex_f[cmeind[0][p]]**2)*np.tan(all_apex_w[cmeind[0][p]]))
   omega=np.sqrt((np.cos(theta)**2)*(all_apex_f[cmeind[0][p]]**2-1)+1)

   #ellipse values, depending on R and lamda and f, from Moestl et al. 2015 Nat. Comm.
   b=(all_apex_r[cmeind[0][p]]*omega*np.sin(all_apex_w[cmeind[0][p]]))/  ( np.cos(all_apex_w[cmeind[0][p]]-theta)+omega*np.sin(all_apex_w[cmeind[0][p]]))
   a=b/all_apex_f[cmeind[0][p]]
   c=all_apex_r[cmeind[0][p]]-b #center distance of ellipse
   
   #print('a,b,c:',a,b,c)
      
   #ellipse apex and center
   [xapex,yapex]=np.array([np.cos(all_apex_lon[cmeind[0][p]]*np.pi/180),np.sin(all_apex_lon[cmeind[0][p]]*np.pi/180)])*all_apex_r[cmeind[0][p]]
   [xc,yc]=np.array([np.cos(all_apex_lon[cmeind[0][p]]*np.pi/180),np.sin(all_apex_lon[cmeind[0][p]]*np.pi/180)])*c
   
   #convert only apex to show
   #now convert to polar coordinates
   rapex=np.sqrt(xapex**2+yapex**2)
   longapex=np.arctan2(yapex,xapex)
   #print(rapex,longapex*180/np.pi)
   #ax.scatter(longapex,rapex,c='k',s=20)
   
   #rc=np.sqrt(xc**2+yc**2)
   #lc=np.arctan2(yc,xc)
   #print(rc,lc*180/np.pi)
   #ax.scatter(lc,rc,c='r',s=20)
   #point at x=1 y=1
   #r1=np.sqrt(0**2+1**2)
   #l1=np.arctan2(0,1)
   #ax.scatter(l1,r1,c='b',s=50)
   
   #make points on ellipse
   circ_ang = ((np.arange(111)*2-110)*np.pi/180)
   
   xe =  b*np.cos(circ_ang)              #Parameterized equation of ellipse
   ye =  a*np.sin(circ_ang)
     
   #rotation angle
   cosang = np.cos(all_apex_lon[cmeind[0][p]]*np.pi/180)#-np.deg2rad(90))
   sinang = np.sin(all_apex_lon[cmeind[0][p]]*np.pi/180)#-np.deg2rad(90))    
      
   xell = xc + xe*cosang - ye*sinang   	#Rotate to desired position angle
   yell = yc + xe*sinang + ye*cosang
   
   #now convert to polar coordinates
   rell=np.sqrt(xell**2+yell**2)
   longell=np.arctan2(yell,xell)
   
   #plot in correct color
   if all_apex_s[cmeind[0][p]] == 'A':    
    #make alpha dependent on distance to solar equatorial plane
    ax.plot(longell,rell, c='red', alpha=1-abs(all_apex_lat[cmeind[0][p]]/50), lw=1.5) 
   if all_apex_s[cmeind[0][p]] == 'B':
    #ax.plot(longell,rell, c='royalblue', alpha=1-abs(all_apex_lat[cmeind[0][p]]/50), lw=1.5) 
    
    #alpha should depend on colorflag
    if all_apex_flag[cmeind[0][p]] == 0:
       #ax.plot(longell,rell, c='grey', alpha=0.6, lw=1,zorder=1)     
       ax.plot(longell,rell, c='silver', alpha=0.6, lw=1,zorder=1)     
       
       #if all_apex_flag[cmeind[0][p]] ==1:
       #ax.plot(longell,rell, c='silver', alpha=0.8, lw=1,zorder=2)     

    if all_apex_flag[cmeind[0][p]] ==1:
       ax.plot(longell,rell, c='silver', alpha=0.6, lw=1,zorder=1)     
       
    if all_apex_flag[cmeind[0][p]] ==2:
       ax.plot(longell,rell, c='black', alpha=1, lw=1,zorder=3) 
    
   
###############################plot elongation
#  #difference array of current frame time frame_time_num+k to position time frame_time_num
#  elondt=frame_time_num+k-et_time_num_interp
#  #get indices where difference is less than half the time resolution
#  elonind=np.where(abs(elondt) < dayjump / 2.0)
#  
#  #print( 'elonind', cmeind)

#  if np.size(elonind) >0: 
#    #for ElEvoHI2 paper Amerstorfer et al. 2017   
#    ################## add tangent from STEREO-B to ellipse using the time elongation profile
#    #this is the currently active epsilon for the active CME
#    angletox=np.deg2rad(180-et_elon_interp[elonind[0]]-abs(np.rad2deg(pos.stb[1,timeind])))#+np.pi/2
#    tangent_size=1 #AU
#    #make x y coordinates of tangent vector from 0/0
#    vecx1=tangent_size*np.cos(angletox) 
#    vecy1=tangent_size*np.sin(angletox)  
#    stbx=pos.stb[0,timeind]*np.cos(pos.stb[1,timeind])
#    stby=pos.stb[0,timeind]*np.sin(pos.stb[1,timeind])
#    elonx1=stbx+vecx1
#    elony1=stby+vecy1   
#    elonr=np.sqrt(elonx1**2+elony1**2)
#    elonlong=np.arctan2(elony1,elonx1)
#    
#    #end of fit at AU 0.4557, this is h_time_num[26] = 734081.4166666657
#    #before this time plot elongation as straight line
#    if frame_time_num+k <  h_time_num[26]:
#      ax.plot([pos.stb[1,timeind],elonlong], [pos.stb[0,timeind],elonr], c='royalblue', alpha=1, lw=1)
#    else:   #afterwards dashed line
#      ax.plot([pos.stb[1,timeind],elonlong], [pos.stb[0,timeind],elonr], c='royalblue', alpha=1, lw=1,  ls='--')
#    
   
 
 ################## plot positions
 
 #index 1 is longitude, 0 is rdist
 ax.scatter(pos.venus[1,timeind], pos.venus[0,timeind], s=50, c='orange', alpha=1, lw=0, zorder=3)
 ax.scatter(pos.mercury[1,timeind], pos.mercury[0,timeind], s=50, c='dimgrey', alpha=1,lw=0, zorder=3)
 ax.scatter(pos.messenger[1,timeind], pos.messenger[0,timeind], s=25, c='dimgrey',marker='s', alpha=1,lw=0,zorder=3)
 ax.scatter(pos.sta[1,timeind], pos.sta[0,timeind], s=25, c='red', alpha=1,marker='s',lw=0,zorder=3)
 ax.scatter(pos.stb[1,timeind], pos.stb[0,timeind], s=25, c='royalblue', alpha=1,marker='s',lw=0,zorder=3)
 ax.scatter(pos.earth[1,timeind], pos.earth[0,timeind], s=50, c='mediumseagreen', alpha=1,lw=0,zorder=3)
 ax.scatter(pos.mars[1,timeind], pos.mars[0,timeind], s=50, c='orangered', alpha=1,lw=0,zorder=3)
 ax.scatter(pos.msl[1,timeind], pos.msl[0,timeind], s=25, c='magenta', marker='s',alpha=1,lw=0,zorder=3)
 ax.scatter(pos.maven[1,timeind], pos.maven[0,timeind], s=25, c='steelblue',marker='s', alpha=1,lw=0,zorder=3)
 ax.scatter(pos.rosetta[1,timeind], pos.rosetta[0,timeind], s=25, c='black', marker='s', alpha=1,lw=0,zorder=3)
 ax.scatter(pos.ulysses[1,timeind], pos.ulysses[0,timeind], s=25, c='darkolivegreen', marker='s', alpha=1,lw=0,zorder=3)

 #######################  plot ICME detections 
 ######## for each frame time, check active ICMEs looking into ICMECAT:
 for m in range(0,len(icme_start_time_num)):
 
  #calculate difference in arrival_time_num_time to current frame
  icme_diff_to_frame=(frame_time_num+k)-icme_start_time_num[m]
 
  #for all arrival_time_num_times that are later than the current frame, 
  #make them active for fadedays (fading) or infinite (keeping).
  if  icme_diff_to_frame > 0 and icme_diff_to_frame < fadedays:
     #check if this active icme belongs to a spacecraft
     #in1d compares to arrays; true or 1 if m is contained in vexind
      if np.in1d(m,vexind) == 1:
         active_icme_vex[m]=icme_diff_to_frame
     #same for the other spacecraft    
      if np.in1d(m,stbind) == 1:
         active_icme_stb[m]=icme_diff_to_frame
      if np.in1d(m,staind) == 1:
         active_icme_sta[m]=icme_diff_to_frame
      if np.in1d(m,winind) == 1:
         active_icme_win[m]=icme_diff_to_frame
      if np.in1d(m,mesind) == 1:
         active_icme_mes[m]=icme_diff_to_frame
      if np.in1d(m,ulyind) == 1:
         active_icme_uly[m]=icme_diff_to_frame
  else:
     #if no detection, set the index to 0
     active_icme_vex[m]=0
     active_icme_stb[m]=0
     active_icme_sta[m]=0
     active_icme_win[m]=0
     active_icme_mes[m]=0
     active_icme_uly[m]=0
   
 
 #look which ICMEs are active
 active_index_vex=np.where(active_icme_vex > 0)
 active_index_stb=np.where(active_icme_stb > 0)
 active_index_sta=np.where(active_icme_sta > 0)
 active_index_win=np.where(active_icme_win > 0)
 active_index_mes=np.where(active_icme_mes > 0) 
 active_index_uly=np.where(active_icme_uly > 0)
 
  
 #fader style plot alpha dependent on time difference - for this loop over each element:

 for y in range(0,np.size(active_index_vex)):
   z=active_index_vex[0][y] #access elements in tuple that is produced by where
   fadealpha=1-active_icme_vex[z]/(fadedays)  #fadedays is maximum difference in time, and alpha from 0 to 1
   ax.scatter(long[z], rdist[z], s=bmean[z], c='orange', alpha=fadealpha,zorder=4)

 for y in range(0,np.size(active_index_sta)):
   z=active_index_sta[0][y]
   fadealpha=1-active_icme_sta[z]/(fadedays)  #30 days is maximum difference in time, and alpha from 0 to 1
   ax.scatter(long[z], rdist[z], s=bmean[z], c='red', alpha=fadealpha,zorder=4)

 for y in range(0,np.size(active_index_stb)):
   z=active_index_stb[0][y]
   fadealpha=1-active_icme_stb[z]/(fadedays)  #30 days is maximum difference in time, and alpha from 0 to 1
   ax.scatter(long[z], rdist[z], s=bmean[z], c='royalblue', alpha=fadealpha,zorder=4)

 for y in range(0,np.size(active_index_win)):
   z=active_index_win[0][y]
   fadealpha=1-active_icme_win[z]/(fadedays)  #30 days is maximum difference in time, and alpha from 0 to 1
   ax.scatter(long[z], rdist[z], s=bmean[z], c='mediumseagreen', alpha=fadealpha,zorder=4)

 for y in range(0,np.size(active_index_mes)): 
   z=active_index_mes[0][y]
   fadealpha=1-active_icme_mes[z]/(fadedays)  #30 days is maximum difference in time, and alpha from 0 to 1
   ax.scatter(long[z], rdist[z], s=bmean[z], c='dimgrey', alpha=fadealpha,zorder=4)

 for y in range(0,np.size(active_index_uly)): 
   z=active_index_uly[0][y]
   fadealpha=1-active_icme_uly[z]/(fadedays)  #30 days is maximum difference in time, and alpha from 0 to 1
   ax.scatter(long[z], rdist[z], s=bmean[z], c='darkolivegreen', alpha=fadealpha,zorder=4)


 ###################### legend and additional text
 
 plt.suptitle('ElEvoHI ensemble simulation ')	
 
 #Sun
 ax.scatter(0,0,s=100,c='yellow',alpha=0.8, edgecolors='yellow')
 plt.figtext(0.51,0.5,'Sun', fontsize=10, ha='center')
 
 #Earth
 plt.figtext(0.51,0.28,'Earth', fontsize=10, ha='center')
 
 if HEEQ == 1: plt.figtext(0.525,0.0735,'HEEQ longitude', fontsize=10, ha='left')
 if HEE == 1: plt.figtext(0.525,0.0735,'HEE longitude', fontsize=10, ha='left')
 #	plt.figtext(0.64,0.213,'AU', fontsize=10, ha='center')

 plt.figtext(0.1-0.02,0.02,'Mercury', color='dimgrey', ha='center',fontsize=labelfontsize)
 plt.figtext(0.2-0.02,0.02,'MESSENGER', color='dimgrey', ha='center', fontsize=labelfontsize)
 plt.figtext(0.3-0.02	,0.02,'Venus', color='orange', ha='center',fontsize=labelfontsize)
 plt.figtext(0.4-0.02,0.02,'STEREO-A', color='red', ha='center',fontsize=labelfontsize)
 plt.figtext(0.53-0.02,0.02,'STEREO-B', color='royalblue', ha='center',fontsize=labelfontsize)
 plt.figtext(0.62-0.02,0.02,'Earth', color='mediumseagreen', ha='center',fontsize=labelfontsize)
 plt.figtext(0.68-0.02,0.02,'Mars', color='orangered', ha='center',fontsize=labelfontsize)
 plt.figtext(0.78-0.02,0.02,'Maven', color='steelblue', ha='center', fontsize=labelfontsize)
 plt.figtext(0.73-0.02,0.02,'MSL', color='magenta', ha='center', fontsize=labelfontsize)
 plt.figtext(0.84-0.02,0.02,'Rosetta', color='black', ha='center', fontsize=labelfontsize)
 plt.figtext(0.90-0.02,0.02,'Ulysses', color='darkolivegreen', ha='center', fontsize=labelfontsize)

 #add legend for bmean
 bleg=np.array([10,50,100])*bscale
 blegstr=['10 nT','50','100']

 blegr=np.zeros(len(bleg))+1.6
 blegt=np.radians(range(170,195,10))
 ax.scatter(blegt, blegr,s=bleg,c='violet', edgecolor='violet')

 for p in range(0,len(bleg)):
   ax.annotate(blegstr[p],xy=(blegt[p],blegr[p]-0.2), ha='center', va='center', fontsize=8)
 
 #set axes
 plt.thetagrids(range(0,360,45),(u'0\u00b0',u'45\u00b0',u'90\u00b0',u'135\u00b0',u'+/- 180\u00b0',u'- 135\u00b0',u'- 90\u00b0',u'- 45\u00b0'), fmt='%d')#, frac = 1.05)
 ax.set_theta_zero_location('S')
 plt.rgrids((0.25,0.5,0.75, 1.0,1.25, 1.5, 1.75, 2.0),('0.25','0.5','0.75','1.0','1.25','1.5','1.75','2.0 AU'),angle=150, fontsize=8)
 ax.set_ylim(0, 2.1)
 
 #plot text for date extra so it does not move 
 #year
 plt.figtext(0.47,0.85,frame_time_str[0:4], fontsize=13, ha='center')
 #month
 plt.figtext(0.51,0.85,frame_time_str[5:7], fontsize=13, ha='center')
 #day
 plt.figtext(0.54,0.85,frame_time_str[8:10], fontsize=13, ha='center')
  #hours
 plt.figtext(0.57,0.85,frame_time_str[11:13], fontsize=13, ha='center')

 #signature
 plt.figtext(0.95,0.01/2,r'$C. M\ddot{o}stl, T. Amerstorfer$', fontsize=4, ha='center')



 ###################### save frame
 plt.savefig(current_event_dir+'/frames/elevo_'+framestr+'.png', dpi=300)
 #clears plot window
 plt.clf()

############ end of loop


################################################### MAKE MOVIE

#convert to jpg
os.system(os.getcwd()+'/ffmpeg -i '+current_event_dir+'frames/elevo_%04d.png '+current_event_dir+'frames/elevo_%04d.jpg ')
#make mp4
os.system(os.getcwd()+'/ffmpeg -r 20 -i '+current_event_dir+'frames/elevo_%04d.jpg -b:v 5000k -r 20 movies/'+current_event+'_ensemble_movie.mp4 -y')
#make gif
os.system(os.getcwd()+'/ffmpeg -r 20 -i movies/'+current_event+'_ensemble_movie.mp4 -b:v 5000k -r 20 movies/'+current_event+'_ensemble_final.gif -y')


plt.close('all')

print( 'Made movie.')
print( 'End ElEvoHI animation program.')


########################### MIT license


#Copyright 2018 Mag. Dr. Christian Moestl
#Permission is hereby granted, free of charge, to any person obtaining a copy of this 
#software and associated documentation files (the "Software"), to deal in the Software 
#without restriction, including without limitation the rights to use, copy, modify, merge, 
#publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons
#to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included 
#in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
#PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
#FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
#TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
#OTHER DEALINGS IN THE SOFTWARE.




