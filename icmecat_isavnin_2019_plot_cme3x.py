#animation ICMECAT rose scatter plot with HI CME circles and in situ data	

from scipy import stats
import scipy.io
from matplotlib import cm
import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import sunpy.time
import time
import pickle
import seaborn as sns
import math


#for reading catalogues  
def getcat(filename):
  print('reading CAT '+filename)
  cat=scipy.io.readsav(filename, verbose='true')  
  print('done reading CAT')
  return cat  
  
  
def getpositions(filename):  
  print( 'reading positions in '+filename)
  pos=scipy.io.readsav(filename, verbose='true')  
  print( 'done reading positions')
  return pos

  
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
   
   
   
   print(time_str[j])
   year=int(time_str[j][0:4])
   month=int(time_str[5:7])
   #convert time to sunpy friendly time and to matplotlibdatetime
   #only for valid times so 9999 in year is not converted
   #pdb.set_trace()
   #if year < 2100:
   
   #time_num[j]=mdates.date2num(datetime.datetime(year,month,day,hour)
   
   # time_num[j]=mdates.date2num(sunpy.time.parse_time(time_str[j]))
   j=j+1  
   #the date format in matplotlib is e.g. 735202.67569444
   #this is time in days since 0001-01-01 UTC, plus 1.
   
   #return time_num which is already an array and convert the list of strings to an array
  return time_num, np.array(time_str)


def decode_array(bytearrin):
 #for decoding the strings from the IDL .sav file to a list of python strings, not bytes 
 #make list of python lists with arbitrary length
 bytearrout= ['' for x in range(len(bytearrin))]
 for i in range(0,len(bytearrin)-1):
  bytearrout[i]=bytearrin[i].decode()
 #has to be np array so to be used with numpy "where"
 bytearrout=np.array(bytearrout)
 return bytearrout  

  
  
def IDL_time_to_num(time_in):  
 #convert IDL time to matplotlib datetime
 time_num=np.zeros(np.size(time_in))
 for ii in np.arange(0,np.size(time_in)):
   time_num[ii]=mdates.date2num(sunpy.time.parse_time(time_in[ii]))
   
 return time_num 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    

######################################################
#main program


plt.close('all')



sns.set_context("talk")     
sns.set_style("darkgrid")  
################## CONTROLS

#how much time is between frames
dayjump=0.25

#either keep or fade detections
fade=1
keep=0

#if keep is selected, the alpha for plotting each dot
keepalpha=0.7

#how long an ARRIVAL stays visible in fade mode
fadedays=30

#how big the circles are on the plot
bscale=4

#half width of the circles
lamda=30



################################


print( 'start icmecat animation program.')

#get ICMECAT
filename_icmecat='../catpy/ALLCATS/HELCATS_ICMECAT_v10_SCEQ.sav'
i=getcat(filename_icmecat)

#get parameters
bmean=i.icmecat['MO_BMEAN']*bscale #bscale makes circles larger in movie
long=i.icmecat['SC_LONG_HEEQ']*np.pi/180 #heeq longitude converted to radians
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


##################################### read in situ

print( 'read MESSENGER')
#get insitu data
mes= pickle.load( open( "../catpy/DATACAT/MES_2007to2015_SCEQ_removed.p", "rb" ) )
#time conversion
#mes_time=IDL_time_to_num(mes.time)
print( 'read MESSENGER done.')



print ('read VEX')
#get insitu data
vex= pickle.load( open( "../catpy/DATACAT/VEX_2007to2014_SCEQ_removed.p", "rb" ) )
#time conversion
#vex_time=IDL_time_to_num(vex.time)
print( 'read VEX done.')



print( 'read Wind')
#get insitu data
wind= pickle.load( open( "../catpy/DATACAT/WIND_2007to2016_HEEQ.p", "rb" ) )
#time conversion
#wind_time=IDL_time_to_num(wind.time)
print( 'read Wind done.')




print( 'read STEREO-A')
#get insitu data
sta= pickle.load( open( "../catpy/DATACAT/STA_2007to2015_SCEQ.p", "rb" ) )
#time conversion
#sta_time=IDL_time_to_num(sta.time)
print( 'read STA done.')




print( 'read STEREO-B')
#get insitu data
stb= pickle.load( open( "../catpy/DATACAT/STB_2007to2014_SCEQ.p", "rb" ) )

#time conversion
#stb_time=IDL_time_to_num(stb.time)
print( 'read STB done.')

#save times
#pickle.dump([vex_time,wind_time,sta_time,stb_time,mes_time], open( "DATACAT/Insitu_times_mdates_2.p", "wb" ) )

#quicker when just reloading times
[vex_time,wind_time,sta_time,stb_time,mes_time]=pickle.load( open( "../catpy/DATACAT/Insitu_times_mdates_2.p", "rb" ) )
#print 'loaded in situ times'
######################################



#initiate plot
plt.figure(1, figsize=(6, 6), dpi=100, facecolor='w', edgecolor='w')


#full  movie April 2014 Jan 1 until end of November 2014
frame_time_num=mdates.date2num(sunpy.time.parse_time('2011-Jun-4').datetime)

dayjump=1

################################### plot over all frames
for k in np.arange(0,1,dayjump):  
#3169 is time in days
  
  
   
 start=time.time()
 #to current frame time, the days need to be added, so +k is done 
 #save frame time as string to write on plot
 frame_time_str=str(mdates.num2date(frame_time_num+k))
 print( 'current frame_time_num+k', frame_time_str)
 

 ############# 5 in situ data plots 

 plotstartdate=mdates.num2date(frame_time_num+k)
 plotenddate=mdates.num2date(frame_time_num+k+4)
 
 
 
 #slicing
 
 #take only those indices where the difference to frame_time_num+k is less than 3
 mes_ind_plot=np.where(abs(mes_time-(frame_time_num+k)) < 4)
 vex_ind_plot=np.where(abs(vex_time-(frame_time_num+k)) < 4)
 stb_ind_plot=np.where(abs(stb_time-(frame_time_num+k)) < 4)
 sta_ind_plot=np.where(abs(sta_time-(frame_time_num+k)) < 4)
 wind_ind_plot=np.where(abs(wind_time-(frame_time_num+k)) < 4)
  
 
 lineweit=0.2
 
 #rows - columns
 #MESSENGER
 ax2 = plt.subplot2grid((3,1), (0, 0))
 ax2.plot_date(mes_time[mes_ind_plot],mes.btot[mes_ind_plot],'-k', lw=lineweit)
 ax2.plot_date(mes_time[mes_ind_plot],mes.bx[mes_ind_plot], '-r',lw=lineweit)
 ax2.plot_date(mes_time[mes_ind_plot],mes.by[mes_ind_plot],'-g',lw=lineweit)
 ax2.plot_date(mes_time[mes_ind_plot],mes.bz[mes_ind_plot],'-b',lw=lineweit)
 #current time
 plt.yticks(fontsize=10)
 plt.ylabel('B SCEQ [nT]', fontsize=10)
 #ax2.plot_date([mdates.num2date(frame_time_num+k),mdates.num2date(frame_time_num+k)], [-120,120],'-k', lw=0.5, alpha=0.8)
 plt.xlim((plotstartdate, plotenddate))
 plt.ylim((-300,300))
 locs, labels = plt.xticks() 
 plt.xticks(locs,labels=[])   
 plt.title('Mercury',color='orange', fontsize=10)  


 
 #VEX
 ax3 = plt.subplot2grid((3,1), (1, 0))
 ax3.plot_date(vex_time[vex_ind_plot],vex.btot[vex_ind_plot],'-k', lw=lineweit)
 ax3.plot_date(vex_time[vex_ind_plot],vex.bx[vex_ind_plot], '-r',lw=lineweit)
 ax3.plot_date(vex_time[vex_ind_plot],vex.by[vex_ind_plot],'-g',lw=lineweit)
 ax3.plot_date(vex_time[vex_ind_plot],vex.bz[vex_ind_plot],'-b',lw=lineweit)
  
 plt.yticks(fontsize=10)
 plt.ylabel('B SCEQ [nT]', fontsize=10)
 #ax3.plot_date([mdates.num2date(frame_time_num+k),mdates.num2date(frame_time_num+k)], [-50,50],'-k', lw=0.5, alpha=0.8)
 plt.xlim((plotstartdate, plotenddate))
 plt.ylim((-120, 120))
 locs, labels = plt.xticks() 
 plt.xticks(locs,labels=[])  
 plt.title('Venus',color='black', fontsize=10)  
 

 #STA
 ax5 = plt.subplot2grid((3,1), (2, 0))
 ax5.plot_date(sta_time[sta_ind_plot],sta.btot[sta_ind_plot],'-k', lw=lineweit)
 ax5.plot_date(sta_time[sta_ind_plot],sta.bx[sta_ind_plot], '-r',lw=lineweit)
 ax5.plot_date(sta_time[sta_ind_plot],sta.by[sta_ind_plot],'-g',lw=lineweit)
 ax5.plot_date(sta_time[sta_ind_plot],sta.bz[sta_ind_plot],'-b',lw=lineweit)

 plt.ylabel('B SCEQ [nT]', fontsize=10)
 #ax5.plot_date([mdates.num2date(frame_time_num+k),mdates.num2date(frame_time_num+k)], [-50,50],'-k', lw=0.5, alpha=0.8)
 plt.xlim((plotstartdate, plotenddate))
 plt.ylim((-30, 30))
 plt.title('STEREO-A',color='red', fontsize=10)  
  
 myformat = mdates.DateFormatter('%h-%d %H:00') 
 ax5.xaxis.set_major_formatter(myformat) 
 
 ax5.tick_params(axis='x', rotation=15)
 
 plt.xlabel('time [UTC] in year 2011' , fontsize=10)
  
 plt.yticks(fontsize=10)
 plt.xticks(fontsize=10)
 '''
 #Earth
 ax4 = plt.subplot2grid((5,1), (3, 0))
 ax4.plot_date(wind_time[wind_ind_plot],wind.btot[wind_ind_plot],'-k', lw=0.3)
 ax4.plot_date(wind_time[wind_ind_plot],wind.bx[wind_ind_plot], '-r',lw=0.3)
 ax4.plot_date(wind_time[wind_ind_plot],wind.by[wind_ind_plot],'-g',lw=0.3)
 ax4.plot_date(wind_time[wind_ind_plot],wind.bz[wind_ind_plot],'-b',lw=0.3)
 
 plt.yticks(fontsize=9)
 plt.ylabel('B SCEQ [nT]', fontsize=9)
 #ax4.plot_date([mdates.num2date(frame_time_num+k),mdates.num2date(frame_time_num+k)], [-50,50],'-k', lw=0.5, alpha=0.8)
 plt.xlim((plotstartdate, plotenddate))
 plt.ylim((-35, 35))
 locs, labels = plt.xticks() 
 plt.xticks(locs,labels=[])   


 #STB
 ax6 = plt.subplot2grid((5,1), (4, 0))
 ax6.plot_date(stb_time[stb_ind_plot],stb.btot[stb_ind_plot],'-k', lw=0.3)
 ax6.plot_date(stb_time[stb_ind_plot],stb.bx[stb_ind_plot], '-r',lw=0.3)
 ax6.plot_date(stb_time[stb_ind_plot],stb.by[stb_ind_plot],'-g',lw=0.3)
 ax6.plot_date(stb_time[stb_ind_plot],stb.bz[stb_ind_plot],'-b',lw=0.3)  
 plt.xlim((plotstartdate, plotenddate))
   
 myformat = mdates.DateFormatter('%d-%Hh')
 ax6.xaxis.set_major_formatter(myformat)  
 plt.yticks(fontsize=8)
 plt.ylabel('B SCEQ [nT]', fontsize=9)
 #ax6.plot_date([mdates.num2date(frame_time_num+k),mdates.num2date(frame_time_num+k)], [-50,50],'-k', lw=0.5, alpha=0.8)
 plt.ylim((-35, 35))
 plt.xticks(fontsize=10)
 '''
 
 #labeling of spacecraft and longitude in HEEQ
 #plt.figtext(0.55,0.96,'Mercury',color='orange', fontsize=10, ha='center')  

 #plt.figtext(0.55,0.64,'Venus',color='black', fontsize=10, ha='center')  
 
 #plt.figtext(0.55,0.32,'STEREO-A',color='red', fontsize=10, ha='center')  

 
 #plt.figtext(0.92,0.82-0.165,'VEX',color='orange', fontsize=10, ha='left')
 #plt.figtext(0.92,0.82-0.165*2,'STEREO-A',color='red', fontsize=10, ha='left')

 #plt.figtext(0.92,0.82-0.165*3,'Wind',color='mediumseagreen', fontsize=10, ha='left')


 #plt.figtext(0.92,0.82-0.165*4,'STEREO-B',color='royalblue', fontsize=10, ha='left')

 #labeling in situ components
 plt.figtext(0.75,0.87,'Bx',color='red', fontsize=10, ha='left')
 plt.figtext(0.8,0.87,'By',color='green', fontsize=10, ha='left')
 plt.figtext(0.85,0.87,'Bz',color='blue', fontsize=10, ha='left')

 
 plt.tight_layout()
 
 #save figure for frame - this starts with zero at the start time
 framestr = '%04i' % (k*4)  
 #framenr=framenr+1
 print( 'frame nr.', framestr) 
 #plt.show()
 plt.savefig('events/june_2011/icmecat_'+framestr+'.eps',  dpi=300)


 end=time.time()
 print( 'took time in seconds:', (end-start) ,'for this frame')

 plt.show()     


