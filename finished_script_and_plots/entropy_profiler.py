
# INPUT DATA HERE
"""
252:
ID:  11 , Mass:  12.754855574221434 , Contamination:  0.7026546540631837
ID:  12 , Mass:  12.758119784129255 , Contamination:  0.0
ID:  13 , Mass:  12.759746139553343 , Contamination:  0.0061259334279152285

3224:
ID:  1 , Mass:  12.907166210396783 , Contamination:  0.01475454746052201
ID:  2 , Mass:  12.978288210012412 , Contamination:  0.008205398211531788
ID:  3 , Mass:  12.976965557061892 , Contamination:  0.0
"""

# Simulation data
zoom_ID = 3224 #ie. 252, 3224, etc...
halo_ID_arr = [3] # caesar halo ID

# Observational data
low_mass_range = 12.6
high_mass_range = 13.1

# Plotting/saving
#savefilename = 'Snap' + str(zoom_ID) + '_Halo_contamination.png' # .png name to save plot under
savefilename = 'RECENT_TEST'
thetitle = 'Entropy - L1 Hyenas Halo Contamination Comparision' # title of graph

#---------------------------------------------------------------------------------------------------------------------
import sys
from pathlib import Path
# sys.path.append(str(Path(__file__).absolute().parent))
sys.path.append('../')

#import yt
import numpy as np
import pynbody as pnb
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import mpl_toolkits.axes_grid1 as axgrid
import scipy.stats as scist
import h5py
from astropy import constants as const
import os.path
# import home.b.babul.aspadawe.scripts.oppenheimer21.group_data
# import group_data

# import numpy as np
# import matplotlib as mpl
# import matplotlib.pyplot as plt
# import mpl_toolkits.axes_grid1 as axgrid
# import os.path

import scaling_relation_functions as scale
# from scaling_relation_functions import h_sun09, h_simba
# from oppenheimer21.groups_profiles_all import run_plotting_script as opp21_plot

mu_e = 1.14
mu = 0.59
mp = 1.6726e-24
G = 6.672e-08
kboltz = 1.38066e-16
kboltz_v2 = 8.617e-8 # keV/K
cmperkpc = 3.086e+21
cmperMpc = 3.086e+24
keVperK = 8.616e-08
gperMsol = 1.989e+33
cmperkpc = 3.086e+21
cmperpc = 3.086e+18

# observational data located in:  /project/b/babul/aspadawe/data/groups-and-clusters
#home_name = "/home/ariboppe/groups_review/"
home_name = "/project/b/babul/aspadawe/data/groups-and-clusters/Oppenheimer+21/"

# group_data.py is for reading and loading in the data


## log(M500*E(z)/Msun) - log(Tspec [keV]) scaling relation for the 512/50 Simba-C box
# tspecs = np.linspace(-0.6,0.48, 20)
# y=1.935*tspecs+13.32
# y1=1.935*tspecs+13.62
# y2=1.935*tspecs+13.02

# def M500_Tspec_simbac(Tspec):
#     ## tspec -> log(Tspec [keV])
#     ## Returns: M500 -> log(M500/Msun)
#     return 1.935*Tspec+13.32

# def M500_T500_sun09(kT500, h):
#     ## kT500 in keV
#     ## M500 in Msun
    
#     M3 = 1.27e14 * h**(-1)  # Msun
#     alpha = 1.67
#     r3 = 0.602 * h**(-1)  # Mpc
    
#     return M3 * (kT500/3.)**alpha

# def R500_T500_sun09(kT500, h):
#     ## kT500 in keV
#     ## R500 in Mpc
    
#     M3 = 1.27e14 * h**(-1)  # Msun
#     alpha = 1.67
#     r3 = 0.602 * h**(-1)  # Mpc
    
#     return r3 * (kT500/3.)**(alpha/3.)
    





#home_name = "/home/ariboppe/groups_review/"
#home_name = "~/project/data/groups-and-clusters/Oppenheimer+21/"
# home_name = "/home/b/babul/aspadawe/project/data/groups-and-clusters/Oppenheimer+21/"
# cm = plt.get_cmap('plasma')
# cm = plt.get_cmap('viridis')


#plot_individual_observations = 0


def load_Sun2009(group_ID, property_base):
    fname = home_name + "Sun2009/" + property_base + "/" + group_ID + "_" + property_base
    #print("fname= ", fname)
    if(os.path.isfile(fname)):
        R,y = np.loadtxt(fname,usecols=(0,1),unpack=True)
        if property_base == 'temperature':
            y /= kboltz_v2
#         print("y[-1]= ", y[-1])
    else:
#         R = [-99,-99]
#         y = [-99,-99]
#         R = [np.nan, np.nan]
#         y = [np.nan, np.nan]
        R = []
        y = []

    return(R,y)

def load_Lovisari2015(obs_Lovisari_ID,property):
    
    fname = home_name + "Lovisari2015/" + obs_Lovisari_ID + "_sb.dat"
    
    if(os.path.isfile(fname)):
        R_kpc,y,yerr = np.loadtxt(fname,usecols=(2,4,5),unpack=True)
    else:
#         R_kpc = [-99,-99]
#         y = [-99,-99]
#         yerr = [0,0]
#         R = [np.nan, np.nan]
#         y = [np.nan, np.nan]
#         yerr = [np.nan, np.nan]
        R_kpc = []  #??
        y = []
        yerr = []

    return(R_kpc,y,yerr)

def load_Sun2009_1550_Z():

    fname = home_name + "Sun2009/N1550_Z_3d"

    R_500_N1550_kpc = 465 

    R_kpc,R_kpc_err,Z_3d,Z_3d_upper,Z_3d_lower = np.loadtxt(fname,usecols=(0,1,2,3,4),unpack=True)

    lR_500 = np.log10(R_kpc/R_500_N1550_kpc)
    
    return(lR_500,Z_3d,Z_3d_upper,Z_3d_lower)
    

def load_and_plot_Sun2009(property, ax_coll, cbar_lo, cbar_hi, lhalo_lo, lhalo_hi, 
                          plot_individual_observations=True, plot_x_log=True, plot_y_log=True, 
                          use_cbar=False, log_cbar=False, cm=plt.get_cmap('plasma'), obs_plot_kwargs=None):

    Sun_lR_coll = []
    Sun_ly_coll = []
    Sun_lM500_coll = []
    
    Sungroups_list = home_name + "Sun2009/Sun_groups.cat.sort"

    Sungroups = open(Sungroups_list,"r").readlines()
    nSungroups = len(Sungroups)

    Sun_label = 0 

    for i in range(nSungroups):
        Sun_ID = Sungroups[i].split()[0]
        Sun_T500 = float(Sungroups[i].split()[1]) ## keV
        Sun_R500 = float(Sungroups[i].split()[2])
        Sun_M500_13 = float(Sungroups[i].split()[3])
        if(Sun_M500_13<0):
            Sun_M500_13 = 3e+13*(Sun_T500/1.0)**(1.5)/1e+13
        
#         Sun_M500_13 = scale.M500_Tspec_simbac_N1024L100_renier(pnb.array.SimArray(Sun_T500, units='keV'))
        Sun_lM500 = np.log10(Sun_M500_13)+13
        
        if (Sun_lM500 >= lhalo_lo and Sun_lM500 <= lhalo_hi):
        
            if(property=="entropy_norm"):
                lR,ly = load_Sun2009(Sun_ID,"K")
            if(property=="density"):
                lR,ly = load_Sun2009(Sun_ID,"den")
                lR = lR + np.log10(0.471) # In R2500, got to get to R500
            if(property=="temperature"):
                lR,y = load_Sun2009(Sun_ID,"T_3d")
    #             ly = np.log10(y) + np.log10(Sun_T500)
                ly = np.log10(y) + np.log10(Sun_T500) - np.log10(kboltz_v2)  # put in K
    #             ly = np.log10(y) #- np.log10(kboltz_v2)  # put in K, and keep in T/T500 (instead of just T [K])
            if(property=="ktemperature"):
                lR,y = load_Sun2009(Sun_ID,"T_3d")
                ly = np.log10(y) + np.log10(Sun_T500)  # in keV
#                 ly = np.log10(y) + np.log10(Sun_T500) - np.log10(kboltz_v2)  # put in K
    #             ly = np.log10(y) #- np.log10(kboltz_v2)  # put in K, and keep in T/T500 (instead of just T [K])
            if(property=="temperature_norm"):
                lR,y = load_Sun2009(Sun_ID,"T_3d")
    #             ly = np.log10(y) + np.log10(Sun_T500)
    #             ly = np.log10(y) + np.log10(Sun_T500) - np.log10(kboltz_v2)  # put in K
                ly = np.log10(y) #- np.log10(kboltz_v2)  # put in K, and keep in T/T500 (instead of just T [K])
            if((property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
                lR,ly = load_Sun2009(Sun_ID,"den")
                lR = lR + np.log10(0.471) # In R2500, got to get to R500
                R = 10**lR*Sun_R500
                if(property=="m100kpc"):
                    m100kpc = IntegrateDensity(R,ly,100)/gperMsol
                    ly = np.log10(m100kpc)
                if(property=="m200kpc"):
                    m200kpc = IntegrateDensity(R,ly,200)/gperMsol
                    ly = np.log10(m200kpc)
                if(property=="m400kpc"):
                    m400kpc = IntegrateDensity(R,ly,400)/gperMsol
                    ly = np.log10(m400kpc)
                    
                    
            if not plot_x_log:
                lR = 10**np.array(lR)
            if not plot_y_log:
                ly = 10**np.array(ly)
            

            if((property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
                ax_all.plot(np.log10(Sun_M500_13*1e+13),ly,color="black",marker="x",zorder=0)
                if(i==0):
                    ax_coll.plot(np.log10(Sun_M500_13*1e+13),ly,color="black",marker="x",ls="",ms=12,
                                 zorder=obs_plot_kwargs['Oppenheimer+21']['zorder'],label=" 219")
                else:
                    ax_coll.plot(np.log10(Sun_M500_13*1e+13),ly,color="black",marker="x",ls="",ms=12,
                                 zorder=obs_plot_kwargs['Oppenheimer+21']['zorder'])

            else:

                Sun_lR_coll.extend(lR)
                Sun_ly_coll.extend(ly)
                Sun_lM500_coll.extend([Sun_lM500])
                
                if use_cbar and obs_plot_kwargs['Oppenheimer+21']['use_cbar']:
                    if log_cbar:
#                         sun_color = cm((Sun_lM500-np.log10(cbar_lo))/(np.log10(cbar_hi)-np.log10(cbar_lo)))
                        sun_color = cm((Sun_lM500-cbar_lo)/(cbar_hi-cbar_lo))
                    else:
                        sun_color = cm((10**Sun_lM500-cbar_lo)/(cbar_hi-cbar_lo))
                else:
                    sun_color = obs_plot_kwargs['Oppenheimer+21']['sun_lovisari_color']
                    
                
                if plot_individual_observations and not obs_plot_kwargs['Oppenheimer+21']['plot_sun_median']:
                    if(Sun_label == 0):
                        if((Sun_lM500 >= lhalo_lo) & (Sun_lM500 <= lhalo_hi)):  # Unnecessary?
                            #if(multi == 2 ):  # BDO 5/13/21, trying to replace inferno colorbar.   
                            ax_coll.plot(lR,ly,color=sun_color,lw=2,ls=":",
                                         alpha=obs_plot_kwargs['Oppenheimer+21']['sun_lovisari_alpha'],
                                         zorder=obs_plot_kwargs['Oppenheimer+21']['zorder'],label="Sun+2009")
                            #else:
                            #    ax_coll.plot(lR,ly,color='k',lw=1,ls="--",zorder=0,label="Sun+ (2009)")
                            Sun_label = 1
                    else:
                        if((Sun_lM500 >= lhalo_lo) & (Sun_lM500 <= lhalo_hi)):  # Unnecessary?
                            #if(multi == 2):# BDO 5/13/21, trying to replace viridis colorbar.
                            ax_coll.plot(lR,ly,color=sun_color,lw=2,ls=":",
                                         alpha=obs_plot_kwargs['Oppenheimer+21']['sun_lovisari_alpha'],
                                         zorder=obs_plot_kwargs['Oppenheimer+21']['zorder'])
                            #else:
                            #    ax_coll.plot(lR,ly,color='k',lw=1,ls="--",zorder=0)

                elif plot_individual_observations and obs_plot_kwargs['Oppenheimer+21']['plot_sun_median']:
                    ax_coll.plot(lR,ly,color=sun_color,lw=2,ls=":",
                                 alpha=obs_plot_kwargs['Oppenheimer+21']['sun_lovisari_alpha'],
                                 zorder=obs_plot_kwargs['Oppenheimer+21']['zorder'])

    if (obs_plot_kwargs['Oppenheimer+21']['plot_sun_median'] and len(Sun_lR_coll)>0):
        
#         if not plot_x_log:
#             Sun_lR_coll = 10**Sun_lR_coll
#         if not plot_y_log:
#             Sun_ly_coll = 10**Sun_ly_coll
        
        Sun_lR,Sun_lRlo,Sun_lRhi,Sun_ly,Sun_lylo,Sun_lyhi=flexbin_med_and_onesigmaspread(np.asarray(Sun_lR_coll),np.asarray(Sun_ly_coll),15)
#         print("Sun_lR,Sun_lRlo,Sun_lRhi,Sun_ly,Sun_lylo,Sun_lyhi,Median_M500=",Sun_lR,Sun_lRlo,Sun_lRhi,Sun_ly,Sun_lylo,Sun_lyhi,np.median(np.asarray(Sun_lM500_coll)))
#         Sun_fractional_color = (np.median(10**np.asarray(Sun_lM500_coll))-cbar_lo)/(cbar_hi-cbar_lo)
    
        if use_cbar and obs_plot_kwargs['Oppenheimer+21']['use_cbar']:
            if log_cbar:
#                 Sun_fractional_color = cm((np.median(np.asarray(Sun_lM500_coll))-np.log10(cbar_lo))/(np.log10(cbar_hi)-np.log10(cbar_lo)))
                Sun_fractional_color = cm((np.median(np.asarray(Sun_lM500_coll))-cbar_lo)/(cbar_hi-cbar_lo))
            else:
                Sun_fractional_color = cm((np.median(10**np.asarray(Sun_lM500_coll))-cbar_lo)/(cbar_hi-cbar_lo))
        else:
            Sun_fractional_color = obs_plot_kwargs['Oppenheimer+21']['sun_lovisari_color']
    
        ax_coll.plot(Sun_lR,Sun_ly,color=Sun_fractional_color,lw=3,ls=":",
                     zorder=obs_plot_kwargs['Oppenheimer+21']['zorder']+5,label="Sun+2009")
        ax_coll.fill_between(Sun_lR, Sun_lyhi, Sun_lylo, 
                             color=Sun_fractional_color,
                             alpha=obs_plot_kwargs['Oppenheimer+21']['sun_lovisari_alpha'], 
                             zorder=obs_plot_kwargs['Oppenheimer+21']['zorder'])
#         ax_coll.plot(Sun_lR,Sun_ly,color=cm(Sun_fractional_color),lw=4,ls=":",zorder=0,label="Sun+ (2009)")
#         for si in range(len(Sun_lR)):
#             print("Sun2009:", Sun_lR[si], Sun_ly[si], Sun_lylo[si], Sun_lyhi[si])


def load_and_plot_Lovisari_Sun_groups(property, ax_coll, cbar_lo, cbar_hi, lhalo_lo, lhalo_hi, 
                                      plot_x_log=True, plot_y_log=True, 
                                      use_cbar=False, log_cbar=False, cm=plt.get_cmap('plasma'), obs_plot_kwargs=None):
#     print('hi 2')

    obs_lR_coll = []
    obs_ly_coll = []
    obs_lM500_coll = []

    obsgroups_list = home_name + "Lovisari2015/Sun_Lovisari_groups.txt"
    
    obsgroups = open(obsgroups_list,"r").readlines()
    nobsgroups = len(obsgroups)

    obs_label = 0 
    
#     print('hi')

    for i in range(nobsgroups):
        obs_Sun_ID = obsgroups[i].split()[0]
        obs_Lovisari_ID = obsgroups[i].split()[1]
        obs_T500 = float(obsgroups[i].split()[2])
        obs_R500 = float(obsgroups[i].split()[3])
#         obs_M500_13 = float(obsgroups[i].split()[4])
#         if(obs_M500_13<0):
#             obs_M500_13 = 3e+13*(obs_T500/1.0)**(1.5)/1e+13
        
        obs_M500_13 = scale.M500_Tspec_simbac_N1024L100_renier(pnb.array.SimArray(obs_T500, units='keV'))
        obs_lM500 = np.log10(obs_M500_13)+13
        
        if (obs_lM500 >= lhalo_lo and obs_lM500 <= lhalo_hi):
        
            if(property=="XraySB_XMMcounts"):
                print("obs_Lovisari_ID= ", obs_Lovisari_ID)
                R_kpc,y,yerr = load_Lovisari2015(obs_Lovisari_ID,"XraySB_XMMcounts")
                lR = np.log10(R_kpc/obs_R500)
                ly = np.log10(y)
                lyerr = np.log10(yerr)
            if(property=="entropy_norm"):
                lR,ly = load_Sun2009(obs_Sun_ID,"K")
            if(property=="density"):
    #             print('something')
                lR,ly = load_Sun2009(obs_Sun_ID,"den")
                lR = lR + np.log10(0.471) # In R2500, got to get to R500
            if(property=="temperature"):
                lR,y = load_Sun2009(obs_Sun_ID,"T_3d")
    #             ly = np.log10(y) + np.log10(obs_T500)
                ly = np.log10(y) + np.log10(obs_T500) - np.log10(kboltz_v2)  # put in K
    #             ly = np.log10(y) #- np.log10(kboltz_v2)  # put in K, and keep in T/T500 (instead of just T [K])
            if(property=="ktemperature"):
                lR,y = load_Sun2009(obs_Sun_ID,"T_3d")
                ly = np.log10(y) + np.log10(obs_T500)  # keV
#                 ly = np.log10(y) + np.log10(obs_T500) - np.log10(kboltz_v2)  # put in K
    #             ly = np.log10(y) #- np.log10(kboltz_v2)  # put in K, and keep in T/T500 (instead of just T [K])
            if(property=="temperature_norm"):
                lR,y = load_Sun2009(obs_Sun_ID,"T_3d")
    #             ly = np.log10(y) + np.log10(obs_T500)
    #             ly = np.log10(y) + np.log10(obs_T500) - np.log10(kboltz_v2)  # put in K
                ly = np.log10(y) #- np.log10(kboltz_v2)  # put in K, and keep in T/T500 (instead of just T [K])
            if((property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
                lR,ly = load_Sun2009(obs_Sun_ID,"den")
                lR = lR + np.log10(0.471) # In R2500, got to get to R500
                R = 10**lR*obs_R500
                if(property=="m100kpc"):
                    m100kpc = IntegrateDensity(R,ly,100)/gperMsol
                    ly = np.log10(m100kpc)
                if(property=="m200kpc"):
                    m200kpc = IntegrateDensity(R,ly,200)/gperMsol
                    ly = np.log10(m200kpc)
                if(property=="m400kpc"):
                    m400kpc = IntegrateDensity(R,ly,400)/gperMsol
                    ly = np.log10(m400kpc)
                    
                    
            if not plot_x_log:
                lR = 10**np.array(lR)
            if not plot_y_log:
                ly = 10**np.array(ly)


            obs_lR_coll.extend(lR)
            obs_ly_coll.extend(ly)
            obs_lM500_coll.extend([obs_lM500])
            
            if use_cbar and obs_plot_kwargs['Oppenheimer+21']['use_cbar']:
                if log_cbar:
#                     color_ = cm((obs_lM500-np.log10(cbar_lo))/(np.log10(cbar_hi)-np.log10(cbar_lo)))
                    color_ = cm((obs_lM500-cbar_lo)/(cbar_hi-cbar_lo))
                else:
                    color_ = cm((10**obs_lM500-cbar_lo)/(cbar_hi-cbar_lo))
            else:
                color_ = obs_plot_kwargs['Oppenheimer+21']['sun_lovisari_color']

    #         print('bye')
            ax_coll.plot(lR,ly,color=color_,lw=2,ls="-",
                         alpha=obs_plot_kwargs['Oppenheimer+21']['sun_lovisari_alpha'],
                         zorder=obs_plot_kwargs['Oppenheimer+21']['zorder'])#,label='Lovisari and Sun Groups')
    #         ax_coll.legend()










def loadAndersonWang():
    fname = home_name + "data/anderson2015_wang2016.M500_LX.dat"
    if(os.path.isfile(fname)):
        LX,LX_err,M500,M500_err = np.loadtxt(fname,usecols=(6,7,8,9),unpack=True)

    return(M500,M500_err,LX,LX_err)
    
def loadSun2009(group_ID, property_base):
    fname = home_name + "Sun2009/" + property_base + "/" + group_ID + "_" + property_base
    #print("fname= ", fname)
    if(os.path.isfile(fname)):
        R,y = np.loadtxt(fname,usecols=(0,1),unpack=True)
        if property_base == 'temperature':
            y /= kboltz_v2
#         print("y[-1]= ", y[-1])
    else:
#         R = [-99,-99]
#         y = [-99,-99]
#         R = [np.nan, np.nan]
#         y = [np.nan, np.nan]
        R = []
        y = []

    return(R,y)

def loadCLoGS_entropy(group_ID):
    r_in,r_out,K, Kmin, Kmax = np.loadtxt(home_name + "CLoGS_entropy_profiles/LGG" + group_ID + "_entropy.dat", usecols=(0,1,2,3,4), unpack=True)
    r_mid = (r_in+r_out)/2.
    return(r_mid,K)


def loadACCEPTdata(group_ID,property):

    group_ID_prof = np.loadtxt(home_name + "ACCEPT/ACCEPT.profiles.dat", usecols=(0), dtype=(str), unpack=True)
    r_in,r_out,ne,ne_err,K_itpl,K_flat,K_err  = np.loadtxt(home_name + "ACCEPT/ACCEPT.profiles.dat", usecols=(1,2,3,4,5,6,7),unpack=True)

    #np.genfromtxt(home_name + "ACCEPT/ACCEPT.profiles.dat",dtype=(str,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float))#,unpack=True)
    
    r_mid = (r_in+r_out)/2.
    #print("group_ID,group_ID_prof=",group_ID,group_ID_prof[0],group_ID_prof[1],group_ID_prof[2],r_in[0],r_in[1],r_ins[2])
    indexes = np.char.equal(group_ID_prof,group_ID)
    #print(indexes)

    r_mid = r_mid[indexes]
    K_flat = K_flat[indexes]
    #print("rmid, Kitpl= ", r_mid, K_itpl)
    return(r_mid, K_flat)

  #  if(property=="K"):
    #    return(r_mid[indexes],K_flat[indexes])        
    

def loadSun2011pressure():
    fname = home_name + "Sun2011/Sun2011_normalizedpressure.groups.dat"
    R, P_med, P_lo, P_hi = np.loadtxt(fname,usecols=(0,1,2,3),unpack=True)
    return(R,P_med,P_lo,P_hi)

def loadLovisari2015dens():
    fname = home_name + "Lovisari2015/ne_groups.txt"
    R, ne_med, ne_lo, ne_hi = np.loadtxt(fname,usecols=(0,1,2,3),unpack=True)
    return(R,ne_med,ne_lo,ne_hi)

def IntegrateGas(Radius,Mass,OuterRadius):

    TotMass = 0.
    for i in range(len(Radius)):
        if(i>0):
            if(Radius[i]<OuterRadius):
                TotMass += Mass[i]
            if((Radius[i]>OuterRadius) & (Radius[i-1]<OuterRadius)):
                TotMass += Mass[i]*(OuterRadius-Radius[i-1])/(Radius[i]-Radius[i-1])

    return(TotMass)

def IntegrateDensity(Radius,Density,OuterRadius):

    TotMass = 0.
    Router = -99
    for i in range(len(Radius)):
        if((i>0) & (i<len(Radius)-1)):
            Router = (Radius[i] + (Radius[i+1]-Radius[i])/2.)
            Rinner = (Radius[i] - (Radius[i]-Radius[i-1])/2.)
            Volume = 4*np.pi/3*((Router)**3 - (Rinner)**3)*(cmperkpc)**3 # cm^3            
            if(Router < OuterRadius):
                TotMass += 10**Density[i]*(mu_e*mp)*Volume
            if((Router > OuterRadius) & (Rinner < OuterRadius)):
                TotMass += 10**Density[i]*(mu_e*mp)*Volume*(OuterRadius-Rinner)/(Router-Rinner)
            #print("Intdens= ", i,Router,Rinner,OuterRadius,Density[i],Volume,TotMass)
            
    if(Router>OuterRadius):    
        return(TotMass)
    else:
        return(-99)
            
def AbelTransform(Radius,Mass):

    SurfaceDensity = np.zeros(len(Radius))
        
    for i in range(len(Radius)):    
        for j in range(len(Radius)):
            if((Radius[j]>=Radius[i]) & (Radius[j]<Radius[-1])):
                path = (Radius[j+1]-Radius[j])*Radius[j+1]/np.sqrt(Radius[j+1]**2-Radius[i]**2)
                delta_R = Radius[j+1]-Radius[j]
                mean_R = (Radius[j+1]+Radius[j])/2.
                Volume = 4*np.pi/3*((mean_R+delta_R/2.)**3-(mean_R-delta_R/2.)**3)
                #print(i,j,Volume, Mass[j], path)
                SurfaceDensity[i] += 2*path*Mass[j]/Volume
    
    return(SurfaceDensity)

def flexbin_med_and_onesigmaspread(x,y,nbins):

    xmed = np.zeros(nbins)
    xlow = np.zeros(nbins)
    xhi  = np.zeros(nbins)

    ymed = np.zeros(nbins)
    ylow = np.zeros(nbins)
    yhi  = np.zeros(nbins)
    
    for i in range(nbins):
        xmed[i] = np.percentile(x,(i+0.5)/nbins*100.)
        xlow[i] = np.percentile(x,(i+0.0)/nbins*100.)
        xhi[i] = np.percentile(x,(i+1.0)/nbins*100.)
        indexes_x = np.where((x>=xlow[i]) & (x<=xhi[i]))
        #print("xlo,xhi,index= ",xlow[i],xhi[i],indexes_x)
        ymed[i] = np.percentile(y[indexes_x],50)
        ylow[i] = np.percentile(y[indexes_x],16)
        yhi[i] = np.percentile(y[indexes_x],84)

    return(xmed,xlow,xhi,ymed,ylow,yhi)

def return_most_massive_topbin(x,y,nbins):

    xmin = np.percentile(x,(nbins-0.5)/nbins*100.)
    massive_indexes = np.where(x>=xmin)
#     print("xmin,massive_index= ",xmin,massive_indexes)
    return(x[massive_indexes],y[massive_indexes])


def calc_med_and_flexspread(x, y, nbins,bin_lo,bin_hi,percent_lo,percent_hi):
    med = scist.binned_statistic(np.asarray(x), np.asarray(y),statistic=lambda y: np.percentile(y, 50), bins=nbins, range=[(bin_lo,bin_hi)])[0]
    range_low = scist.binned_statistic(np.asarray(x), np.asarray(y),statistic=lambda y: np.percentile(y, percent_lo), bins=nbins, range=[(bin_lo,bin_hi)])[0]
    range_high = scist.binned_statistic(np.asarray(x), np.asarray(y),statistic=lambda y: np.percentile(y, percent_hi), bins=nbins, range=[(bin_lo,bin_hi)])[0]
    mean = scist.binned_statistic(np.asarray(x), np.asarray(y),statistic=np.mean, bins=nbins, range=[(bin_lo,bin_hi)])[0]
    return(med,range_low,range_high,mean)





# CALL run_plotting_script(...) within THIS SCRIPT TO PLOT
# scaling_relation_functions is a script that aviv has



def run_plotting_script(property, lhalo_lo, lhalo_hi, axes_=plt, sm=None, cm=plt.get_cmap('plasma'), plot_individual_observations=0, plot_only_Lovisari_Sun_groups=True, plot_x_log=True, plot_y_log=True, use_cbar=False, log_cbar=False, obs_plot_kwargs=None):

    # print(sys.argv)
    # file_list = sys.argv[1]
#     property = sys.argv[1]
#     coloring = "mass"
    # legend = int(sys.argv[3])
    # multi = int(sys.argv[2])
    
#     if use_cbar:
#     if sm==None:
#         cbar_lo = 10**lhalo_lo
#         cbar_hi = 10**lhalo_hi
#     else:
#         cbar_lo = sm.get_clim()[0]
#         cbar_hi = sm.get_clim()[1]
        
    try:
        cbar_lo = sm.get_clim()[0]
        cbar_hi = sm.get_clim()[1]
        if obs_plot_kwargs['Oppenheimer+21']['verbose']:
            print('Oppenheimer+21 cbar_lo:', cbar_lo)
            print('Oppenheimer+21 cbar_hi:', cbar_hi)
    except:
        cbar_lo = 0
        cbar_hi = 0
        if obs_plot_kwargs['Oppenheimer+21']['verbose']:
            print('No sm for colorbar, likely because not using colorbar')
        
    
    

    # base_name = "/data1/extboppe/"
    # base_name = "/hpcdata5/users/extboppe/"

    # title_label = file_list.split(".")[0]

#     cm = plt.get_cmap('plasma') 
    # if(multi==1):
    #     cm = plt.get_cmap('viridis') 

# -----------------------------------------------------------------------------------------------------------for r^1.1 & r^0.7
    lRlo = -2.0
    # lRhi = 0.6 # 0.3
    lRhi = np.log10(3.0)
    nrbins = 26
    lRbins = np.linspace(lRlo,lRhi,nrbins)

    #if(len(sys.argv)>5):
#     lhalo_lo = float(sys.argv[2])
#     lhalo_hi = float(sys.argv[3])
    # nhalo_bins = int(sys.argv[5])
    # mhalo_str = '.MH%d_%d_%dbins'%(lhalo_lo*10,lhalo_hi*10,nhalo_bins) 


    # fig, axes_ = plt.subplots(nrows=1, ncols=1)
    # axes_ = eval(sys.argv[5])


    #else:
    #    nhalo_bins = 1
    #    mhalo_str = ''


#     plot_individual_observations = 0
    # if(multi==1):
    #     plot_individual_observations = 0
    # else:
    #     plot_individual_observations = 0

#     topborder = 0.93
#     leftborder = 0.16
#     rightborder = 0.98
    # if(multi==2):
    #     topborder = 0.98

#     if(coloring=="mass"):
#         cbarlabel = "log$M_{500}/M_{\odot}$"

    # if(property=="entropy_norm"):
    #     lylo = -1.5
    #     lyhi = 1.0
    #     if((multi==0) | (multi==2)):
    #         lylo = -2.0
    #         lyhi = 0.5
    #         if(multi==2):
    #             lylo = -2.2

    #     ylabel = "log $K/K_{500}$" # [keV cm$^{2}$]"
    #     leftborder = 0.14
    # if(property=="pressure_norm"):
    #     lylo = -1.0
    #     lyhi = 2.0
    #     ylabel = "log $P/P_{500}$" # [keV cm$^{2}$]"
    #     leftborder = 0.14
    # if(property=="pressure"):
    #     lylo = 1.0
    #     lyhi = 6.0
    #     ylabel = "log $P$" # [keV cm$^{2}$]"
    #     leftborder = 0.14

    # if(property=="density"):
    #     lylo = -5.0 #-4.0
    #     lyhi = -1.0
    #     ylabel = "log $n_e$ [cm$^{-3}$]"
    # if(property=="temperature"):
    #     lylo = -1.0
    #     lyhi = 1.0
    #     ylabel = "log $T$ [keV]"
    # if(property=="metallicity"):
    #     lylo = -1.0
    #     lyhi = 0.0
    #     ylabel = "log $Z/Z_{\odot}$"
    #     Zsol = 0.0134
    # if(property=="XraySB_XMMcounts"):
    #     lylo = -4.0
    #     lyhi = 0.0
    #     ylabel = "log SB [cts s$^{-1}$ arcmin$^{-2}$]"
    #     XMM_erg_to_count_conversion = -39.5
    #     lRlo = -2.0
    #     lRhi = 0.6 # 0.3
    #     nrbins = 9
    #     lRbins = np.linspace(lRlo,lRhi,nrbins)
    # if(property=="coolgas"):
    #     lylo = 3.0
    #     lyhi = 8.0
    #     ylabel = "log $\Sigma_{\mathrm{cool}}$ [$M_{\odot}$ kpc$^{-2}$]"
    # if(property=="warmhotgas"):
    #     lylo = 4.0
    #     lyhi = 9.0
    #     ylabel = "log $\Sigma_{\mathrm{warm-hot}}$ [$M_{\odot}$ kpc$^{-2}$]"
    # if(property=="dispersionmeasure"):
    #     ylo = 0
    #     yhi = 1000
    #     ylabel = "log DM$_{\mathrm{CGM}}$ [pc cm$^{-3}$]"
    #     lRlo = 1.0
    #     lRhi = 3.0 # 0.3
    #     nrbins = 25
    #     lRbins = np.linspace(lRlo,lRhi,nrbins)
    # if(property=="fgas"):
    #     lylo = 0.00
    #     lyhi = 0.16
    #     lMlo = 11.9
    #     lMhi = 14.7
    #     cbar_lo = lMlo
    #     cbar_hi = lMhi
    #     ylabel = "$f_{\mathrm{gas}}$"
    # if((property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
    #     lylo = 10.5
    #     lyhi = 13.0
    #     lMlo = 12.9
    #     lMhi = 14.6
    #     cbar_lo = lMlo
    #     cbar_hi = lMhi
    #     if(property=="m100kpc"):
    #         ylabel = "log $M_{\mathrm{IGrM,<100 kpc}}$ [$M_{\odot}$]" 
    #     if(property=="m200kpc"):
    #         ylabel = "log $M_{\mathrm{IGrM,<200 kpc}}$ [$M_{\odot}$]" 
    #     if(property=="m400kpc"):
    #         ylabel = "log $M_{\mathrm{IGrM,<400 kpc}}$ [$M_{\odot}$]" 
    # if(property=="LX"):
    #     lylo = 37.0
    #     lyhi = 44.5
    #     lMlo = 12.0
    #     lMhi = 15.0
    #     cbar_lo = lMlo
    #     cbar_hi = lMhi
    #     ylabel = "log $L_{\mathrm{X}}$ [erg s$^{-1}$]" 
    # if((property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
    #     fig_all = plt.figure(figsize=(4.5,5.0))
    #     ax_all = fig_all.add_subplot(111)

    #     fig_coll = plt.figure(figsize=(4.5,5.0))
    #     ax_coll = fig_coll.add_subplot(111)
    #     leftborder = 0.18
    # else:
    #     fig_all = plt.figure(figsize=(6.0,5.0))
    #     ax_all = fig_all.add_subplot(111)

    #     fig_coll = plt.figure(figsize=(6.0,5.0))
    #     ax_coll = fig_coll.add_subplot(111)

    # if(property=="LX"):
    #     fig_all = plt.figure(figsize=(5.0,5.0))
    #     ax_all = fig_all.add_subplot(111)

    #     fig_coll = plt.figure(figsize=(5.0,5.0))
    #     ax_coll = fig_coll.add_subplot(111)
    #     leftborder = 0.18


    f_b = 0.04825/0.307


    # if(multi):
    #     sim_list = file_list
    #     sims = open(sim_list,"r").readlines()
    #     nsims = len(sims)
    # else:
    #     nsims = 1

    # file_base = file_list

    # for isim in range(nsims):

    #     lM500_coll = []
    #     mass_coll = []
    #     R_coll = [[1e-99 for i in range(1)] for i in range(nhalo_bins)]
    #     y_coll = [[1e-99 for i in range(1)] for i in range(nhalo_bins)]
    #     y_median = [[1e-99 for i in range(1)] for i in range(nhalo_bins)]
    #     y_low = [[1e-99 for i in range(1)] for i in range(nhalo_bins)]
    #     y_high = [[1e-99 for i in range(1)] for i in range(nhalo_bins)]
    #     y_mean = [[1e-99 for i in range(1)] for i in range(nhalo_bins)]
    #     count_halos = np.zeros(nhalo_bins)

    #     if(multi):
    #         file_list = sims[isim].split()[0]
    #         sim_color = sims[isim].split()[1]
    #         sim_label = sims[isim].split()[2]
    #         sim_linestyle = sims[isim].split()[3]
    #         sim_zorder = int(sims[isim].split()[4])

    #         cbar_lo = lhalo_lo
    #         cbar_hi = lhalo_hi
    #         if(property=="entropy_norm"):
    #             title_label = "Entropy"
    #         if((property=="pressure_norm") | (property=="pressure")):
    #             title_label = "Pressure"
    #         if(property=="density"):
    #             title_label = "Density"
    #         if(property=="temperature"):
    #             title_label = "Temperature"            
    #         if(property=="metallicity"):
    #             title_label = "Metallicity"            
    #         if(property=="XraySB_XMMcounts"):
    #             title_label = "Xray Surface Brightness"            
    #         if(property=="fgas"):
    #             title_label = "Gas Fraction inside $R_{500}$"
    #         if(property=="m100kpc"):
    #            title_label = "Hot Gas Inside 100 kpc"
    #         if(property=="m200kpc"):
    #            title_label = "Hot Gas Inside 200 kpc"
    #         if(property=="m400kpc"):
    #            title_label = "Hot Gas Inside 400 kpc"
    #         if(property=="LX"):
    #            title_label = "X-ray Luminosity 0.15-1.0$R_{500}$"
    #         if(property=="dispersionmeasure"):
    #             title_label = "Dispersion Measure"            


    #     else:
    #         sim_color = "blue"
    #         sim_label = title_label
    #         sim_linestyle = "-"
    #         sim_zorder = 5

    #     halos = open(file_list,"r").readlines()

    #     nhalos = len(halos)

    #     for i in range(nhalos):

    #         file_type = halos[i].split()[0].split(".")[-1]
    #         if(file_type == "hdf5"):
    #             halo_hdf5_file = h5py.File(base_name + halos[i].split()[0])

    #             print(halo_hdf5_file)

    #             if(property=="LX"):
    #                 #if not "/EventFileAnalysis".attrs["LSoft_01R200"] in halo_hdf5_file:
    #                 if not "/EventFileAnalysis/fR200_LSoftDens" in halo_hdf5_file:
    #                     continue
    #             else:
    #                 if not "/RadialProfiles/R_Mass" in halo_hdf5_file:
    #                     continue

    #             redshift = 0.0
    #             lM200 = halo_hdf5_file['Halo'].attrs['lM200']
    #             lMStar = halo_hdf5_file['Halo'].attrs['lMStar']
    #             SFR = halo_hdf5_file['Halo'].attrs['SFR']
    #             if(lM200>20): # SIMBA got is wrong
    #                 lM200= np.log10(lM200)
    #             if(lMStar>20): # SIMBA is wrong
    #                 lMStar= np.log10(lMStar)
    #             M500 = 10**(lM200-0.16) # Msol
    #             lM500 = np.log10(M500)
    #         else:
    #             home_name =  "/home/ariboppe/groups_review/"
    #             ascii_file = home_name + halos[i].split()[0]
    #             #print(halos[i].split("_")[2][:-5])
    #             if(property=="entropy_norm"):
    #                 print("halos[i]= ", halos[i])
    #                 if(file_type == "hdf5"):
    #                     redshift = float(halos[i].split("_")[2][:-5])
    #                 else:
    #                     redshift = float(halos[i].split("_")[1][:-5])

    #                 M500,entropy,R = np.genfromtxt(ascii_file, usecols=(0,8,9),unpack=True)
    #                 M500 = M500[2]
    #                 entropy = entropy[2:]
    #                 R = R[2:]
    #                 lM500 = np.log10(M500)
    #             if((property=="density") | (property=="temperature") | (property=="pressure_norm") | (property=="pressure")):
    #                 print(halos[i])
    #                 redshift = float(halos[i].split("_")[1][:-5])
    #                 M500,TkeV,density,R = np.genfromtxt(ascii_file, usecols=(0,8,9,10),unpack=True)
    #                 M500 = M500[2]
    #                 temperature = TkeV[2:]/keVperK
    #                 density = density[2:]
    #                 pressure = density*temperature
    #                 R = R[2:]
    #                 lM500 = np.log10(M500)

    #                 lM200 = lM500+0.16
    #                 lMStar = -99
    #                 #print("M500= ", M500)
    #                 #print("redshift= ", redshift)


    #         #print("hdf5name, SFR,lMStar,LM200=",halo_hdf5_file, SFR,lMStar,lM200)

    #         if((property=="entropy_norm") | (property=="density") | (property=="temperature") | (property=="pressure_norm") | (property=="pressure") | (property=="metallicity") | (property=="XraySB_XMMcounts") |  (property=="coolgas") |  (property=="warmhotgas") | (property=="dispersionmeasure")):
    #             if( (lM500>=lhalo_lo) & (lM500<lhalo_hi)):
    #                 hbin =int((lM500-lhalo_lo)/((lhalo_hi-lhalo_lo)/nhalo_bins))
    #             else:
    #                 continue 

    #         rhocrit = 8.52e-30*(1+redshift)**3 # g/cm^3  
    #         ne500 = 500*f_b*rhocrit/(mu_e*mp) # cm^-3

    #         rho500 = 500*rhocrit  # g/cm^3
    #         R500 = (M500*1.989e+33/rho500*(3/(4*np.pi)))**0.3333   # cm

    #         T500 = (G*M500*1.989e+33*mu*mp/(2*R500*kboltz))
    #         K500 = T500*ne500**-0.6666*keVperK

    #         P500 = ne500*T500

    #         print ("ID,M500,R500,T500,ne500,K500= ",halos[i].split()[0],np.log10(M500),np.log10(R500),np.log10(T500),np.log10(ne500),np.log10(K500))

    #         if(property=="entropy_norm"):
    #             if(file_type == "hdf5"):
    #                 R = np.asarray(halo_hdf5_file["RadialProfiles/R"]) # kpc
    #                 entropy = np.asarray(halo_hdf5_file["RadialProfiles/Entropy_hot"])
    #                 entropy = np.where(np.isnan(entropy), 1e-99, entropy)
    #             #print(R, entropy)
    #             R = R*cmperkpc/R500
    #             y = entropy/K500
    #             #print(R, y)
    #         if((property=="pressure_norm") | (property=="pressure")):
    #             if(file_type == "hdf5"):
    #                 R = np.asarray(halo_hdf5_file["RadialProfiles/R"]) # kpc
    #                 density = np.asarray(halo_hdf5_file["RadialProfiles/Density_hot"])
    #                 temperature = np.asarray(halo_hdf5_file["RadialProfiles/Temperature_hot"])
    #                 pressure = density*temperature/(mp*mu_e)
    #                 pressure = np.where(np.isnan(pressure), 1e-99, pressure)
    #                 #pressure = np.asarray(halo_hdf5_file["RadialProfiles/Pressure_hot"])
    #             R = R*cmperkpc/R500
    #             if(property=="pressure_norm"):
    #                 y = pressure/P500
    #             if(property=="pressure"):
    #                 y = pressure
    #             #print("pressure = " , y)
    #             #print(R, y)
    #         if(property=="density"):
    #             if(file_type == "hdf5"):
    #                 R = np.asarray(halo_hdf5_file["RadialProfiles/R"]) # kpc
    #                 density = np.asarray(halo_hdf5_file["RadialProfiles/Density_hot"])
    #                 density = np.where(np.isnan(density), 1e-99, density)
    #             R = R*cmperkpc/R500
    #             if(file_type == "hdf5"):
    #                 y = density/(mp*mu_e)
    #             else:
    #                 y = density
    #         if(property=="temperature"):
    #             if(file_type == "hdf5"):
    #                 R = np.asarray(halo_hdf5_file["RadialProfiles/R"]) # kpc
    #                 temperature = np.asarray(halo_hdf5_file["RadialProfiles/Temperature_hot"])
    #                 temperature = np.where(np.isnan(temperature), 1e-99, temperature)
    #             R = R*cmperkpc/R500
    #             y = temperature*keVperK
    #         if(property=="metallicity"):
    #             if(file_type == "hdf5"):
    #                 R = np.asarray(halo_hdf5_file["RadialProfiles/R"]) # kpc
    #                 metallicity = np.asarray(halo_hdf5_file["RadialProfiles/Metallicity_hot"])
    #                 metallicity = np.where(np.isnan(metallicity), 1e-99, metallicity)
    #             R = R*cmperkpc/R500
    #             y = metallicity/Zsol
    #         if(property=="XraySB_XMMcounts"):
    #             XraySB = np.asarray(halo_hdf5_file["EventFileAnalysis/LSoftDens"])
    #             y = 10**(np.log10(XraySB) + XMM_erg_to_count_conversion)
    #             R_edges = np.asarray(halo_hdf5_file["EventFileAnalysis/R"])
    #             R_kpc = np.asarray([10**((np.log10(R_edges[i+1])+np.log10(R_edges[i]))/2.) for i in range(len(R_edges)-1)])
    #             print("R_kpc,y= ", R_kpc,y)
    #             R = R_kpc*cmperkpc/R500            
    #         if(property=="coolgas"):
    #             if(file_type == "hdf5"):            
    #                 R = np.asarray(halo_hdf5_file["RadialProfiles/R_Mass"]) # kpc
    #                 coolgas = np.asarray(halo_hdf5_file["RadialProfiles/Mass_cool"]) # g
    #             SurfaceDensity = AbelTransform(R,coolgas)/gperMsol
    #             R = R*cmperkpc/R500
    #             y = SurfaceDensity
    #         if((property=="warmhotgas") | (property=="dispersionmeasure")):
    #             if(file_type == "hdf5"):
    #                 R = np.asarray(halo_hdf5_file["RadialProfiles/R_Mass"]) # kpc
    #                 warmhotgas = np.asarray(halo_hdf5_file["RadialProfiles/Mass_warmhot"]) # g
    #             SurfaceDensity = AbelTransform(R,warmhotgas)
    #             if(property=="warmhotgas"):
    #                 R = R*cmperkpc/R500
    #                 y = SurfaceDensity/gperMsol # Msol/kpc^2
    #             if(property=="dispersionmeasure"):
    #                 R = R #kpc
    #                 y = SurfaceDensity/(mp*mu_e) #electrons/kpc^2
    #                 y = y/cmperkpc**2/cmperpc
    #         if(property=="fgas"):
    #             R = np.asarray(halo_hdf5_file["RadialProfiles/R"]) # kpc
    #             if(file_type == "hdf5"):
    #                 R = np.asarray(halo_hdf5_file["RadialProfiles/R_Mass"]) # kpc
    #                 allgas = np.asarray(halo_hdf5_file["RadialProfiles/Mass"]) # g

    #             Mgas500 = IntegrateGas(R,allgas,R500/cmperkpc)/gperMsol
    #             fgas500 = Mgas500/M500
    #             y = fgas500
    #         if((property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
    #             R = np.asarray(halo_hdf5_file["RadialProfiles/R"]) # kpc
    #             if(file_type == "hdf5"):
    #                 R = np.asarray(halo_hdf5_file["RadialProfiles/R_Mass"]) # kpc
    #                 hotgas = np.asarray(halo_hdf5_file["RadialProfiles/Mass_hot"]) # g

    #             if(property=="m100kpc"):
    #                 y = IntegrateGas(R,hotgas,100)/gperMsol
    #             if(property=="m200kpc"):
    #                 y = IntegrateGas(R,hotgas,200)/gperMsol
    #             if(property=="m400kpc"):
    #                 y = IntegrateGas(R,hotgas,400)/gperMsol

    #             y = np.log10(y)

    #         if(property=="LX"):
    #             y = np.log10(halo_hdf5_file['EventFileAnalysis'].attrs['LSoft_01R200'])
    #             ax_all.plot(lM500,y,marker='o',color=cm((lM500-cbar_lo)/(cbar_hi-cbar_lo)),zorder=int(lM500*100))

    #         if((property=="fgas") | (property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc") | (property=="LX")):
    #             ax_all.plot(lM500,y,marker='o',color=cm((lM500-cbar_lo)/(cbar_hi-cbar_lo)),zorder=int(lM500*100))
    #             print("LM500,y= ",lM500,y)
    #             lM500_coll.extend([lM500])
    #             mass_coll.extend([y])

    #             #print("M,fgas=", lM500, fgas500, lM_coll,fgas_coll)

    #         else:     # radial profile               
    #             if((lM500>=cbar_lo) & (lM500<cbar_hi)):
    #                 if(property=="dispersionmeasure"):
    #                     ax_all.plot(np.log10(R),y,color=cm((lM500-cbar_lo)/(cbar_hi-cbar_lo)),zorder=int(lM500*100))
    #                 else:
    #                     print("R= ", R, len(R))
    #                     print("y= ", y, len(y))
    #                     ax_all.plot(np.log10(R),np.log10(y),color=cm((lM500-cbar_lo)/(cbar_hi-cbar_lo)),zorder=int(lM500*100))

    #             if((hbin>=0) & (hbin<nhalo_bins)):
    #                 R_coll[hbin].extend(R)
    #                 y_coll[hbin].extend(y)
    #                 count_halos[hbin] += 1

    #     if((property=="fgas") | (property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc") | (property=="LX")):
    #         nflexbins = int(len(mass_coll)/20)
    #         if(nflexbins>30):
    #             nflexbins = 30
    #         M,Mlo,Mhi,y,ylo,yhi=flexbin_med_and_onesigmaspread(np.asarray(lM500_coll),np.asarray(mass_coll),nflexbins)
    #         M_massive,y_massive = return_most_massive_topbin(np.asarray(lM500_coll),np.asarray(mass_coll),nflexbins)

    #         ax_coll.plot(M, y, color = sim_color ,ls=sim_linestyle,lw=4,label=sim_label,zorder=sim_zorder) #max(y)*100)
    #         ax_coll.fill_between(M,ylo,yhi, facecolor=sim_color, alpha=0.2,zorder=sim_zorder) #max(y)*100-1)
    #         ax_coll.plot(M_massive,y_massive,color=sim_color,marker='.',ls='',zorder=sim_zorder) #max(y)*100)
    #     else:
    #         for h in range(nhalo_bins):

    #             lMhalo  = lhalo_lo + (lhalo_hi-lhalo_lo)*(h+0.5)/nhalo_bins
    #             lMhalo_min = lhalo_lo + (lhalo_hi-lhalo_lo)*(h)/nhalo_bins
    #             lMhalo_max = lhalo_lo + (lhalo_hi-lhalo_lo)*(h+1)/nhalo_bins

    #             if(property=="dispersionmeasure"):
    #                 y_median[h], y_low[h], y_high[h], y_mean[h] = calc_med_and_flexspread(np.log10(R_coll[h]),y_coll[h],nrbins,lRlo,lRhi,16,84)
    #                 print("y_median= ",y_median[h])
    #             else:
    #                 y_median[h], y_low[h], y_high[h], y_mean[h] = calc_med_and_flexspread(np.log10(R_coll[h]),y_coll[h],nrbins,lRlo,lRhi,16,84)
    #                 print("y_median= ",y_median[h])

    #             #print("calc_med_and_flexspread= ", y_median[h], y_low[h], y_high[h],len(y_coll[h]), sim_label, nrbins, count_halos[h])
    #             if(multi):
    #                 if(multi==1):
    #                     if(property=="dispersionmeasure"):                    
    #                         ax_coll.plot(lRbins, y_median[h], color = sim_color,lw=4,label=sim_label)
    #                         ax_coll.fill_between(lRbins,y_low[h],y_high[h], facecolor=sim_color, alpha=0.2)
    #                         ###ax_coll.plot(lRbins, y_mean[h], color = sim_color,lw=4, ls=':')                    
    #                     else:
    #                         ax_coll.plot(lRbins, np.log10(y_median[h]), color = sim_color,lw=4,label=sim_label)
    #                         ax_coll.fill_between(lRbins,np.log10(y_low[h]),np.log10(y_high[h]), facecolor=sim_color, alpha=0.2)
    #                         ###ax_coll.plot(lRbins, np.log10(y_mean[h]), color = sim_color,lw=4, ls=':')                    
    #                 else: #multi==2
    #                     if(sim_linestyle=="-"):
    #                         ax_coll.plot(lRbins, np.log10(y_median[h]), color = cm((lMhalo-cbar_lo)/(cbar_hi-cbar_lo)),zorder=int(lMhalo*100),ls=sim_linestyle, lw=4) ###,label="log[$M_{500}$]=%4.1f-%4.1f"%(lMhalo_min,lMhalo_max)) #BDO 5/14/21, replacing with colorbar
    #                         ax_coll.fill_between(lRbins,np.log10(y_low[h]),np.log10(y_high[h]), facecolor=cm((lMhalo-cbar_lo)/(cbar_hi-cbar_lo)), alpha=0.2)
    #                         ###ax_coll.plot(lRbins, np.log10(y_mean[h]), color = cm((lMhalo-cbar_lo)/(cbar_hi-cbar_lo)),zorder=int(lMhalo*100),ls=":", lw=4)
    #                     else:
    #                         ax_coll.plot(lRbins, np.log10(y_median[h]), color = cm((lMhalo-cbar_lo)/(cbar_hi-cbar_lo)),zorder=int(lMhalo*100),ls=sim_linestyle, lw=4)
    #                         # No fill_between for 2nd simulation.  

    #             else:
    #                 if(property=="dispersionmeasure"):                    
    #                     ax_coll.plot(lRbins, y_median[h], color = cm((lMhalo-cbar_lo)/(cbar_hi-cbar_lo)),zorder=int(lMhalo*100),lw=4)
    #                     ax_coll.fill_between(lRbins,y_low[h],y_high[h], facecolor=cm((lMhalo-cbar_lo)/(cbar_hi-cbar_lo)), alpha=0.2)
    #                     ###ax_coll.plot(lRbins, y_mean[h], color = cm((lMhalo-cbar_lo)/(cbar_hi-cbar_lo)),zorder=int(lMhalo*100),ls=":",lw=4)
    #                 else:
    #                     ax_coll.plot(lRbins, np.log10(y_median[h]), color = cm((lMhalo-cbar_lo)/(cbar_hi-cbar_lo)),zorder=int(lMhalo*100),lw=4) ###,label="log[$M_{500}$]=%4.1f-%4.1f"%(lMhalo_min,lMhalo_max)) #BDO 5/14/21, replacing with colorbar
    #                     ax_coll.fill_between(lRbins,np.log10(y_low[h]),np.log10(y_high[h]), facecolor=cm((lMhalo-cbar_lo)/(cbar_hi-cbar_lo)), alpha=0.2)
    #                     ###ax_coll.plot(lRbins, np.log10(y_mean[h]), color = cm((lMhalo-cbar_lo)/(cbar_hi-cbar_lo)),zorder=int(lMhalo*100),ls=":",lw=4)


    #             ###

    #             data_output_filename_rvir = 'profile_coll_' + sims[isim].split()[2] + '.' + property  + mhalo_str +  '.dat'

    #             data_output_file_rvir = open(data_output_filename_rvir, "w")

    #             data_output_file_rvir.write("#R/R500 median 16percentile 84percentile nhalos  \n")

    #             for k in range(len(lRbins)):
    #                 data_output_file_rvir.write("%5.3e %5.2f %5.2f %5.2f %4d\n"%(lRbins[k],np.log10(y_median[h][k]),np.log10(y_low[h][k]),np.log10(y_high[h][k]),count_halos[h]))

    #             data_output_file_rvir.close()
    #             ###

    # if(multi==2):
    #     for isim in range(nsims):
    #         sim_label = sims[isim].split()[2]
    #         sim_linestyle = sims[isim].split()[3]
    #         ax_coll.plot([-99,-99],[-99,-99], color = 'k', ls=sim_linestyle, lw=4, label=sim_label)


    if(property=="entropy_norm"):
         # r^1.1
         axes_.plot(lRbins, np.log10(1.32/1.09 * (10**(lRbins))**1.1),color="gray",lw=2,ls=":", label = r'K ∝ $R^{1.1}$ (Voit+2005)') # from (Lewis+2000, Babul+2002, Voit+2005)
         #axes_.plot(lRbins, np.log10(1.32/1.09 * (10**(lRbins))**1.1),color="gray",lw=2,ls=":",label="Voit+2005")

         # r^0.7
         axes_.plot(lRbins, np.log10(95.4 * (10**lRbins / 100)**0.7),color="black",lw=2,ls="--",label=r'K ∝ $R^{0.7}$ (Panagoulia+2014)')


    if((property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
        lR = np.linspace(-1.0,3.3,87) # log kpc
        lM500 = np.linspace(12.9,15.1,23) # log Msolar

        Mbase_all = []
        lM500_all = []
        for i in range(len(lM500)):        
            R500 = (10**lM500[i]*1.989e+33/rho500*(3/(4*np.pi)))**0.3333   # cm
            v_c = np.sqrt(G*10**lM500[i]*gperMsol/(R500))
            Density = f_b*v_c**2/(4*np.pi*G*10**lR*cmperkpc) #/(mu_e*mp)
            #Density = (1.32* (((10**lR)/(R500/cmperkpc))**1.1)*K500/T500)**-1.5 # cm**-3
#             print("R500,v_c, Density= " , R500,v_c, Density)
            if(property=="m100kpc"):
                Mbase = IntegrateDensity(10**lR,Density,100)/gperMsol
            if(property=="m200kpc"):
                Mbase = IntegrateDensity(10**lR,Density,200)/gperMsol
            if(property=="m400kpc"):
                Mbase = IntegrateDensity(10**lR,Density,400)/gperMsol

#             print("LM500,Mbase= ",lM500[i],Mbase)
            Mbase_all.extend([Mbase])
            #lM500_all.extend([lM500])

        if not plot_x_log:
            lM500 = 10**lM500
        if not plot_y_log:
            Mbase_all = 10**Mbase_all
        
#         print("lM500,Mbase= ",lM500,np.log10(Mbase_all))
        axes_.plot(lM500,np.log10(Mbase_all),color="gray",lw=2,ls=":")


        #load Sun 2009 data:

#     plot_only_Lovisari_Sun_groups = True
    if((property=="entropy_norm") | (property=="density") | (property=="temperature") | (property=="ktemperature") | (property=="temperature_norm") | (property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
        if(plot_only_Lovisari_Sun_groups):
#             print('hi')
            load_and_plot_Lovisari_Sun_groups(property, axes_, cbar_lo, cbar_hi, lhalo_lo, lhalo_hi, 
                                              cm=cm, plot_x_log=plot_x_log, plot_y_log=plot_y_log, 
                                              use_cbar=use_cbar, log_cbar=log_cbar, obs_plot_kwargs=obs_plot_kwargs)
        else:
            load_and_plot_Sun2009(property, axes_, cbar_lo, cbar_hi, lhalo_lo, lhalo_hi, cm=cm, 
                                  plot_individual_observations=plot_individual_observations, 
                                  plot_x_log=plot_x_log, plot_y_log=plot_y_log, 
                                  use_cbar=use_cbar, log_cbar=log_cbar, obs_plot_kwargs=obs_plot_kwargs)

    #Sun_lR_coll = []
    #Sun_ly_coll = []
    #Sun_lM500_coll = []
    #if((property=="entropy_norm") | (property=="density") | (property=="temperature") | (property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
    #    Sungroups_list = home_name + "Sun2009/Sun_groups.cat.sort"

    #    Sungroups = open(Sungroups_list,"r").readlines()
    #    nSungroups = len(Sungroups)

    #    Sun_label = 0 

    #    for i in range(nSungroups):
    #        Sun_ID = Sungroups[i].split()[0]
    #        Sun_T500 = float(Sungroups[i].split()[1])
    #        Sun_R500 = float(Sungroups[i].split()[2])
    #        Sun_M500_13 = float(Sungroups[i].split()[3])
    #        if(Sun_M500_13<0):
    #            Sun_M500_13 = 3e+13*(Sun_T500/1.0)**(1.5)/1e+13

    #        Sun_lM500 = np.log10(Sun_M500_13)+13
    #        if(property=="entropy_norm"):
    #            lR,ly = loadSun2009(Sun_ID,"K")
    #        if(property=="density"):
    #            lR,ly = loadSun2009(Sun_ID,"den")
    #            lR = lR + np.log10(0.471) # In R2500, got to get to R500
    #        if(property=="temperature"):
    #            lR,y = loadSun2009(Sun_ID,"T_3d")
    #            ly = np.log10(y) + np.log10(Sun_T500)
    #        if((property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
    #            lR,ly = loadSun2009(Sun_ID,"den")
    #            lR = lR + np.log10(0.471) # In R2500, got to get to R500
    #            R = 10**lR*Sun_R500
    #            if(property=="m100kpc"):
    #                m100kpc = IntegrateDensity(R,ly,100)/gperMsol
    #                ly = np.log10(m100kpc)
    #            if(property=="m200kpc"):
    #                m200kpc = IntegrateDensity(R,ly,200)/gperMsol
    #                ly = np.log10(m200kpc)
    #            if(property=="m400kpc"):
    #                m400kpc = IntegrateDensity(R,ly,400)/gperMsol
    #                ly = np.log10(m400kpc)

    #        if((property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
    #            ax_all.plot(np.log10(Sun_M500_13*1e+13),ly,color="black",marker="x",zorder=0)
    #            if(i==0):
    #                ax_coll.plot(np.log10(Sun_M500_13*1e+13),ly,color="black",marker="x",ls="",ms=12,zorder=0,label="Sun+ (2009)")
    #            else:
    #                ax_coll.plot(np.log10(Sun_M500_13*1e+13),ly,color="black",marker="x",ls="",ms=12,zorder=0)

    #        else:
    #            Sun_lR_coll.extend(lR)
    #            Sun_ly_coll.extend(ly)
    #            Sun_lM500_coll.extend([Sun_lM500])
    #            if(plot_individual_observations):
    #                 if(Sun_label == 0):
    #                    if((Sun_lM500 > lhalo_lo) & (Sun_lM500 < lhalo_hi)):
    #                        #if(multi == 2 ):  # BDO 5/13/21, trying to replace inferno colorbar.   
    #                        ax_coll.plot(lR,ly,color=cm((Sun_lM500-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls=":",zorder=0,label="Sun+ (2009)")
    #                        #else:
    #                        #    ax_coll.plot(lR,ly,color='k',lw=1,ls="--",zorder=0,label="Sun+ (2009)")
    #                        Sun_label = 1
    #                else:
    #                    if((Sun_lM500 > lhalo_lo) & (Sun_lM500 < lhalo_hi)):
    #                        #if(multi == 2):# BDO 5/13/21, trying to replace viridis colorbar.
    #                            ax_coll.plot(lR,ly,color=cm((Sun_lM500-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls=":",zorder=0)
    #                        #else:
    #                        #    ax_coll.plot(lR,ly,color='k',lw=1,ls="--",zorder=0)
    #                            
    #            else:
    #                ax_coll.plot(lR,ly,color=cm((Sun_lM500-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls=":",zorder=0)

    #    if(plot_individual_observations==0):

    #        Sun_lR,Sun_lRlo,Sun_lRhi,Sun_ly,Sun_lylo,Sun_lyhi=flexbin_med_and_onesigmaspread(np.asarray(Sun_lR_coll),np.asarray(Sun_ly_coll),15)
    #        print("Sun_lR,Sun_lRlo,Sun_lRhi,Sun_ly,Sun_lylo,Sun_lyhi,Median_M500=",Sun_lR,Sun_lRlo,Sun_lRhi,Sun_ly,Sun_lylo,Sun_lyhi,np.median(np.asarray(Sun_lM500_coll)))
    #        Sun_fractional_color = (np.median(np.asarray(Sun_lM500_coll))-cbar_lo)/(cbar_hi-cbar_lo)
    #        ax_coll.plot(Sun_lR,Sun_ly,color=cm(Sun_fractional_color),lw=4,ls=":",zorder=0,label="Sun+ (2009)")
    #        for si in range(len(Sun_lR)):
    #            print("Sun2009:", Sun_lR[si], Sun_ly[si], Sun_lylo[si], Sun_lyhi[si])
    #        #ax_coll.plot(Sun_lR,Sun_lylo, color=cm(Sun_fractional_color),ls=":",zorder=0)
    #        #ax_coll.plot(Sun_lR,Sun_lyhi, color=cm(Sun_fractional_color),ls=":",zorder=0)
    #        #ax_coll.fill_between(Sun_lR,Sun_lylo,Sun_lyhi, facecolor=cm(Sun_fractional_color),alpha=0.2,zorder=0)


#     if(lhalo_hi>14.5):

#         CLoGS_lR_coll = []
#         CLoGS_ly_coll = []
#         CLoGS_lM500_coll = []
#         if((property=="entropy_norm")):
#             CLoGS_list = home_name + "CLoGS_entropy_profiles/CLoGS_Table4.dat"

#             CLoGS = open(CLoGS_list,"r").readlines()
#             nCLoGS = len(CLoGS)

#             CLoGS_label = 0

#             for i in range(nCLoGS):
#                 CLoGS_ID = CLoGS[i].split()[0]
#                 CLoGS_M500 = float(CLoGS[i].split()[2])*1e+13
#                 if(property=="entropy_norm"):
#                     CLoGS_entropy_exists = int(CLoGS[i].split()[3])
#                     if(CLoGS_entropy_exists):
#                         R,y = loadCLoGS_entropy(CLoGS_ID)
#                         CLoGS_redshift = 0.0
#                         CLoGS_lM500 = np.log10(CLoGS_M500)
#                         CLoGS_rhocrit = 8.52e-30*(1+CLoGS_redshift)**3 # g/cm^3  
#                         CLoGS_ne500 = 500*f_b*CLoGS_rhocrit/(mu_e*mp) # cm^-3
#                         CLoGS_rho500 = 500*CLoGS_rhocrit  # g/cm^3

#                         CLoGS_R500 = (CLoGS_M500*1.989e+33/CLoGS_rho500*(3/(4*np.pi)))**0.3333   # cm
#                         CLoGS_ne500 = 500*f_b*CLoGS_rhocrit/(mu_e*mp) # cm^-3

#                         CLoGS_T500 = (G*CLoGS_M500*1.989e+33*mu*mp/(2*CLoGS_R500*kboltz))
#                         CLoGS_K500 = CLoGS_T500*CLoGS_ne500**-0.6666*keVperK
#                         lR = np.log10(R/(CLoGS_R500/cmperkpc))
#                         ly = np.log10(y/CLoGS_K500)

#                         CLoGS_lR_coll.extend(lR)
#                         CLoGS_ly_coll.extend(ly)
#                         CLoGS_lM500_coll.extend([CLoGS_lM500])
#                         if(plot_individual_observations):                
#                             if(CLoGS_label == 0):
#                                 if(CLoGS_lM500 > lhalo_lo):
#                                     axes_.plot(lR,ly,color=cm((CLoGS_lM500-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls="--",zorder=5,label="CLoGS Groups")
#                                     CLoGS_label = 1
#                             else:
#                                 if(CLoGS_lM500 > lhalo_lo):
#                                     axes_.plot(lR,ly,color=cm((CLoGS_lM500-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls="--",zorder=5)
#                         else:
#                             axes_.plot(lR,ly,color=cm((CLoGS_lM500-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls="--",zorder=0)

#                     else:
#                         continue

#             if(plot_individual_observations==0):

#                 CLoGS_lR,CLoGS_lRlo,CLoGS_lRhi,CLoGS_ly,CLoGS_lylo,CLoGS_lyhi=flexbin_med_and_onesigmaspread(np.asarray(CLoGS_lR_coll),np.asarray(CLoGS_ly_coll),10)
# #                 print("CLoGS_lR,CLoGS_lRlo,CLoGS_lRhi,CLoGS_ly,CLoGS_lylo,CLoGS_lyhi,median_M500=",CLoGS_lR,CLoGS_lRlo,CLoGS_lRhi,CLoGS_ly,CLoGS_lylo,CLoGS_lyhi,np.median(np.asarray(CLoGS_lM500_coll)))
#                 CLoGS_fractional_color = (np.median(np.asarray(CLoGS_lM500_coll))-cbar_lo)/(cbar_hi-cbar_lo)
#                 axes_.plot(CLoGS_lR,CLoGS_ly,color=cm(CLoGS_fractional_color),lw=4,ls="--",zorder=0,label="CLoGS Groups")
#                 #ax_coll.plot(CLoGS_lR,CLoGS_lylo, color=cm(CLoGS_fractional_color),ls="--",zorder=0)
#                 #ax_coll.plot(CLoGS_lR,CLoGS_lyhi, color=cm(CLoGS_fractional_color),ls="--",zorder=0)
#                 #ax_coll.fill_between(CLoGS_lR,CLoGS_lylo,CLoGS_lyhi, facecolor=cm(CLoGS_fractional_color),alpha=0.2,zorder=0)
                
                
                
                
                
#     if(lhalo_hi>14.5):

    CLoGS_lR_coll = []
    CLoGS_ly_coll = []
    CLoGS_lM500_coll = []
    if property=="entropy_norm":
        CLoGS_list = home_name + "CLoGS_entropy_profiles/CLoGS_Table4.dat"

        CLoGS = open(CLoGS_list,"r").readlines()
        nCLoGS = len(CLoGS)

        CLoGS_label = 0

        for i in range(nCLoGS):
            CLoGS_ID = CLoGS[i].split()[0]
            CLoGS_T = float(CLoGS[i].split()[1])  ## keV
#             CLoGS_M500 = float(CLoGS[i].split()[2])*1e+13
            CLoGS_M500 = scale.M500_Tspec_simbac_N1024L100_renier(pnb.array.SimArray(CLoGS_T, units='keV'))
            
            if CLoGS_M500>=10*lhalo_lo and CLoGS_M500<=10**lhalo_hi: # there is a halo that is getting through this every time
                CLoGS_entropy_exists = int(CLoGS[i].split()[3])
                if(CLoGS_entropy_exists):
                    R,y = loadCLoGS_entropy(CLoGS_ID)
                    CLoGS_redshift = 0.0
                    CLoGS_lM500 = np.log10(CLoGS_M500)
                    CLoGS_rhocrit = 8.52e-30*(1+CLoGS_redshift)**3 # g/cm^3  
                    CLoGS_ne500 = 500*f_b*CLoGS_rhocrit/(mu_e*mp) # cm^-3
                    CLoGS_rho500 = 500*CLoGS_rhocrit  # g/cm^3

                    CLoGS_R500 = (CLoGS_M500*1.989e+33/CLoGS_rho500*(3/(4*np.pi)))**0.3333   # cm
                    CLoGS_ne500 = 500*f_b*CLoGS_rhocrit/(mu_e*mp) # cm^-3

                    CLoGS_T500 = (G*CLoGS_M500*1.989e+33*mu*mp/(2*CLoGS_R500*kboltz))
                    CLoGS_K500 = CLoGS_T500*CLoGS_ne500**-0.6666*keVperK
                    lR = np.log10(R/(CLoGS_R500/cmperkpc))
                    ly = np.log10(y/CLoGS_K500)
                    
                    
                    if not plot_x_log:
                        lR = 10**np.array(lR)
                    if not plot_y_log:
                        ly = 10**np.array(ly)
                    

                    CLoGS_lR_coll.extend(lR)
                    CLoGS_ly_coll.extend(ly)
                    CLoGS_lM500_coll.extend([CLoGS_lM500])
                    
                    if use_cbar and obs_plot_kwargs['Oppenheimer+21']['use_cbar']:
                        if log_cbar:
#                             clogs_color = cm((CLoGS_lM500-np.log10(cbar_lo))/(np.log10(cbar_hi)-np.log10(cbar_lo)))
                            clogs_color = cm((CLoGS_lM500-cbar_lo)/(cbar_hi-cbar_lo))
                        else:
                            clogs_color = cm((10**CLoGS_lM500-cbar_lo)/(cbar_hi-cbar_lo))
                    else:
                        clogs_color = obs_plot_kwargs['Oppenheimer+21']['clogs_color']
                    
#                     if plot_individual_observations:
                    if plot_individual_observations and not obs_plot_kwargs['Oppenheimer+21']['plot_clogs_median']:
                        if(CLoGS_label == 0):
                            if(CLoGS_lM500 > lhalo_lo):
                                axes_.plot(lR,ly,color=clogs_color,lw=1,ls="--",
                                           alpha=obs_plot_kwargs['Oppenheimer+21']['clogs_alpha'], 
                                           zorder=obs_plot_kwargs['Oppenheimer+21']['zorder'],
                                           label="CLoGS Groups 1")
                                CLoGS_label = 1
                        else:
                            if(CLoGS_lM500 > lhalo_lo):
                                axes_.plot(lR,ly,color=clogs_color,lw=1,ls="--",
                                           alpha=obs_plot_kwargs['Oppenheimer+21']['clogs_alpha'],
                                           zorder=obs_plot_kwargs['Oppenheimer+21']['zorder'])
                    elif plot_individual_observations and obs_plot_kwargs['Oppenheimer+21']['plot_clogs_median']:
                        axes_.plot(lR,ly,color=clogs_color,lw=1,ls="--",
                                   alpha=obs_plot_kwargs['Oppenheimer+21']['clogs_alpha'],
                                   zorder=obs_plot_kwargs['Oppenheimer+21']['zorder'])

                else:
                    continue

#         if (plot_individual_observations==0 and lhalo_lo>13.0):
#         if plot_individual_observations==0:
        if (obs_plot_kwargs['Oppenheimer+21']['plot_clogs_median'] and len(CLoGS_lR_coll)>0):

#             if not plot_x_log:
#                 CLoGS_lR_coll = 10**CLoGS_lR_coll
#             if not plot_y_log:
#                 CLoGS_ly_coll = 10**CLoGS_ly_coll
            CLoGS_lR,CLoGS_lRlo,CLoGS_lRhi,CLoGS_ly,CLoGS_lylo,CLoGS_lyhi=flexbin_med_and_onesigmaspread(np.asarray(CLoGS_lR_coll),np.asarray(CLoGS_ly_coll),10)
#                 print("CLoGS_lR,CLoGS_lRlo,CLoGS_lRhi,CLoGS_ly,CLoGS_lylo,CLoGS_lyhi,median_M500=",CLoGS_lR,CLoGS_lRlo,CLoGS_lRhi,CLoGS_ly,CLoGS_lylo,CLoGS_lyhi,np.median(np.asarray(CLoGS_lM500_coll)))
#             CLoGS_fractional_color = (np.median(10**np.asarray(CLoGS_lM500_coll))-cbar_lo)/(cbar_hi-cbar_lo)
            
            if use_cbar and obs_plot_kwargs['Oppenheimer+21']['use_cbar']:
                if log_cbar:
#                     CLoGS_fractional_color = cm((np.median(np.asarray(CLoGS_lM500_coll))-np.log10(cbar_lo))/(np.log10(cbar_hi)-np.log10(cbar_lo)))
                    CLoGS_fractional_color = cm((np.median(np.asarray(CLoGS_lM500_coll))-cbar_lo)/(cbar_hi-cbar_lo))
                else:
                    CLoGS_fractional_color = cm((np.median(10**np.asarray(CLoGS_lM500_coll))-cbar_lo)/(cbar_hi-cbar_lo))
            else:
                CLoGS_fractional_color = obs_plot_kwargs['Oppenheimer+21']['clogs_color']
            
    
            axes_.plot(CLoGS_lR,CLoGS_ly,color=CLoGS_fractional_color,lw=3,ls="--",
                       zorder=obs_plot_kwargs['Oppenheimer+21']['zorder']+5,label="CLoGS Groups 2") # THIS IS BEING PLOTTED
            
            axes_.fill_between(CLoGS_lR, CLoGS_lyhi, CLoGS_lylo, 
                                 color=CLoGS_fractional_color,
                                 alpha=obs_plot_kwargs['Oppenheimer+21']['clogs_alpha'], 
                                 zorder=obs_plot_kwargs['Oppenheimer+21']['zorder'])
            
#             axes_.plot(CLoGS_lR,CLoGS_ly,color=cm(CLoGS_fractional_color),lw=4,ls="--",zorder=0,label="CLoGS Groups")
            #ax_coll.plot(CLoGS_lR,CLoGS_lylo, color=cm(CLoGS_fractional_color),ls="--",zorder=0)
            #ax_coll.plot(CLoGS_lR,CLoGS_lyhi, color=cm(CLoGS_fractional_color),ls="--",zorder=0)
            #ax_coll.fill_between(CLoGS_lR,CLoGS_lylo,CLoGS_lyhi, facecolor=cm(CLoGS_fractional_color),alpha=0.2,zorder=0)
                
                
                
                


#     if(lhalo_hi>14.5):
#         ACCEPT_lR_coll = []
#         ACCEPT_ly_coll = []
#         ACCEPT_lM500_coll = []
#         if((property=="entropy_norm")):
#             ACCEPT_list = home_name + "ACCEPT/ACCEPT.main.tab"

#             ACCEPT = open(ACCEPT_list,"r").readlines()
#             nACCEPT = len(ACCEPT)

#             ACCEPT_label = 0
#             ACCEPT_plot_frequency = 5

#             for i in range(nACCEPT):
#                 if(i<2):
#                     continue # skip header lines.  
#                 ACCEPT_name = ACCEPT[i].split()[0]
#                 ACCEPT_redshift = float(ACCEPT[i].split()[3])
#                 ACCEPT_Tcl = float(ACCEPT[i].split()[7])
#                 ACCEPT_M500 = 3e+13*(ACCEPT_Tcl/1.0)**(1.5)
#                 ACCEPT_lM500 = np.log10(ACCEPT_M500)
#                 if(ACCEPT_lM500<14.0):
#                     continue
#                 #print("ACCEPT_name= ",ACCEPT_name)
#                 if(property=="entropy_norm"):
#                     R,y = loadACCEPTdata(ACCEPT_name,property)
#                     ACCEPT_rhocrit = 8.52e-30*(1+ACCEPT_redshift)**3 # g/cm^3  
#                     ACCEPT_ne500 = 500*f_b*ACCEPT_rhocrit/(mu_e*mp) # cm^-3
#                     ACCEPT_rho500 = 500*ACCEPT_rhocrit  # g/cm^3

#                     ACCEPT_R500 = (ACCEPT_M500*1.989e+33/ACCEPT_rho500*(3/(4*np.pi)))**0.3333   # cm
#                     ACCEPT_ne500 = 500*f_b*ACCEPT_rhocrit/(mu_e*mp) # cm^-3

#                     ACCEPT_T500 = (G*ACCEPT_M500*1.989e+33*mu*mp/(2*ACCEPT_R500*kboltz))
#                     ACCEPT_K500 = ACCEPT_T500*ACCEPT_ne500**-0.6666*keVperK

#                     #print("ACCEPT name, M500, R500, K500,color= ",ACCEPT_name,ACCEPT_lM500,ACCEPT_R500/cmperMpc,ACCEPT_K500,(ACCEPT_lM500-cbar_lo)/(cbar_hi-cbar_lo))#,R,y)
#                     lR = np.log10(R/(ACCEPT_R500/cmperMpc))
#                     ly = np.log10(y/ACCEPT_K500)

#                 ACCEPT_lR_coll.extend(lR[np.where(np.isnan(ly)==0)])
#                 ACCEPT_ly_coll.extend(ly[np.where(np.isnan(ly)==0)])
#                 ACCEPT_lM500_coll.extend([ACCEPT_lM500])
#                 if(plot_individual_observations):                
#                     if(ACCEPT_label == 0):
#                         if(ACCEPT_lM500 > lhalo_lo):                    
#                             if((i-2)%ACCEPT_plot_frequency==0):
#                                 axes_.plot(lR,ly,color=cm((ACCEPT_lM500-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls="-.",zorder=0,label="ACCEPT Clusters")
#                                 ACCEPT_label = 1
#                     else:
#                         if(ACCEPT_lM500 > lhalo_lo):
#                             if((i-2)%ACCEPT_plot_frequency==0):
#                                 axes_.plot(lR,ly,color=cm((ACCEPT_lM500-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls="-.",zorder=0)
#                 else: # Plot a selection of ACCEPT clusters.  
#                     if((i-2)%ACCEPT_plot_frequency==0):
#                         axes_.plot(lR,ly,color=cm((ACCEPT_lM500-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls="-.",zorder=0)

#             if(plot_individual_observations==0):
#                 ACCEPT_lR,ACCEPT_lRlo,ACCEPT_lRhi,ACCEPT_ly,ACCEPT_lylo,ACCEPT_lyhi=flexbin_med_and_onesigmaspread(np.asarray(ACCEPT_lR_coll),np.asarray(ACCEPT_ly_coll),50)
#                 ACCEPT_fractional_color = (np.median(np.asarray(ACCEPT_lM500_coll))-cbar_lo)/(cbar_hi-cbar_lo)
# #                 print("ACCEPT_lR,ACCEPT_lRlo,ACCEPT_lRhi,ACCEPT_ly,ACCEPT_lylo,ACCEPT_lyhi,median_mass=",ACCEPT_lR,ACCEPT_lRlo,ACCEPT_lRhi,ACCEPT_ly,ACCEPT_lylo,ACCEPT_lyhi,np.median(np.asarray(ACCEPT_lM500_coll)))
#                 axes_.plot(ACCEPT_lR,ACCEPT_ly,color=cm(ACCEPT_fractional_color),lw=4,ls="-.",zorder=10,label="ACCEPT Clusters")
#                 #ax_coll.fill_between(ACCEPT_lR,ACCEPT_lylo,ACCEPT_lyhi, facecolor=cm(ACCEPT_fractional_color),alpha=0.2,zorder=10)

        
        
        
        
        
#     if(lhalo_hi>14.5):
    ACCEPT_lR_coll = []
    ACCEPT_ly_coll = []
    ACCEPT_lM500_coll = []
    if property=="entropy_norm":
        ACCEPT_list = home_name + "ACCEPT/ACCEPT.main.tab"

        ACCEPT = open(ACCEPT_list,"r").readlines()
        nACCEPT = len(ACCEPT)

        ACCEPT_label = 0
        ACCEPT_plot_frequency = obs_plot_kwargs['Oppenheimer+21']['accept_frequency'] #4

        for i in range(nACCEPT):
            if(i<2):
                continue # skip header lines.  
            ACCEPT_name = ACCEPT[i].split()[0]
            ACCEPT_redshift = float(ACCEPT[i].split()[3])
            ACCEPT_Tcl = float(ACCEPT[i].split()[7])  # keV
#             ACCEPT_M500 = 3e+13*(ACCEPT_Tcl/1.0)**(1.5)
            ACCEPT_M500 = scale.M500_Tspec_simbac_N1024L100_renier(pnb.array.SimArray(ACCEPT_Tcl, units='keV'))
            ACCEPT_lM500 = np.log10(ACCEPT_M500)
#             if(ACCEPT_lM500<14.0):
            if (ACCEPT_lM500<lhalo_lo or ACCEPT_lM500>lhalo_hi):
                continue
            #print("ACCEPT_name= ",ACCEPT_name)
            if(property=="entropy_norm"):
                R,y = loadACCEPTdata(ACCEPT_name,property)
                ACCEPT_rhocrit = 8.52e-30*(1+ACCEPT_redshift)**3 # g/cm^3  
                ACCEPT_ne500 = 500*f_b*ACCEPT_rhocrit/(mu_e*mp) # cm^-3
                ACCEPT_rho500 = 500*ACCEPT_rhocrit  # g/cm^3

                ACCEPT_R500 = (ACCEPT_M500*1.989e+33/ACCEPT_rho500*(3/(4*np.pi)))**0.3333   # cm
                ACCEPT_ne500 = 500*f_b*ACCEPT_rhocrit/(mu_e*mp) # cm^-3

                ACCEPT_T500 = (G*ACCEPT_M500*1.989e+33*mu*mp/(2*ACCEPT_R500*kboltz))
                ACCEPT_K500 = ACCEPT_T500*ACCEPT_ne500**-0.6666*keVperK

                #print("ACCEPT name, M500, R500, K500,color= ",ACCEPT_name,ACCEPT_lM500,ACCEPT_R500/cmperMpc,ACCEPT_K500,(ACCEPT_lM500-cbar_lo)/(cbar_hi-cbar_lo))#,R,y)
                lR = np.log10(R/(ACCEPT_R500/cmperMpc))
                ly = np.log10(y/ACCEPT_K500)
                
            if not plot_x_log:
                lR = 10**np.array(lR)
            if not plot_y_log:
                ly = 10**np.array(ly)

            ACCEPT_lR_coll.extend(lR[np.where(np.isnan(ly)==0)])
            ACCEPT_ly_coll.extend(ly[np.where(np.isnan(ly)==0)])
            ACCEPT_lM500_coll.extend([ACCEPT_lM500])
            
            if use_cbar and obs_plot_kwargs['Oppenheimer+21']['use_cbar']:
                if log_cbar:
#                     accept_color = cm((ACCEPT_lM500-np.log10(cbar_lo))/(np.log10(cbar_hi)-np.log10(cbar_lo)))
                    accept_color = cm((ACCEPT_lM500-cbar_lo)/(cbar_hi-cbar_lo))
                else:
                    accept_color = cm((10**ACCEPT_lM500-cbar_lo)/(cbar_hi-cbar_lo))
            else:
                accept_color = obs_plot_kwargs['Oppenheimer+21']['accept_color']
            
#             if plot_individual_observations:
            if plot_individual_observations and not obs_plot_kwargs['Oppenheimer+21']['plot_accept_median']:
                if(ACCEPT_label == 0):
                    if(ACCEPT_lM500 > lhalo_lo):                    
                        if((i-2)%ACCEPT_plot_frequency==0):
                            axes_.plot(lR,ly,color=accept_color,lw=1,ls="-.",
                                       alpha=obs_plot_kwargs['Oppenheimer+21']['accept_alpha'],
                                       zorder=obs_plot_kwargs['Oppenheimer+21']['zorder'],
                                       label="ACCEPT Clusters")
                            ACCEPT_label = 1
                else:
                    if(ACCEPT_lM500 > lhalo_lo):
                        if((i-2)%ACCEPT_plot_frequency==0):
                            axes_.plot(lR,ly,color=accept_color,lw=1,ls="-.",
                                       alpha=obs_plot_kwargs['Oppenheimer+21']['accept_alpha'],
                                       zorder=obs_plot_kwargs['Oppenheimer+21']['zorder'])
#             else: # Plot a selection of ACCEPT clusters.
            elif plot_individual_observations and obs_plot_kwargs['Oppenheimer+21']['plot_clogs_median']:
                if((i-2)%ACCEPT_plot_frequency==0):
                    axes_.plot(lR,ly,color=accept_color,lw=1,ls="-.",
                               alpha=obs_plot_kwargs['Oppenheimer+21']['accept_alpha'],
                               zorder=obs_plot_kwargs['Oppenheimer+21']['zorder'])

#         if (plot_individual_observations==0 and lhalo_lo>13.0):
#         if plot_individual_observations==0:
        if (obs_plot_kwargs['Oppenheimer+21']['plot_accept_median'] and len(ACCEPT_lR_coll)>0):
            
#             if not plot_x_log:
#                 ACCEPT_lR_coll = 10**ACCEPT_lR_coll
#             if not plot_y_log:
#                 ACCEPT_ly_coll = 10**ACCEPT_ly_coll
            ACCEPT_lR,ACCEPT_lRlo,ACCEPT_lRhi,ACCEPT_ly,ACCEPT_lylo,ACCEPT_lyhi=flexbin_med_and_onesigmaspread(np.asarray(ACCEPT_lR_coll),np.asarray(ACCEPT_ly_coll),50)
    
#             ACCEPT_fractional_color = (np.median(10**np.asarray(ACCEPT_lM500_coll))-cbar_lo)/(cbar_hi-cbar_lo)
            
            if use_cbar and obs_plot_kwargs['Oppenheimer+21']['use_cbar']:
                if log_cbar:
#                     ACCEPT_fractional_color = cm((np.median(np.asarray(ACCEPT_lM500_coll))-
#                                                   np.log10(cbar_lo))/(np.log10(cbar_hi)-np.log10(cbar_lo)))
                    ACCEPT_fractional_color = cm((np.median(np.asarray(ACCEPT_lM500_coll))-cbar_lo)/(cbar_hi-cbar_lo))
                else:
                    ACCEPT_fractional_color = cm((np.median(10**np.asarray(ACCEPT_lM500_coll))-cbar_lo)/(cbar_hi-cbar_lo))
            else:
                ACCEPT_fractional_color = obs_plot_kwargs['Oppenheimer+21']['accept_color']
                
            axes_.plot(ACCEPT_lR,ACCEPT_ly,color=ACCEPT_fractional_color,lw=3,ls="-.",
                       zorder=obs_plot_kwargs['Oppenheimer+21']['zorder']+5,label="ACCEPT Clusters")
            axes_.fill_between(ACCEPT_lR, ACCEPT_lyhi, ACCEPT_lylo, 
                                 color=ACCEPT_fractional_color,
                                 alpha=obs_plot_kwargs['Oppenheimer+21']['accept_alpha'], 
                                 zorder=obs_plot_kwargs['Oppenheimer+21']['zorder'])
#                 print("ACCEPT_lR,ACCEPT_lRlo,ACCEPT_lRhi,ACCEPT_ly,ACCEPT_lylo,ACCEPT_lyhi,median_mass=",ACCEPT_lR,ACCEPT_lRlo,ACCEPT_lRhi,ACCEPT_ly,ACCEPT_lylo,ACCEPT_lyhi,np.median(np.asarray(ACCEPT_lM500_coll)))
#             axes_.plot(ACCEPT_lR,ACCEPT_ly,color=cm(ACCEPT_fractional_color),lw=4,ls="-.",zorder=10,label="ACCEPT Clusters")
            #ax_coll.fill_between(ACCEPT_lR,ACCEPT_lylo,ACCEPT_lyhi, facecolor=cm(ACCEPT_fractional_color),alpha=0.2,zorder=10)
        
        
        
        
        
        

    if(property=="density"):
        if(plot_only_Lovisari_Sun_groups is not True):
            R,y,ylo,yhi = loadLovisari2015dens()
            lR = np.log10(R)
            ly_ne = np.log10(y)
            lylo_ne = np.log10(ylo)
            lyhi_ne = np.log10(yhi)
            
            if not plot_x_log:
                lR = 10**np.array(lR)
            if not plot_y_log:
                ly_ne = 10**np.array(ly_ne)
                lylo_ne = 10**np.array(lylo_ne)
                lyhi_ne = 10**np.array(lyhi_ne)
            
            axes_.plot(lR,ly_ne,color='black',lw=2,ls="-.",zorder=1,label="Lovisari+2015")
            axes_.plot(lR,lylo_ne,color='black',lw=1,ls="-.",zorder=1)
            axes_.plot(lR,lyhi_ne,color='black',lw=1,ls="-.",zorder=1)

    if(property=="pressure_norm"):
        lR,ly_P,lylo_P,lyhi_P = loadSun2011pressure()
    #    lR = np.log10(R)
    #    ly_P = np.log10(y)
    #    lylo_P = np.log10(ylo)
    #    lyhi_P = np.log10(yhi)
    
        if not plot_x_log:
            lR = 10**np.array(lR)
        if not plot_y_log:
            ly_P = 10**np.array(ly_P)
            lylo_P = 10**np.array(lylo_P)
            lyhi_P = 10**np.array(lyhi_P)
    
        axes_.plot(lR,ly_P,color='black',lw=2,ls="-.",zorder=1,label="Sun+2011")
        axes_.plot(lR,lylo_P,color='black',lw=1,ls="-.",zorder=1)
        axes_.plot(lR,lyhi_P,color='black',lw=1,ls="-.",zorder=1)

    if((property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
        L15_M500, L15_Mgas100kpc, L15_Mgas200kpc, L15_Mgas400kpc = np.loadtxt("Mg_LL15.dat",usecols=(0,1,2,3),unpack=True)
        if(property=="m100kpc"):
            indexes = np.where(L15_Mgas100kpc>0)
            L15_lM500 = np.log10(L15_M500[indexes])
            ly = np.log10(L15_Mgas100kpc[indexes])
        if(property=="m200kpc"):
            indexes = np.where(L15_Mgas200kpc>0)
            L15_lM500 = np.log10(L15_M500[indexes])
            ly =  np.log10(L15_Mgas200kpc[indexes])
        if(property=="m400kpc"):
            indexes = np.where(L15_Mgas400kpc>0)
            L15_lM500 = np.log10(L15_M500[indexes])
            ly = np.log10(L15_Mgas400kpc[indexes])
            
        if not plot_x_log:
            L15_lM500 = 10**np.array(L15_lM500)
        if not plot_y_log:
            ly = 10**np.array(ly)

#         print("Lov15= ", L15_lM500, ly)
        axes_.plot(L15_lM500,ly,color="black",marker="+",zorder=0)
        axes_.plot(L15_lM500,ly,color="black",marker="+",ls="",ms=12,zorder=0,label="Lovisari+2015")


    if(property=="metallicity"):
        R500_low_Lov, R500_hi_Lov, Arelax_Lov, sigrelax_Lov, Adist_Lov, sigdist_Lov = load_Lovisari2019()

        lR_500_Lov = np.log10((R500_low_Lov+R500_hi_Lov)/2.)
        lArelax_Lov = np.log10(Arelax_Lov)
        lAdist_Lov = np.log10(Adist_Lov)
        
        if not plot_x_log:
            lR_500_Lov = 10**np.array(lR_500_Lov)
        if not plot_y_log:
            lArelax_Lov = 10**np.array(lArelax_Lov)
            lAdist_Lov = 10**np.array(lAdist_Lov)

#         print('Lov= ', lR_500_Lov, lArelax_Lov, lAdist_Lov)
        axes_.plot(lR_500_Lov,lArelax_Lov,color='black',marker='o',ms=5,ls='-',label="Lovisari+2019 relaxed",zorder=100000)
        axes_.plot(lR_500_Lov,lAdist_Lov,color='black',marker='s',ms=5,ls=':',label="Lovisari+2019 disturbed",zorder=100000)

        lR_500, Z_3d, Z_3d_upper, Z_3d_lower = load_Sun2009_1550_Z()
        
        if not plot_x_log:
            lR_500 = 10**np.array(lR_500)
        if plot_y_log:
            Z_3d = np.log10(np.array(Z_3d))
            Z_3d_lower = np.log10(np.array(Z_3d_lower)+np.array(Z_3d))-np.log10(np.array(Z_3d))
            Z_3d_upper = -np.log10(np.array(Z_3d))+np.log10(np.array(Z_3d)+np.array(Z_3d_upper))
            
        if use_cbar and obs_plot_kwargs['Oppenheimer+21']['use_cbar']:
            if log_cbar:
#                 metal_color = cm((np.log10(3.18e+13)-np.log10(cbar_lo))/(np.log10(cbar_hi)-np.log10(cbar_lo)))
                metal_color = cm((np.log10(3.18e+13)-cbar_lo)/(cbar_hi-cbar_lo))
            else:
                metal_color = cm((3.18e+13-cbar_lo)/(cbar_hi-cbar_lo))
        else:
            metal_color = 'black'
            

        axes_.errorbar(lR_500, np.log10(Z_3d),yerr=[Z_3d_lower,Z_3d_upper],color=metal_color,lw=1,ls="--",zorder=10000,label="Sun+2009 NGC 1550 ")


    if(property=="XraySB_XMMcounts"):
        load_and_plot_Lovisari_Sun_groups(property, axes_, cbar_lo, cbar_hi, cm=cm, plot_x_log=plot_x_log, plot_y_log=plot_y_log)

    if((property=="LX")):
        M500,M500_err,LX,LX_err = loadAndersonWang()
        
        if plot_x_log:
            M500 = np.log10(np.array(M500))
            M500_err = [np.abs(np.array(M500) - np.log10(10**np.array(M500)-np.array(M500_err))), np.abs(np.array(M500) - np.log10(10**np.array(M500)+np.array(M500_err)))]
        if plot_y_log:
            LX = np.log10(np.array(LX))
            LX_err = [np.abs(np.array(LX) - np.log10(10**np.array(LX)-np.array(LX_err))), np.abs(np.array(LX) - np.log10(10**np.array(LX)+np.array(LX_err)))]
        
        axes_.plot(M500,LX,marker="o",ls='',color='k',label="Wang+ (2016)")
        axes_.errorbar(M500,LX,xerr=M500_err,yerr=LX_err,ls='',color='k')

    # if((property=="fgas") | (property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc") | (property=="LX")):

    #     if(property=="fgas"):
    #         M500_Xray, fgas500 = np.loadtxt("McCarthy17_fgas.txt", usecols=(0,1), unpack=True)

    #         axes_.plot(M500_Xray, fgas500,color="black",marker="x",ls='',label="Observations")
    #         axes_.plot(M500_Xray, fgas500,color="black",marker="x",ls='',label="Observations")

    #     ax_all.set_title(title_label,fontsize=18)
    #     ax_all.tick_params(axis='both', which='major', labelsize=12)
    #     ax_all.set_xlim(lMlo,lMhi)
    #     ax_all.set_ylim(lylo,lyhi)
    #     ax_all.set_xlabel("log $M_{500}$ [$M_{\odot}$]",fontsize=16)
    #     ax_all.set_ylabel(ylabel,fontsize=16)
    #     if(multi==0):
    #         fig_all.subplots_adjust(left=leftborder, bottom=0.17,top=0.90,right=0.98)
    #         fig_all.savefig('profile_all_' + file_base + '.' + property  + mhalo_str + '.png')
    # else:    
    #     #ax_all.set_xscale('log')
    #     #ax_all.set_yscale('log')
    #     ax_all.set_title(title_label,fontsize=18)
    #     ax_all.tick_params(axis='both', which='major', labelsize=12)

    #     ax_all.set_xlim(lRlo,lRhi)
    #     if(property=="dispersionmeasure"):
    #         ax_all.set_ylim(ylo,yhi)
    #         ax_all.set_xlabel("log $R$ [kpc]",fontsize=16)
    #     else:
    #         ax_all.set_ylim(lylo,lyhi)
    #         ax_all.set_xlabel("log $R/R_{500}$",fontsize=16)
    #     ax_all.set_ylabel(ylabel,fontsize=16)
    #     sm = plt.cm.ScalarMappable(cmap=cm, norm=mpl.colors.Normalize(vmin=cbar_lo, vmax=cbar_hi))
    #     sm._A = []
    #     div = axgrid.make_axes_locatable(ax_all)
    #     cbar = fig_all.colorbar(sm)
    #     cbar.ax.set_xlabel(cbarlabel,labelpad=8,fontsize=12)

    #     #if(multi==0):
    #     fig_all.subplots_adjust(left=leftborder, bottom=0.17,top=0.90,right=0.98)
    #     fig_all.savefig('profile_all_' + file_base + '.' + property  + mhalo_str + '.png')



    #ax_coll.set_xscale('log')
    #ax_coll.set_yscale('log')
    # if(multi < 2):
    #     ax_coll.set_title(title_label,fontsize=18)

    # if(legend>0):
    #     if(multi==2):
    #         ax_coll.legend(loc="best",fontsize=12)
    #     else:
    #         ax_coll.legend(loc="best",fontsize=12)
    #     ax_coll.tick_params(axis='both', which='major', labelsize=12)

    # if(property=="dispersionmeasure"):
    #     ax_coll.set_ylim(ylo,yhi)
    # else:
    #     ax_coll.set_ylim(lylo,lyhi)

    # if((property=="fgas") | (property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc") | (property=="LX")):
    #     plt.locator_params(axis="x", nbins=5)
    #     ax_coll.set_xlim(lMlo,lMhi)
    #     ax_coll.set_xlabel("log $M_{500}$ [$M_{\odot}$]",fontsize=16)
    # else:
    #     ax_coll.set_xlim(lRlo,lRhi)
    #     if(property=="dispersionmeasure"):
    #         ax_coll.set_xlabel("log $R$ [kpc]",fontsize=16)
    #     else:
    #         ax_coll.set_xlabel("log $R/R_{500}$",fontsize=16)



    # if((property=="temperature") | ((property=="dispersionmeasure") & (multi!=1)) | ((property=="entropy_norm") & (multi==2))):
    #     sm = plt.cm.ScalarMappable(cmap=cm, norm=mpl.colors.Normalize(vmin=cbar_lo, vmax=cbar_hi))
    #     sm._A = []
    #     div = axgrid.make_axes_locatable(ax_coll)
    #     cbar = fig_coll.colorbar(sm)
    #     cbar.ax.set_xlabel(cbarlabel,labelpad=8,fontsize=12)

    # ax_coll.set_ylabel(ylabel,fontsize=16)
    # fig_coll.subplots_adjust(left=leftborder, bottom=0.12,top=topborder,right=rightborder)
    # fig_coll.savefig('profile_coll_' + file_base + '.' + property  + mhalo_str + '.png')

    # axes_.legend()
    print("finished")

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



from palettable.cartocolors.sequential import DarkMint_7, BurgYl_7, Purp_7, Sunset_7, agGrnYl_7
COLOR_LIST = [DarkMint_7.hex_colors[1:], BurgYl_7.hex_colors[1:], Purp_7.hex_colors[1:],
              Sunset_7.hex_colors[1:], agGrnYl_7.hex_colors[1:],
              ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', 'darkgreen', 'black', 'pink'],
              ["#3a405a","#7483a3","#aec5eb","#bdc0d9","#ccbac7","#dbb5b5","#e9afa3","#a98074","#685044"]]

# CALL run_plotting_script(...) within THIS SCRIPT TO PLOT
# scaling_relation_functions is a script that aviv has


# plot_individual_observations=0 or 1
# plot_only_Lovisari_Sun_groups = true or False
# above 2 parameters, mess with them. one configuration will plot ALL data, one will plot some, and one will plot less

# MIGHT have to play with obs_plot_kwargs a bit??
# plot_sim_profiles.py in avivs shows how to set

# line 3436 in group_xray_profiles_observations_plotting_functions.py
# log_bin_low = 12.8, log_bin_high = 13.2

# to get what obs_plot_kwargs looks like
# plot_profiles_simba.py
# big cell, line 909-919, might have to adjust to use new matplotlib 
#           make a dictionary that has line 2081-2102 in it. i can remove 2097,2098, 2100, and others

oppenheimer_dict = {
    'Oppenheimer+21': {
        'plot_individual_observations': True, # OG = True
        'plot_only_Lovisari_Sun_groups': False, # OG = False
        'sun_lovisari_color': 'tab:orange',
        'sun_lovisari_alpha': 0.25,
        'plot_sun_median': False, # OG = True
        'clogs_color': 'tab:purple',
        'clogs_alpha': 0.25,
        'plot_clogs_median': False, # OG = True
        'accept_color': 'tab:green',
        'accept_alpha': 0.25,
        'accept_frequency': 1,
        'plot_accept_median': False, # OG = True
        'plot_just_entropy_theory': False,
        'use_cbar': False,
        'entropy_theory_zorder': 120,
        'zorder': 50,
        'verbose': True
    }
}

cm_ = plt.get_cmap('viridis') ## og

plot_kwargs = { 
    'xaxis':{
        'is_log':True,

    },
    'yaxis':{
        'is_log': True,
        },

    'cbar':{
        'plot':False,
        'cm':cm_, 
        'log_cbar':True,
    },
}


# FROM plot_profiles_simba.py

cmap = mpl.colors.ListedColormap([COLOR_LIST[5][1::2][0], COLOR_LIST[5][1::2][1], COLOR_LIST[5][1::2][2]])
bounds = [13, 13.5, 14, 15]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
sm_ = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)

obs_plot_kwargs=oppenheimer_dict
axes=plt

run_plot = run_plotting_script(
    property='entropy_norm',
    lhalo_lo=low_mass_range,
    lhalo_hi=high_mass_range,
    axes_=axes,
    sm=sm_,
    cm=plot_kwargs['cbar']['cm'],
    plot_individual_observations=obs_plot_kwargs['Oppenheimer+21']['plot_individual_observations'], 
    plot_only_Lovisari_Sun_groups=obs_plot_kwargs['Oppenheimer+21']['plot_only_Lovisari_Sun_groups'], 
    plot_x_log=plot_kwargs['xaxis']['is_log'],
    plot_y_log=plot_kwargs['yaxis']['is_log'], 
    use_cbar=plot_kwargs['cbar']['plot'], 
    log_cbar=plot_kwargs['cbar']['log_cbar'],
    obs_plot_kwargs=obs_plot_kwargs
)


# -------------------------------------------------------------------------------------------------------------------------------------- SIM_DATA
import yt
import unyt
from unyt import unyt_quantity, unyt_array, km, s, Mpc, kg, m, Msun, kpc
from unyt import Gyr, Mpc, km, s
from unyt.physical_constants import G as G_unyt, kboltz, mp as mp_unyt  # G, k_B, m_p

import numpy as np
import matplotlib.pyplot as plt
import caesar
import pynbody as pnb
import pynbody.units as u
import pynbody.array as pna

import dill
import os
from astro_constants import NA_no_units, NA, kB_no_units, kB, mH, mp_no_units, mp, mu_e, mu, G_no_units, G
#from astro_constants import NA_no_units, NA, kB_no_units, kB, mH, mp_no_units, mp, mu_e, mu, G_no_units,

code = 'Simba'

snap_file = '/scratch/b/babul/aspadawe/snapshots/Hyenasc_L1_zooms/halo_' + str(zoom_ID) + '/snap_halo_' + str(zoom_ID) + '_151.hdf5'

# ADD CAESAR FILEPATH HERE TOO

snap = yt.load(snap_file)#, unit_system='cgs')
s = pnb.load(snap_file)


## Metal indices
if code.lower() == 'simba':
    # Simba
    metals_idx = {'Z':0, 'He': 1, 'C': 2, 'N': 3, 'O': 4, 'Ne': 5, 'Mg': 6, 'Si': 7, 'S': 8, 'Ca': 9, 'Fe': 10}
    metals_idx_xigrm = {'H':0, 'He': 1, 'C': 2, 'N': 3, 'O': 4, 'Ne': 5, 'Mg': 6, 'Si': 7, 'S': 8, 'Ca': 9, 'Fe': 10}
    elems = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca', 'Fe']
elif code.lower() == 'simba-c':
    # Simba-C
    metals_idx = {'Z':0, 'H':1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al':13,
                  'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19,'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24,
                  'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30}
    metals_idx_xigrm = {'H':0, 'He': 1, 'Li': 2, 'Be': 3, 'B': 4, 'C': 5, 'N': 6, 'O': 7, 'F': 8, 'Ne': 9, 'Na': 10, 'Mg': 11, 
                        'Al':12, 'Si': 13, 'P': 14, 'S': 15, 'Cl': 16, 'Ar': 17, 'K': 18,'Ca': 19, 'Sc': 20, 'Ti': 21, 'V': 22, 'Cr': 23,
                      'Mn': 24, 'Fe': 25, 'Co': 26, 'Ni': 27, 'Cu': 28, 'Zn': 29}
    elems=['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si',
           'P', 'S', 'Cl', 'Ar', 'K','Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']#, 'Ga', 'Ge']\n",
else:
    raise SystemExit('Code is not Simba or Simba-C')


# Get data, and some functions

data = snap.all_data()

def gas_entropy(field, data):
    return data["gas", "kT"] * data["gas", "El_number_density"] ** (-2.0/3.0)

def gas_entropy2(field, data):
    return data["gas", "kT"] * data["PartType0", "El_number_density"] ** (-2.0/3.0)

snap.add_field(
    ("gas", "gas_entropy"),
    units="keV*cm**2",
    function=gas_entropy,
    sampling_type="cell",
)

snap.add_field(
    ("PartType0", "gas_entropy"),
    units="keV*cm**2",
    function=gas_entropy2,
    sampling_type="particle",
    force_override=True
)

gas_entropy = data['gas', 'entropy']


caesarfilepath = '/scratch/b/babul/aspadawe/snapshots/Hyenasc_L1_zooms/halo_' + str(zoom_ID) + '/Groups/Caesar_halo_' + str(zoom_ID) + '_151.hdf5'
obj = caesar.load(caesarfilepath)
halo_data = obj.halos


# ------------------------------------------------------------------------------------------------------------------ Loop!!!
color_count = 0
for target_halo_id in halo_ID_arr:

    halo_center = halo_data[target_halo_id].minpotpos
    halo_center = halo_center.in_units('Mpc')
    print(halo_center)


    igrm = data.include_above(('PartType0', 'temperature'), 5e5, units='K')
    igrm = igrm.include_below(('PartType0', 'H_nuclei_density'), 0.13, units='cm**-3')
    igrm = igrm.include_equal(('PartType0', 'DelayTime'), 0)


    # Radial Binning

    # radial_bins = unyt.unyt_array(np.append([0], np.logspace(-3, np.log10(5), 30)), 'Mpc')
    radial_bins = unyt.unyt_array(np.append([0], np.logspace(-4, np.log10(5), 45)), 'Mpc')

    halo_center.in_units('Mpc')

    igrm_pos_x = igrm['gas', 'particle_position_x']
    igrm_pos_y = igrm['gas', 'particle_position_y']
    igrm_pos_z = igrm['gas', 'particle_position_z']

    # isolate igrm

    igrm_halo_center = ((igrm_pos_x-halo_center[0])**2 + (igrm_pos_y-halo_center[1])**2 + (igrm_pos_z-halo_center[2])**2)**0.5

    igrm_zoom_entropies = igrm['gas', 'entropy']
    igrm_zoom_masses = igrm['gas', 'mass']


    entropy_profile_median = []
    entropy_profile_mean = []
    entropy_profile_mass_weighted = []
    for bin_lo, bin_hi in zip(radial_bins[:-1], radial_bins[1:]):
    #     print(bin_lo, bin_hi)
        radial_filter = np.logical_and(igrm_halo_center.in_units('Mpc').value>=bin_lo.in_units('Mpc').value,
                                    igrm_halo_center.in_units('Mpc').value<=bin_hi.in_units('Mpc').value)
    #     print(np.where(radial_filter))
        entropies = igrm_zoom_entropies[radial_filter]
        masses = igrm_zoom_masses[radial_filter]
        # print(len(entropies))
        # print(len(masses))
    #     print()
        
        entropy_profile_median.append(np.nanmedian(entropies))
        entropy_profile_mean.append(np.nanmean(entropies))
        if sum(masses)>0:
            K_mass = np.average(entropies, weights=masses)
            entropy_profile_mass_weighted.append(K_mass)
        else:
            entropy_profile_mass_weighted.append(np.nan)
        
        # print(K_mass)
        # print()



    from unyt import Msun, kpc, cm, keV, unyt_quantity
    from unyt.physical_constants import mp, kboltz, G  # ✅ use unyt's version


    # Constants
    mu = 0.59        # Mean molecular weight
    mu_e = 1.14      # Mean molecular weight per free electron

    # Load CAESAR catalog
    h = obj.simulation.hubble_constant
    omega_b = obj.simulation.omega_baryon
    omega_m = obj.simulation.omega_matter

    # Update f_b if desired
    f_b = omega_b / omega_m

    # Critical density
    H0_si = (100 * h * 1e3) / (3.086e22)  # H0 in 1/s
    rho_crit_si = 3 * H0_si**2 / (8 * np.pi * 6.6743e-11)  # kg/m^3
    rho_crit = unyt_quantity(rho_crit_si, "kg/m**3").to("Msun/kpc**3")

    # Choose your halo
    halo = obj.halos[target_halo_id]  # Replace with desired index or search method

    M500 = unyt_quantity(halo.virial_quantities['m500c'].in_units("Msun"), "Msun")
    M200 = unyt_quantity(halo.virial_quantities['m200c'].in_units("Msun"), "Msun")

    R500 = unyt_quantity(halo.virial_quantities['r500c'].in_units("m"), "m")
    R200 = unyt_quantity(halo.virial_quantities['r200c'].in_units("m"), "m")

    ne_500 = (500 * f_b * rho_crit / (mu_e * mp)).to("cm**-3")
    ne_200 = (200 * f_b * rho_crit / (mu_e * mp)).to("cm**-3")

    T500 = (G * M500 * mu * mp / R500).to("keV") # removing 1/kb here
    T200 = (G * M200 * mu * mp / R200).to("keV") # removing 1/kb here

    K500 = (T500 * ne_500**(-2/3)).to("keV*cm**2") # removing kb here
    K200 = (T200 * ne_200**(-2/3)).to("keV*cm**2") # removing kb here

    R500 = halo_data[target_halo_id].virial_quantities['r500c'].in_units('Mpc') # redefine r500 and r200 with correct units before binning
    R200 = halo_data[target_halo_id].virial_quantities['r200c'].in_units('Mpc')
    radial_bin_centres = (radial_bins[:-1] + radial_bins[1:])/2.
    radial_bin_centres_logged = np.log10(radial_bin_centres)
    log_radial_bin_centres = (np.log10(radial_bins[:-1]) + np.log10(radial_bins[1:]))/2.
    log_radial_bin_centres_unlogged = 10**log_radial_bin_centres

    normalized_r = radial_bin_centres/R500.value
    normalized_k_mean = np.array(entropy_profile_mean)/K500.value
    normalized_k_median = np.array(entropy_profile_median)/K500.value


    # axes_.plot(CLoGS_lR,CLoGS_ly,color=CLoGS_fractional_color,lw=6,ls="--", zorder=obs_plot_kwargs['Oppenheimer+21']['zorder']+5,label="CLoGS Groups 1401")

    log_r = np.log10(normalized_r)
    log_k = np.log10(normalized_k_median)

    color_arr = ['blue', 'cyan', 'orange', 'red', 'pink']


    plt.plot(
        log_r,
        log_k,
        color=color_arr[color_count],
        lw=2,
        ls='--',
        zorder=60,
        label='Halo of M500: ' + str(round(np.log10(halo_data[target_halo_id].virial_quantities['m500c'].value), 1)) + ', Cont: ' + str(round(halo_data[target_halo_id].contamination, 3))
    )

    color_count+=1


# ------------------------------------------------------------------------------------------------------------------ end Loop!!!



note = 'Observational Data Mass Range [M500]: ' + str(low_mass_range) + ' - ' + str(high_mass_range)

plt.text(
    0.98, 0.02,                # x, y (in axes fraction coordinates)
    note,            # your note
    transform=plt.gca().transAxes,  # place relative to axis (not data)
    fontsize=8, 
    ha='right', 
    va='bottom'
)


plt.title(thetitle)
plt.xlabel(r'log($R / R_{500}$)')
plt.ylabel(r'$\log(K_{\mathrm{mass}} / K_{500})$')
plt.legend(loc='upper left', fontsize=7)



plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.savefig(savefilename)