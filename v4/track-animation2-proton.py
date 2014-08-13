#!/usr/bin/env python
from math import *
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import BasicFunc as func
import Input as para

###########################################################################################
num_of_particles = para.num_of_particles # set no. of particles of a beam for tracking
num_of_turns = para.app4_num_of_turns+1 # track no. of turns 
sigma_dPoP = para.sigma_dPoP # set the sigma of (+/-) Delta_P/P
mean_dPoP = para.mean_dPoP # set the mean of (+/-) Delta_P/P
range_dPoP = para.range_dPoP # define the survival range of Delta_P/P
range_phi1 = para.range_phi1 # define the lower limit of survival range phi (rad) for protons
range_phi2 = para.range_phi2 # define the upper limit of survival range phi (rad) for protons
###########################################################################################

var_t_tmp = 0.0 # set initial ramping time = 0 (s)
var_t = np.zeros(num_of_particles) # start ramping from t=0 (sec) 
var_E = func.E_total_p(var_t_tmp)*np.ones(num_of_particles)

var_E_tmp = func.E_total_p(var_t_tmp)
var_beta2_tmp = func.beta2_p(var_E_tmp)
var_beta2 = func.beta2_p(var_E_tmp)*np.ones(num_of_particles)

# assume the dE distribution is the gaussian with the mean and sigma
np.random.seed(12345)
var_dE = mean_dPoP + sigma_dPoP*np.random.randn(num_of_particles)*var_E_tmp*var_beta2_tmp

# assume the phi is randomly distributed between [-pi, +pi]
np.random.seed(34567) 
var_phi = np.pi*2.0*np.random.random(num_of_particles)-np.pi

show_phi = np.zeros(num_of_particles) 
show_dPoP = np.zeros(num_of_particles) 

#
# animation settings
#
FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='RF phase space', artist='Mark C.C. Chiang', comment=' ')
writer = FFMpegWriter(fps=25, metadata=metadata)

#
# plot settings
#
fig = plt.figure()
l, = plt.plot([], [], 'ro', markeredgecolor = 'none')

##################################################
plt.xlim(para.set_xlim1, para.set_xlim2)
plt.ylim(para.set_ylim1, para.set_ylim2)
##################################################
plt.xlabel('$\phi$ (rad)', fontsize=20)
plt.ylabel('$\Delta P / P $ (%)', fontsize=20)

#
# set text position
#
ax = plt.axes()

#
# make the animation
#
show_eff = 9999.0*np.ones(num_of_turns)
show_turn = 9999*np.ones(num_of_turns)
resolution = 100
with writer.saving(fig, "track-animation2-proton.mp4", resolution):
    for i in range(num_of_turns):
        count = 0
        eff = 0.0
        time = 0.0
        for j in range(num_of_particles):

            show_phi[j] = var_phi[j]
            show_dPoP[j] = var_dE[j]/var_E[j]/var_beta2[j]

            var_dE[j], var_phi[j] = func.iteration_p(var_dE[j], var_phi[j], var_t[j], var_E[j])

            var_t[j] = func.t_p_new(var_t[j], var_E[j])
            time_tmp = var_t[j]

            var_E[j] = func.E_total_p(var_t[j])
            var_beta2[j] = func.beta2_p(var_E[j])

            if (range_phi1<=show_phi[j]<=range_phi2 and abs(show_dPoP[j])<=range_dPoP):
                count +=1

        eff = 100.0*count/num_of_particles
        time = time_tmp
        print 'time= ', time, ' ; ref. capture rate (%)= ', eff

        if (i<=1000 and i%5==0):
            l.set_data(show_phi, 100.0*show_dPoP)
            ax.set_title('$%3.3f$ ms' %(time*1000), fontsize=28)
            writer.grab_frame()
        elif (1000<i<=2000 and i%10==0):
            l.set_data(show_phi, 100.0*show_dPoP)
            ax.set_title('$%3.3f$ ms' %(time*1000), fontsize=28)
            writer.grab_frame()        
        elif (i>2000 and i%200==0):
            l.set_data(show_phi, 100.0*show_dPoP)
            ax.set_title('$%3.3f$ ms' %(time*1000), fontsize=28)
            writer.grab_frame()

