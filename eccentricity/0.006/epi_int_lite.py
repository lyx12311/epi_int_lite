#!/usr/bin/env python

#epi_int_lite.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 27 August 2017.
#
#this simulates the dynamical evolution of self-gravitating planetary rings.
#The majority of this code was written while drinking and singing at
#the Water Tank karaoke bar in northwest Austin TX, so buyer beware.

#read input parameters
import numpy as np
execfile('inputs.py')
print 'number_of_streamlines =', number_of_streamlines
print 'particles_per_streamline =', particles_per_streamline
print 'dt =', dt
print 'timesteps_per_output = ', timesteps_per_output
print 'total_number_of_outputs =', total_number_of_outputs
print 'radial_width =', radial_width
print 'total_ring_mass =', total_ring_mass
print 'ring gravitation constant =', G_ring
print 'shear_viscosity =', shear_viscosity
print "Toomre's Q_ring =", Q_ring
print 'Rp =', Rp
print 'J2 =', J2
print 'initial_orbits =', initial_orbits
print 'output_folder =', output_folder

#initialize orbits
from helper_fns import *
r, t, vr, vt, lambda0, c = initialize_streamline(number_of_streamlines, particles_per_streamline,
    radial_width, total_ring_mass, G_ring, Q_ring, shear_viscosity, J2, Rp, initial_orbits)

#prep for main loop
timestep = 0
number_of_outputs = 0
(rz, tz, vrz, vtz, timestepz) = ([r], [t], [vr], [vt], [timestep])
import time as tm
clock_start = tm.time()

#evolve system...this largely follows Chamber's (1993) 2nd order drift-kick scheme but assumes
#the central mass has negligable motion about system's center-of-mass ie the ring is nearly
#circular and there are no point-mass satellites such that Chamber's exp(tau*C/2)=1
print 'evolving system...'
while (number_of_outputs < total_number_of_outputs):
    #kick velocities forwards by timestep +dt/2
    vr, vt = kick(J2, Rp, lambda0, G_ring, shear_viscosity, c, r, t, vr, vt, dt/2.0)
    timesteps_since_output = 0
    while (timesteps_since_output < timesteps_per_output):
        #convert coordinates to elements
        a, e, wt, M = coords2elem(J2, Rp, r, t, vr, vt)
        #advance mean anomaly during drift step
        M = drift(a, M, J2, Rp, dt)
        #convert orbit elements to coordinates
        r, t, vr, vt = elem2coords(J2, Rp, a, e, wt, M)
        #kick velocities
        vr, vt = kick(J2, Rp, lambda0, G_ring, shear_viscosity, c, r, t, vr, vt, dt)
        #updates
        timestep += 1
        timesteps_since_output += 1
        #print timestep
    #kick velocities backwards by timestep -dt/2
    vr, vt = kick(J2, Rp, lambda0, G_ring, shear_viscosity, c, r, t, vr, vt, -dt/2.0)
    #save output
    number_of_outputs += 1
    rz, tz, vrz, vtz, timestepz = \
        store_system(rz, tz, vrz, vtz, timestepz, r, t, vr, vt, timestep)
    run_time_min = (tm.time() - clock_start)/60.0
    eta_min = int((total_number_of_outputs - number_of_outputs)*run_time_min/number_of_outputs)
    if (20*number_of_outputs%total_number_of_outputs == 0):
        print 'time = ' + str(timestep*dt) + \
            '    number of outputs = ' + str(number_of_outputs) + \
            '    number of orbits = ' + str(int(timestep*dt/2.0/np.pi)) + \
            '    eta (minutes) = ', eta_min

#save results
timez = np.array(timestepz)*dt
save_output(rz, tz, vrz, vtz, timez, lambda0, output_folder)
print 'execution time (minutes) = ', (tm.time() - clock_start)/60.0
