#libration.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 18 March 2018.
#
#these helper functions are used to analyze an evolving ringlet

#imports
import numpy as np
from helper_fns import *

#compute BGT's H(q) function
def H(q):
    q2 = q**2
    q_factor = np.sqrt(1 - q2)
    H = (1 - q_factor)/q2/q_factor
    return H

#calculate da, de, dwt differences at inner & outer streamline's periapse
def calculate_Deltas(r, a, e, wt):

    #calculate delta_a
    total_number_of_outputs = r.shape[0]
    number_of_streamlines = r.shape[1]
    particles_per_streamline = r.shape[2]
    a_inner = a[:, 0, :].mean(axis=1)
    a_outer = a[:, number_of_streamlines-1, :].mean(axis=1)
    a_avg = (a_inner + a_outer)/2
    delta_a = a_outer - a_inner
    
    #compute delta_e and e_prime = a_mid*delta_e/delta_a versus time
    e_inner = e[:, 0, :].mean(axis=1)
    e_outer = e[:, number_of_streamlines-1, :].mean(axis=1)
    e_avg = (e_inner + e_outer)/2
    delta_e = e_outer - e_inner
    e_prime = a_avg*delta_e/delta_a
    delta_e_avg = (delta_e.max() - delta_e.min())/2
    e_prime_avg = a_avg*delta_e_avg/delta_a
    
    #compute delta_wt, wt_prime = a_mid*e_mid*delta_wt/delta_a, and q
    delta_wt_list = []
    for t_idx in range(total_number_of_outputs):
        s_idx = 0
        r0 = r[t_idx, s_idx]
        theta_idx = np.argmin(r0)
        wt_inner = wt[t_idx, s_idx, theta_idx]
        s_idx = number_of_streamlines - 1
        r0 = r[t_idx, s_idx]
        theta_idx = np.argmin(r0)
        wt_outer = wt[t_idx, s_idx, theta_idx]
        delta_wt_list += [wt_outer - wt_inner]
    delta_wt = adjust_angle(np.array(delta_wt_list))
    delta_wt_avg = (delta_wt.max() + delta_wt.min())/2
    wt_prime = a_avg*e_avg*delta_wt/delta_a
    wt_prime_avg = a_avg*e_avg*delta_wt_avg/delta_a

    #compute q and H(q)
    q = np.sqrt(e_prime**2 + wt_prime**2)
    Hq = H(q)
    
    return a_inner, a_outer, a_avg, delta_a, e_inner, e_outer, e_avg, \
        delta_e, delta_e_avg, e_prime, e_prime_avg, delta_wt, delta_wt_avg, \
        wt_prime, wt_prime_avg, q, Hq
