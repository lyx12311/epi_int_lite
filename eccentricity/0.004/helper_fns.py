#helper_fns.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 28 September 2017.
#
#helper functions used by epi_int_lite.py

#angular frequency
def Omega(J2, Rp, a, Ar=0.0):
    GM = 1.0
    a2 = a*a
    Ra2 = (Rp*Rp)/a2
    Omega2 = (GM/a2/a)*(   1.0 + (1.5*J2)*Ra2 - Ar*(a2/GM)   )
    return np.sqrt(Omega2)

#epicyclic frequency
def Kappa(J2, Rp, a, Ar=0.0, kappa_squared=False):
    GM = 1.0
    a2 = a*a
    Ra2 = (Rp*Rp)/a2
    Kappa2 = (GM/a2/a)*(   1.0 - (1.5*J2)*Ra2 - Ar*(a2*(3.0/GM))   )
    if (kappa_squared):
        return Kappa2
    else:
        return np.sqrt(Kappa2)

#adjust angles to live between -Pi and Pi
import numpy as np
def adjust_angle(angle):
    idx = angle > np.pi
    angle[idx] -= 2.0*np.pi
    idx = angle < -np.pi
    angle[idx] += 2.0*np.pi    
    return angle

#unwrap angle that lives between -Pi and Pi
def unwrap_angle(angle):
    angle_unw = angle.copy()
    twopi = 2.0*np.pi
    for idx in range(1, len(angle)):
        delta = angle_unw[idx] - angle_unw[idx - 1]
        if (delta < -np.pi):
            angle_unw[idx:] += twopi 
        if (delta > np.pi):
            angle_unw[idx:] -= twopi 
    return angle_unw

#drift step advances M
def drift(a, M, J2, Rp, dt):
    return M + Kappa(J2, Rp, a)*dt

#use lagrange polynomial to evaluate function f that is evaluated n adjacent
#streamlines away and sampled at longitude t
def interpolate_fn(t, f, n, interpolate=True):
    if (interpolate):
        t0 = np.roll(t, (-n,  1), axis=(0,1))
        t1 = np.roll(t, (-n,  0), axis=(0,1))
        t2 = np.roll(t, (-n, -1), axis=(0,1))
        f0 = np.roll(f, (-n,  1), axis=(0,1))
        f1 = np.roll(f, (-n,  0), axis=(0,1))
        f2 = np.roll(f, (-n, -1), axis=(0,1))
        f_n = lagrange_poly_fit(t0, t1, t2, f0, f1, f2, t)
    else:
        #skip lagrange interpolation
        f_n = np.roll(f, (-n, 0), axis=(0,1))
    return f_n

#fit 2nd order lagrange polynomial to data (x0,y0),(x1,y1),(x2,y2) & interpolate y(x)
def lagrange_poly_fit(x0, x1, x2, y0, y1, y2, x):
    l0 = ((x - x1)/(x0 - x1))*((x - x2)/(x0 - x2))
    l1 = ((x - x0)/(x1 - x0))*((x - x2)/(x1 - x2))
    l2 = ((x - x0)/(x2 - x0))*((x - x1)/(x2 - x1))
    y = y0*l0 + y1*l1 + y2*l2
    return y

#wrap the ring's coordinate array about in longitude
def wrap_ring(c, longitude=False):
    Nr, Nt = c.shape
    left = c[:, -1].copy().reshape(Nr, 1)
    right = c[:, 0].copy().reshape(Nr, 1)
    if (longitude):
        twopi = 2*np.pi
        left -= twopi
        right += twopi
    cw = np.concatenate((left, c, right), axis=1)
    return cw

#compute ring surface density
def surface_density(lambda0, dr):
    sd = lambda0/dr
    return sd

#calculate radial distance between exterior streamline and interior streamline
def delta_f(f, t):
    #in ring interior
    f_plus  = interpolate_fn(t, f, -1)
    f_minus = interpolate_fn(t, f,  1)
    if (f.shape[0] > 2):
        df = f_plus - f_minus
    else:
        df = np.zeros_like(f)
    #along ring edges
    df[0] = f_plus[0] - f[0]
    df[-1] = f[-1] - f_minus[-1]
    return df

#calculation derivative df/dr
def df_dr(delta_f, delta_r):
    return delta_f/delta_r

#radial acceleration due to ring self-gravity
def ring_gravity(lambda0, G_ring, r, t):
    two_G_lambda = 2.0*G_ring*lambda0
    Ar = np.zeros_like(r)
    Nr, Nt = r.shape
    for shft in range(1, Nr):
        dr = interpolate_fn(t, r, -shft) - r
        Ar += two_G_lambda/dr
    return Ar

#acceleration due to pressure P
def A_P(lambda0, sd, P, t, delta_P, delta_r):
    dPdr = df_dr(delta_P, delta_r)
    #acceleration in ring interior
    A = -dPdr/sd
    #at inner streamline
    A[0] = -P[0]/lambda0[0]
    #at outer streamline, interpolated from neighbor streamline
    P_outer = interpolate_fn(t[-2:], P[-2:], 1)[-1]
    A[-1] = P_outer/lambda0[-1]
    return A

#radial acceleration due to ring pressure
def ring_pressure(c, lambda0, sd, r, t, delta_r):
    P = (c*c)*sd
    delta_P = delta_f(P, t)
    Ar = A_P(lambda0, sd, P, t, delta_P, delta_r)
    return Ar

#tangential acceleration due to ring viscosity
def ring_viscosity(shear_viscosity, lambda0, sd, r, t, vt, delta_r):
    w = vt/r
    delta_w = delta_f(w, t)
    dw_dr = df_dr(delta_w, delta_r)
    #viscous pseudo-pressure
    P = (-shear_viscosity*sd)*r*dw_dr
    delta_P = delta_f(P, t)
    At = A_P(lambda0, sd, P, t, delta_P, delta_r)
    return At

#calculate radial and tangential accelerations due to ring gravity, pressure, visocisty
def accelerations(lambda0, G_ring, shear_viscosity, c, r, t, vt):
    #wrap ring around in longitude
    rw = wrap_ring(r, longitude=False)
    tw = wrap_ring(t, longitude=True)
    lw = wrap_ring(lambda0, longitude=False)
    Ar = np.zeros_like(rw)
    At = np.zeros_like(rw)
    #radial acceleration due to streamline gravity
    if (G_ring > 0.0):
        Ar += ring_gravity(lw, G_ring, rw, tw)
    #radial acceleration due to streamline pressure and viscosity
    if ((c > 0.0) or (shear_viscosity > 0.0)):
        delta_rw = delta_f(rw, tw)
        sdw = surface_density(lw, delta_rw)
        if (c > 0.0):
            Ar += ring_pressure(c, lw, sdw, rw, tw, delta_rw)
        if (shear_viscosity > 0.0):
            vtw = wrap_ring(vt, longitude=False)
            At += ring_viscosity(shear_viscosity, lw, sdw, rw, tw, vtw, delta_rw)
    #drop left and right edges from Ar,At
    Ar = Ar[:, 1:-1]
    At = At[:, 1:-1]
    return Ar, At

#velocity kicks due to ring gravity and viscosity
def kick(J2, Rp, lambda0, G_ring, shear_viscosity, c, r, t, vr, vt, dt):
    #radial acceleration due to ring gravity and pressure
    Ar, At = accelerations(lambda0, G_ring, shear_viscosity, c, r, t, vt)
    #kick velocity
    vr += Ar*dt
    vt += At*dt
    return vr, vt

#convert orbit elements to coordinates
def elem2coords(J2, Rp, a, e, wt, M, Ar=0.0, sort_particle_longitudes=True):
    e_sin_M = e*np.sin(M)
    e_cos_M = e*np.cos(M)
    r = a/np.sqrt(1.0 + 2.0*e_cos_M)
    Omg = Omega(J2, Rp, a, Ar=Ar)
    Kap = Kappa(J2, Rp, a, Ar=Ar)
    t = adjust_angle(   (Omg/Kap)*(M + 2.0*e_sin_M) + wt   )
    ra = r/a
    ra3 = ra*ra*ra
    vr = (a*Kap*ra3)*e_sin_M
    vt = a*a*Omg/r
    #sort each streamline's particles by longitude as needed
    if (sort_particle_longitudes):
        r, t, vr, vt = sort_particles(r, t, vr, vt)
    return r, t, vr, vt

#convert coordinates to orbit elements
def coords2elem(J2, Rp, r, t, vr, vt, Ar=0.0):
    GM = 1.0
    h = r*vt
    c = (h*h)/(2.0*GM*Rp)
    a = Rp*(   c + np.sqrt(c*c - 1.5*J2)   )
    Omg = Omega(J2, Rp, a, Ar=Ar)
    Kap = Kappa(J2, Rp, a, Ar=Ar)
    ar = a/r
    ar2 = ar*ar
    ar3 = ar2*ar
    e_cos_M = (ar2 - 1.0)/2.0
    aK = a*Kap
    e_sin_M = vr*(ar3/aK)
    e = np.sqrt(e_sin_M*e_sin_M + e_cos_M*e_cos_M)
    M = np.arctan2(e_sin_M, e_cos_M)
    wt = adjust_angle(   t - (Omg/Kap)*(M + 2.0*e_sin_M)   )
    return a, e, wt, M

#order particles in each streamline by their longitudes
def sort_particles(r, t, vr, vt):
    for streamline_idx in range(len(t)):
        longitude_idx = t[streamline_idx].argsort()
        r[streamline_idx] = r[streamline_idx][longitude_idx]
        t[streamline_idx] = t[streamline_idx][longitude_idx]
        vr[streamline_idx] = vr[streamline_idx][longitude_idx]
        vt[streamline_idx] = vt[streamline_idx][longitude_idx]
    return r, t, vr, vt

#append current r,t,vr,vt,a,timestep to lists rz,tz etc
def store_system(rz, tz, vrz, vtz, timestepz, r, t, vr, vt, timestep):
    rz.append(r)
    tz.append(t)
    vrz.append(vr)
    vtz.append(vt)
    timestepz.append(timestep)
    return rz, tz, vrz, vtz, timestepz

#save orbit element arrays in files
def save_output(r, t, vr, vt, times, lambda0, output_folder):
    import os
    cmd = 'mkdir -p ' + output_folder
    q = os.system(cmd)
    np.save(output_folder + '/r.npy', r)
    np.save(output_folder + '/t.npy', t)
    np.save(output_folder + '/vr.npy', vr)
    np.save(output_folder + '/vt.npy', vt)
    np.save(output_folder + '/times.npy', times)
    np.save(output_folder + '/lambda0.npy', lambda0)

#restore orbit elements from files
def restore_output(output_folder):
    r = np.load(output_folder + '/r.npy')
    t = np.load(output_folder + '/t.npy')
    vr = np.load(output_folder + '/vr.npy')
    vt = np.load(output_folder + '/vt.npy')
    times = np.load(output_folder + '/times.npy')
    lambda0 = np.load(output_folder + '/lambda0.npy')
    return r, t, vr, vt, times, lambda0

#initialize streamlines
def initialize_streamline(number_of_streamlines, particles_per_streamline, radial_width,
    total_ring_mass, G_ring, Q_ring, shear_viscosity, J2, Rp, initial_orbits):
    
    #initialize particles in circular orbits
    a_streamlines = np.linspace(1.0, 1.0 + radial_width, num=number_of_streamlines)
    a_list = []
    for sma in a_streamlines:
        a_list.append(np.zeros(particles_per_streamline) + sma)
    a = np.array(a_list)
    e = np.zeros_like(a)
    #particles anomalies are uniformly spaced
    M_streamline = np.linspace(-np.pi, np.pi, num=particles_per_streamline, endpoint=False)
    M_list = [M_streamline]*number_of_streamlines
    M = np.array(M_list)
    #tweak longitude of periapse away from zero so that any eccentric streamlines are closed loops
    Omg = Omega(J2, Rp, a)
    Kap = Kappa(J2, Rp, a)
    wt = -(Omg/Kap - 1)*M
    
    #modify initial orbits as needed
    if (initial_orbits['shape'] == 'circular'):
        pass
    if (initial_orbits['shape'] == 'eccentric'):
        e_init = initial_orbits['e']
        adeda = initial_orbits['adeda']
        e = e_init + adeda*(a - a[0])/a[0]
    if (initial_orbits['shape'] == 'breathing mode'):
        e_init = initial_orbits['e']
        e[:] = e_init
        M[:] = 0.0
    if (initial_orbits['shape'] == 'log-e'):
        #streamlines' e is lograthmically distributed between initial_e[0] & initial_e[1] with random M,wt
        wt = np.random.uniform(low=-np.pi, high=np.pi, size=a.shape)
        M = np.random.uniform(low=-np.pi, high=np.pi, size=a.shape)
        e = np.zeros_like(M)
        e_init = initial_orbits['e']
        for idx in range(number_of_streamlines):
            e[idx] += np.exp( np.random.uniform(low=np.log(e_init[0]), high=np.log(e_init[1])) )
    
    #lambda0=streamline mass-per-lenth
    mass_per_streamline = total_ring_mass/number_of_streamlines
    twopi = 2.0*np.pi
    lambda0 = mass_per_streamline/(twopi*a)
    if (total_ring_mass > 0):
        print('this lambda-check should equal one = ', \
            (lambda0[:,0]*twopi*a_streamlines).sum()/total_ring_mass)
    
    #calculate ring sound speed c
    r, t, vr, vt = elem2coords(J2, Rp, a, e, wt, M)
    delta_r = delta_f(r, t)
    sd = surface_density(lambda0, delta_r)
    G = 1.0
    c = (Q_ring*np.pi*G*sd/Omg).mean()
    
    #convert elements to coordinates
    Ar, At = accelerations(lambda0, G_ring, shear_viscosity, c, r, t, vt)
    #r, t, vr, vt = elem2coords(J2, Rp, a, e, wt, M, Ar=Ar) #causes jitter in librating ringlets
    r, t, vr, vt = elem2coords(J2, Rp, a, e, wt, M)
    
    return r, t, vr, vt, lambda0, c

#recompute coordinates in coordinate system that co-rotates with ringlet's middle streamline's peri
def peri_corotate(r, t, vr, vt, wt):
    number_of_streamlines = r.shape[0]
    s_idx = (number_of_streamlines - 1)/2
    r_middle_streamline = r[s_idx]
    t_idx = np.argmin(r_middle_streamline)
    wt_middle_streamline = wt[s_idx]
    wt_min = wt_middle_streamline[t_idx]
    tw = adjust_angle(t - wt_min)
    wts = adjust_angle(wt - wt_min)
    rs, ts, vrs, vts = sort_particles(r, tw, vr, vt)
    return rs, ts, vrs, vts, wts
