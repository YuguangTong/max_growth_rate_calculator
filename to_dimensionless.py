import os
import numpy as np
from py_vlasov.util import (pmass, emass, echarge, 
                   permittivity, permeability, cspeed, boltzmann)

data_dir = '../chris_small_data'
full_data_path = os.path.abspath(data_dir)

b = np.loadtxt(full_data_path + '/B.txt')
ne = np.loadtxt(full_data_path + '/NE.txt')
np2 = np.loadtxt(full_data_path + '/NP2.txt')
np1 = np.loadtxt(full_data_path + '/NP1.txt')
na = np.loadtxt(full_data_path + '/NA.txt')
te = np.loadtxt(full_data_path + '/TE.txt')
tp2 = np.loadtxt(full_data_path + '/TP2.txt')
tp1 = np.loadtxt(full_data_path + '/TP1.txt')
ta = np.loadtxt(full_data_path + '/TA.txt')
ve = np.loadtxt(full_data_path + '/VE.txt')
vp2 = np.loadtxt(full_data_path + '/VP2.txt')
vp1 = np.loadtxt(full_data_path + '/VP1.txt')
va = np.loadtxt(full_data_path + '/VA.txt')

# b in tesla
def f_beta_p1(b, np1, tp1):
    b_si = b * 1e-9
    np1_si = np1 * 1e6
    if type(tp1) == np.ndarray:
        kb_tp1_par_si = tp1[:, 0] * echarge
    else:
        kb_tp1_par_si = tp1 * echarge
    return b_si**2 / (2 * permeability * np1_si * kb_tp1_par_si)

# beta_p1
beta_sample = f_beta_p1(b, np1, tp1)

def f_va1(b, np1):
    b_si = b * 1e-9
    np1_si = np1 * 1e6
    return b_si/np.sqrt(permeability * pmass * np1_si)

# Alfven speed of the dominant proton component
va1 = f_va1(b, np1)

# drift in the frame of the dominant proton component,
# normalized by alfven speed
def f_vd(vp1, vs, ns):
    vds = (vs-vp1)
    vds[ns==0] = [0, 0, 0]
    vds_norm = np.array([np.linalg.norm(v) for v in vds])
    return vds_norm * 1e3

norm_vd_p1 = f_vd(vp1, vp1, np1)/va1
norm_vd_p2 = f_vd(vp1, vp2, np2)/va1
norm_vd_e = f_vd(vp1, ve, ne)/va1
norm_vd_ap = f_vd(vp1, va, na)/va1

# temperature anisotropy 
def f_a(t):
    tperp = t[:, 1] 
    tpar = t[:, 0]
    a = tperp/tpar
    a[np.isnan(a)] = 1
    return a

a_e = f_a(te)
a_p1 = f_a(tp1)
a_p2 = f_a(tp2)
a_ap = f_a(ta)

# temperature ratio
def f_tratio(tp1, ts):
    tp1_par = tp1[:, 0]
    ts_par = ts[:, 0]
    tratio = ts_par/tp1_par
    tratio[np.isnan(tratio)] = 1
    return tratio

t_sample = np.array([f_tratio(tp1, tp1),
                     f_tratio(tp1, tp2),
                     f_tratio(tp1, ta),
                     f_tratio(tp1, te)]).T

# density
n_sample = np.array([np1/np1, np2/np1, na/np1, ne/np1]).T
m_sample = [1, 1, 4, 1/1836.]
q_sample = [1, 1, 2, -1]
a_sample = np.array([a_p1, a_p2, a_ap, a_e]).T
v_sample = np.array([norm_vd_p1, norm_vd_p2, norm_vd_ap, norm_vd_e]).T

np.savez('dimensionless_input', beta = beta_sample, m = m_sample, q = q_sample,
        t=t_sample, a = a_sample, n = n_sample, v = v_sample)
