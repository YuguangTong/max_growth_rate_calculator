# provide function(s) to find a start frequency for a specified mode
# in a plasma

import numpy as np
import scipy
from py_vlasov.wrapper import disp_det
from py_vlasov.util import real_imag, list_to_complex
from py_vlasov.new_follow_parameter import (follow_kz, 
                                            follow_beta,
                                            follow_temperature,
                                            follow_anisotropy,
                                            follow_density,
                                            follow_drift)
import matplotlib.pyplot as plt

def find_freq(param, wave_mode='pfh', show_plot=False):
    """
    parameters
    ----------
    param = [target_kz, target_kp, target_beta,
             target_t, target_a, target_n,
             target_q, target_m, target_v,
             n, method, target_aol]
    wave_mode:
              'pfh'--> parallel firehose (right hand)
              'ic'--> parallel ion-cyclotron mode (left hand)
              'ofh'--> oblique firehose
              'm'--> mirror

    return
    ------
    frequency of WAVE_MODE in the plasma specified by PARAM

    N.B. 
    1, This wrapper assumes that there are four species: proton core,
    proton halo, alpha particles and electrons.
    2, For now, waves are assumed to propagate in the positive z direction.
    To loosen this constraint, we can flip the drift velocity.
    """
    # unpack plasma parameters for target plasmas
    (target_kz, target_kp, target_beta, target_t, target_a, target_n,
     target_v, target_aol) = param

    
    # default initial parameters
    betap = 1.
    t_list=[1., 1., 1., 1.]
    a_list=[1., 1., 1., 1.]
    n_list=[1.,0., 0.0, 1.0] 
    q_list=[1., 1, 2, -1.]
    m_list=[1., 1., 4., 1/1836]
    v_list=[0, 0, 0, 0]
    method = 'numpy'
    aol = 1/5000.
    n = 10

    if wave_mode == 'pfh': # parallel firehose
        kp, kz = 0, 0.1
        seed_freq = 0.1
    elif wave_mode == 'ic': # ion cyclotron
        kp, kz = 0, 0.1
        seed_freq = 0.1 # MODIFY HERE
    elif wave_mode == 'ofh': # ion cyclotron
        kp, kz = 0.1, 0.1
        seed_freq = 0.1 # MODIFY HERE
    elif wave_mode == 'm': # ion cyclotron
        kp, kz = 0.1, 0.1
        seed_freq = 0.1 # MODIFY HERE
    else:
        raise Exception('Do not recognize the wave mode:{0}'.format(wave_mode))

    # follow density
    param = [kz, kp, betap, t_list, a_list, n_list, q_list,
             m_list, v_list, n, method, aol]
    lin_inc = 0.1
    freq1, param1, freq_lst = follow_density(
        seed_freq, target_n, param,
        lin_incrmt=lin_inc, show_plot = show_plot)
    
    # follow along beta
    seed_freq = freq1
    log_inc = 0.05
    freq2, param2, freq_lst = follow_beta(
        seed_freq, target_beta, param1,
        log_incrmt=log_inc, incrmt_method='log', show_plot = show_plot)

    # follow along anisotropy
    seed_freq = freq2
    log_inc = 0.05
    freq3, param3, freq_lst = follow_anisotropy(
        seed_freq, target_a, param2, log_incrmt=log_inc,
        incrmt_method = 'log', show_plot = show_plot)

    # follow along temperature ratio
    seed_freq = freq3
    log_inc = 0.05
    freq4, param4, freq_lst = follow_temperature(
        seed_freq, target_t, param3, log_incrmt=log_inc,
        incrmt_method = 'log', show_plot = show_plot)

    # follow along drift
    seed_freq = freq4
    lin_inc = 0.1
    freq5, param5, freq_lst = follow_drift(
        seed_freq, target_v, param4,
        lin_incrmt=lin_inc, show_plot = show_plot)

    # follow along kz
    seed_freq = freq5
    log_inc = 0.05
    freq6, param6, freq_lst= follow_kz(
        seed_freq, target_kz, param, show_plot=show_plot,
        log_incrmt=log_inc, incrmt_method = 'log')


    return freq6
