#!/usr/bin/python

"""
Using Nick's "NORMALIZED" SHW Model. Please see the writeup.
"""

import random as rd
import time
import math

rd.seed(time.time())

RUNS    = 10        # Number of runs
W       = 24        # Periodicity
KMAX    = 10000     # Maximum iterations
STEPMEAN= 1e-3      # Mean step size
TEMP0   = 1e10      # Initial temperature
COOLING = 0.9999    # Cooling ratio

def calc_fcast(t, alpha, beta, gamma):
    if len(t) <= 2 * W:
        print "TS too short for any forecast!"
        return None

    l = sum(t[W:2*W]) / float(W)
    b = sum([t[i+W] - t[i] for i in range(W)]) / (float(W) ** 2)
    s = [t[i+W] - l for i in range(W-1)]      # is a list of length (W-1)
    
    fcast = []

    for i in range(len(t) - 2*W):
        y = t[2*W+i]
        # forecast
        f = l + b + s[0]
        fcast.append(f)
        # l, b, s evolution
        l_new = (1.0 - alpha) * l   + (1.0 - alpha) * b     + (-1.0 * alpha) * s[0]     + alpha * y
        b_new = (-1.0 * beta ) * l  + (1.0 - beta) * b      + (-1.0 * beta) * s[0]      + beta * y
        s_new = [0.0] * (W-1)
        for j in range(W-2):
            s_new[j] = (gamma/W) * l + (gamma/W) * b + (gamma/W) * s[0] + s[j+1] - (gamma/W) * y
        s_new[W-2] = (gamma/W) * l + (gamma/W) * b + (gamma/W - 1) * s[0] - sum(s[1:])
        # Update
        l = l_new
        b = b_new
        s = (s_new)

    return fcast
    



def calc_mse(ts1, ts2):
    if len(ts1) != len(ts2):
        print "calc_mse: ts1 and ts2 not aligned!"
        return -1

    mse = 0.0
    ts_len = len(ts1)
    for t1, t2 in zip(ts1, ts2):
        mse += ((t1 - t2) ** 2) / ts_len
    
    return mse




def get_neighbor(alpha, beta, gamma):
    valid = 0
    ret_a = alpha
    ret_b = beta
    ret_g = gamma

    STEP = rd.random() * 2 * STEPMEAN

    while valid == 0:
        r = rd.randint(1,6) # Roll a dice
        if r == 1:
            if ret_a + STEP < 1.0:
                ret_a += STEP
                valid = 1
        if r == 2:
            if ret_b + STEP < 1.0:
                ret_b += STEP
                valid = 1
        if r == 3:
            if ret_g + STEP < 1.0:
                ret_g += STEP
                valid = 1
        if r == 4:
            if ret_a - STEP > 0.0:
                ret_a -= STEP
                valid = 1
        if r == 5:
            if ret_b - STEP > 0.0:
                ret_b -= STEP
                valid = 1
        if r == 6:
            if ret_g - STEP > 0.0:
                ret_g -= STEP
                valid = 1

    return ret_a, ret_b, ret_g




def sim_anneal(ts, alpha, beta, gamma):
    ret_a = alpha
    ret_b = beta
    ret_g = gamma
    temp  = TEMP0
    update = 0

    fcast = calc_fcast(ts, ret_a, ret_b, ret_g)
    energy = calc_mse(ts[2*W:], fcast)

    for k in range(KMAX):
        #if k % int(KMAX/20.0) == 0:
        #    print "Iteration %d: energy = %E, alpha = %f, beta = %f, gamma = %f" %(k, energy, ret_a, ret_b, ret_g)
        temp = COOLING * temp
        a_nbr, b_nbr, g_nbr = get_neighbor(ret_a, ret_b, ret_g) # Get neightbor state
        fcast_nbr = calc_fcast(ts, a_nbr, b_nbr, g_nbr)
        energy_nbr = calc_mse(ts[2*W:], fcast_nbr)

        if energy_nbr < energy:
            energy = energy_nbr
            fcast = fcast_nbr
            ret_a, ret_b, ret_g = a_nbr, b_nbr, g_nbr
            update = k+1
        else:
            transition_prob = math.exp((energy - energy_nbr) / temp)
            rd_num = rd.random()
            if rd_num <= transition_prob:
                energy = energy_nbr
                fcast = fcast_nbr
                ret_a, ret_b, ret_g = a_nbr, b_nbr, g_nbr
                update = k+1

    return energy, fcast, ret_a, ret_b, ret_g, update 


ts = []
in_file = open("./ts.txt")
for line in in_file:    ts.append(float(line.rstrip()))
out_file = open("./opt.csv", "w")

lowest_energy = float("inf")
best_fcast = None

for i in range(RUNS):
    alpha   = rd.random()   # Random initial values. Subject to change
    beta    = rd.random()
    gamma   = rd.random()
    print "%-2d-th run:" %(i),

    st_time = time.time()
    energy, fcast, a_opt, b_opt, g_opt, update = sim_anneal(ts, alpha, beta, gamma)
    el_time = time.time() - st_time

    print "energy = %.6e, alpha = %.6f, beta = %.6f, gamma = %.6f, iterations = %-6d, elapsed time = %.3f" %(energy, a_opt, b_opt, g_opt, update, el_time) 
    out_file.write("%d,%.20f,%.20f,%.20f,%.20f,%d\n" %(i, energy, a_opt, b_opt, g_opt, update))

    if energy < lowest_energy:  best_fcast = fcast

out_file.close()

out_file = open("./compare.csv", "w")
cnt = 0
for ts1, ts2 in zip(ts[2*W:], best_fcast):
    out_file.write("%d,%.20e,%.20e,%.20e\n" %(cnt, ts1, ts2, ts1-ts2))
    cnt += 1
