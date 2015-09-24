#!/usr/bin/python

import pyvaropt as pv
import sys
import os
import time
import numpy as np
import netaddr as na

RSVR_SIZE= 30000
FN_PREF  = "./test_data/prefix_all.txt"
DATA_DIR = "/home/users/cgi0911/Data/Waikato_5/hourly_flowbin/"
RES_DIR  = "/home/users/cgi0911/Results/Waikato_5/hourly_flowbin/%s/" %(time.strftime("%Y%m%d-%H%M%S", time.localtime()))
INTERVAL = 3600         # Seconds in a time slot
TS_START = 1181088000   # Starting timestamp (in seconds)
TS_END   = TS_START + INTERVAL * 16 * 24    # Ending timestamp
FILETYPE = "flowbin"
PERIOD   = 24    # # of time slots in a period
R        = 1    # Forecast # of steps
ALPHA    = 0.2
BETA     = 0.2
GAMMA    = 0.2

# ---------- Global variables and objects ----------
TS_CURR  = 0                # Current timestamp
x_vec    = None             # State vector. Make it globally accessible (shared dictionary)
y        = None             # Current time slot's observation
m_mat    = None             # Transition matrix
u_vec    = None             # U-Vector for transition

RES_DIR_TRUEFCAST = os.path.join(RES_DIR, "true_fcast")
if not os.path.exists(RES_DIR_TRUEFCAST):   os.makedirs(RES_DIR_TRUEFCAST)




def make_trans_matrix():
    ret = []
    for i in range(PERIOD + 1):
        ret.append([0.0] * (PERIOD+1))

    ret[0][0] = 1.0 - ALPHA
    ret[0][1] = 1.0 - ALPHA
    ret[0][2] = -1.0 * ALPHA
    ret[1][0] = -1.0 * BETA
    ret[1][1] = 1.0 - BETA
    ret[1][2] = -1.0 * BETA

    for i in range(PERIOD - 2):
        ret[2+i][0] = (GAMMA / float(PERIOD))
        ret[2+i][1] = (GAMMA / float(PERIOD))
        ret[2+i][2] = (GAMMA / float(PERIOD))
        ret[2+i][3+i] = 1.0

    ret[PERIOD][0] = (GAMMA / float(PERIOD))
    ret[PERIOD][1] = (GAMMA / float(PERIOD))
    ret[PERIOD][2] = (GAMMA / float(PERIOD)) - 1.0
    for i in range(PERIOD + 1 - 3):
        ret[PERIOD][3+i] = -1.0

    ret = np.array(ret)     # Make it a numpy array
    return ret




def make_u_vec():
    ret = [-1.0 * GAMMA / float(PERIOD)] * (PERIOD+1)
    ret[0] = ALPHA
    ret[1] = BETA
    
    ret = np.array(ret).reshape(PERIOD+1, 1)
    return ret




def read_rec(ts, pfx):
    table = pv.KWTable(filetype=FILETYPE, fn=os.path.join(DATA_DIR, str(ts)+".rec"))
    pfx.query(table)
    pfx_val = pfx.get_values()
    return np.array(pfx_val).reshape(1, len(pfx_val))
    



def forecast():
    # r is the # of steps to look ahead.
    one_pos = 1 + (R % PERIOD)
    f_vec = [1.0, float(R)] + [0.0] * (PERIOD-1)
    f_vec[one_pos] = 1.0
    f_mat = np.array(f_vec).reshape(1, PERIOD+1)
    #print f_mat.shape, x_vec.shape
    ret = np.dot(f_mat, x_vec).flatten()
    el_time = time.time() - st_time
    print "Making %d-step forecast. Time stamp = %d. Elapsed time = %f" %(R, TS_CURR + INTERVAL * (R-1), el_time)
    return ret




def transition():
    #print m_mat.shape, x_vec.shape, u_vec.shape, y.shape
    print x_vec
    ret = np.dot(m_mat, x_vec) + np.dot(u_vec, y)
    return ret




def to_txt_file(vec, pfx, fn):
    keys = pfx.get_keys()
    pfx_size = len(pfx)
    vec = vec.flatten()
    outfile = open(fn, "w")
    for i in range(pfx_size):
        key_fields = keys[i].split(' ')
        nw = key_fields[0]
        sd = key_fields[1]         
        outfile.write("%s,%s,%f\n" %(sd, nw, vec[i]))
    outfile.close()




if __name__ == "__main__":
    pfx = pv.PrefixQueryTable(fn=FN_PREF)
    pfx_size = len(pfx)
    l0 = np.zeros(pfx_size).reshape(1, pfx_size)
    b0 = np.zeros(pfx_size).reshape(1, pfx_size)
    s0_list = [np.zeros(pfx_size).reshape(1, pfx_size) for i in range(PERIOD)]
                                                        # create (w-1) empty vectors
                                                        # note that s0_list[0] is not used
                                                        # indices 1..(PERIOD-1) are used

    print
    print "-" * 80
    print "Using Python interpreter:", sys.executable
    print "TS_START = %d, TS_END = %d" %(TS_START, TS_END)
    print "INTERVAL = %d, PERIOD = %d intervals" %(INTERVAL, PERIOD)
    print "(ALPHA, BETA, GAMMA) = (%f, %f, %f)" %(ALPHA, BETA, GAMMA)
    print "-" * 80
    print

    # ---------- Sample-based Seasonal Holt-Winters Initialization ----------
    # (1) l0 = \sum\limits_{i=0}^{W-1} Y_{-i}/W
    # (2) b0 = \sum\limits_{i=0}^{w-1} (Y_{-i} - Y_{-W-i}) / W^2
    # (3) S_{0,i} = Y_{i-W} - l0


    # ---------- First training period ----------
    print "---------- Initialization: First training period ----------"
    TS_CURR = TS_START
    for i in range(PERIOD):
        st_time = time.time()
        y = read_rec(TS_CURR, pfx)
        b0 = b0 - y/float(PERIOD)**2
        print "Read time slot", TS_CURR, "   b0 -= y/W^2",
        TS_CURR += INTERVAL
        el_time = time.time() - st_time
        print "   Elapsed time =", el_time


    # ---------- Second training period ----------
    print
    print "---------- Initialization: Second training period ----------"
    for i in range(PERIOD):
        st_time = time.time()
        y = read_rec(TS_CURR, pfx)
        l0 = l0 + y/float(PERIOD)
        b0 = b0 + y/float(PERIOD)**2
        print "Read time slot", TS_CURR, "   l0 += y/W",
        print "   b0 += y/W^2",
        if i < (PERIOD-1):   # Filling S_{0,1}, S_{0,2}, ... , S_{0,W-1}
            s0_list[i+1] = s0_list[i+1] + y
            print "   s0_list[%d] += y" %(i+1),
        TS_CURR += INTERVAL
        el_time = time.time() - st_time
        print "   Elapsed time =", el_time


    for i in range(1, PERIOD):
        st_time = time.time()
        s0_list[i] = s0_list[i] - l0
        print "s0_list[%d] -= l0" %(i),
        el_time = time.time() - st_time
        print "   Elapsed time =", el_time

    print l0.shape
    print b0.shape
    for s0 in s0_list: print s0.shape

    # ---------- Form state vector ----------
    print
    print "--------- Form initial state vector and sample each element ----------"
    x_vec = np.array([l0.flatten(), b0.flatten()] + [s0.flatten() for s0 in s0_list[1:]])
    print "---------- End of initialization ----------"

    # ---------- Recurrence: Transition and forecast ----------
    print
    print "---------- Start of recursion ----------"
    print "Make transition matrix..."
    m_mat = make_trans_matrix()
    print m_mat
    print "Make U-vector..."
    u_vec = make_u_vec()
    print u_vec

    while TS_CURR < TS_END:
        st_time_recur = time.time()
        # ---------- Forecast based on existing state vector ----------
        print
        print "## Iteration of time %d ##" %(TS_CURR)
        fcast = forecast()     # Make one step forecast.
        fcast_fn = os.path.join(RES_DIR_TRUEFCAST, "%d.txt" %(TS_CURR + (R-1)*INTERVAL))
        to_txt_file(fcast, pfx, fcast_fn)

        # ---------- Transition ----------
        # First read in current time slot's records
        y = read_rec(TS_CURR, pfx)
        print "Read time slot %d, size = %d" %(TS_CURR, len(y))

        # Then do the transition of current state vector
        print "Transition of x_vec..."
        x_vec = transition()

        TS_CURR += INTERVAL
