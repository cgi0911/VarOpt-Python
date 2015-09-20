#!/usr/bin/python

import pyvaropt as pv
import sys
import os
import time
import multiprocessing as mp

RSVR_SIZE= 30000
DATA_DIR = "/home/users/cgi0911/Data/Waikato_5/hourly_flowbin/"
RES_DIR  = "/home/users/cgi0911/Results/Waikato_5/hourly_flowbin/%s/" %(time.strftime("%Y%m%d-%H%M%S", time.localtime()))
INTERVAL = 3600         # Seconds in a time slot
TS_START = 1181088000   # Starting timestamp (in seconds)
TS_END   = TS_START + INTERVAL * 60    # Ending timestamp
FILETYPE = "flowbin"
PERIOD   = 24   # # of time slots in a period
R        = 1    # Forecast # of steps
ALPHA    = 0.2
BETA     = 0.2
GAMMA    = 0.2
N_WORKERS= 4

# ---------- Global variables and objects ----------
TS_CURR  = 0                # Current timestamp
x_vec    = []               # State vector. Make it globally accessible
y        = pv.KWTable()     # Current time slot's observation
m_mat    = []               # Transition matrix
u_vec    = []               # U-Vector for transition

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

    return ret




def make_u_vec():
    ret = [-1.0 * GAMMA / float(PERIOD)] * (PERIOD+1)
    ret[0] = ALPHA
    ret[1] = BETA
    return ret




def read_rec(ts):
    return pv.KWTable(filetype=FILETYPE, fn=os.path.join(DATA_DIR, str(ts)+".rec"))




def forecast(r):
    # r is the # of steps to look ahead.
    st_time = time.time()
    ret = pv.KWTable()
    l = x_vec[0]
    b = x_vec[1]
    s = x_vec[r+1]  # r-th seasonal component
    ret.aggr_inplace(l, 1.0, 1.0)
    ret.aggr_inplace(b, 1.0, float(r))
    ret.aggr_inplace(s, 1.0, 1.0)
    ret.rsvr_sample(RSVR_SIZE, in_place=True)
    el_time = time.time() - st_time
    print "Making %d-step forecast. Time stamp = %d. Elapsed time = %f" %(R, TS_CURR, el_time)
    return ret




def worker_samp(i): # Worker function of sampling
    # Working on x_vec[i] as assigned by parent
    st_time = time.time()
    print "Worker is working on x_vec[%d]: KWTable(%-12x). Current size = %d" %(i, id(x_vec[i]), len(x_vec[i])),
    ret = x_vec[i].rsvr_sample(RSVR_SIZE, in_place=False)
        # To do multiprocessing, I have no choice but turn off in_place option. But I still get large speedup.
    el_time = time.time() - st_time
    print "   Sampled size = %d    Elapsed time = %f" %(len(ret), el_time)
    return ret




def worker_row(i):   # Worker function of row operation in transition
    # Working on row #i as assigned by parent
    st_time = time.time()
    print "Worker is working on transition: row #%d." %(i),
    ret = pv.KWTable()
    row_vec = m_mat[i]
    u       = u_vec[i]

    for j in range(len(row_vec)):               # Vector multiplication
        ret.aggr_inplace(x_vec[j], 1.0, row_vec[j])
    ret.aggr_inplace(y, 1.0, u)                 # Add the observation
    ori_size = len(ret)
    ret = ret.rsvr_sample(RSVR_SIZE, in_place=False)   # Sample tht return
    print "   Size = %d -> %d" %(ori_size, len(ret)),
    el_time = time.time() - st_time
    print "   Elapsed time = %f" %(el_time)
    return ret




if __name__ == "__main__":
    l0 = pv.KWTable()
    b0 = pv.KWTable()
    s0_list = [pv.KWTable() for i in range(PERIOD)]     # create (w-1) empty tables
                                                        # note that s0_list[0] is not used
                                                        # indices 1..(PERIOD-1) are used
    
    print
    print "-" * 80
    print "Using Python interpreter:", sys.executable
    print "N_WORKERS =", N_WORKERS
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
        y = read_rec(TS_CURR)
        b0.aggr_inplace(y, 1.0, -1.0/PERIOD**2)
        print "Read time slot", TS_CURR, "   b0 -= y/W^2",
        TS_CURR += INTERVAL
        el_time = time.time() - st_time
        print "   Elapsed time =", el_time
        
    # ---------- Second training period ----------
    print
    print "---------- Initialization: Second training period ----------"
    for i in range(PERIOD):
        st_time = time.time()
        y = read_rec(TS_CURR)
        l0.aggr_inplace(y, 1.0, 1.0/PERIOD)
        b0.aggr_inplace(y, 1.0, 1.0/PERIOD**2)
        print "Read time slot", TS_CURR, "   l0 += y/W",
        print "   b0 += y/W^2",
        if i < (PERIOD-1):   # Filling S_{0,1}, S_{0,2}, ... , S_{0,W-1} 
            s0_list[i+1].aggr_inplace(y, 1.0, 1.0)
            print "   s0_list[%d] += y" %(i+1),
        TS_CURR += INTERVAL
        el_time = time.time() - st_time
        print "   Elapsed time =", el_time


    for i in range(PERIOD-1):
        st_time = time.time()
        s0_list[i].aggr_inplace(l0, 1.0, -1.0)
        print "s0_list[%d] -= l0" %(i),
        el_time = time.time() - st_time
        print "   Elapsed time =", el_time
 
    # ---------- Form state vector ----------
    print
    print "--------- Form initial state vector and sample each element ----------"

    x_vec = [l0, b0] + s0_list[1:]

    # ---------- Sample each element ----------
    for i in range(len(x_vec)):
        print "x_vec[%d]: KWTable(%-12x)    size = %d" %(i, id(x_vec[i]), len(x_vec[i]))

    wk_pool = mp.Pool(N_WORKERS)
    new_x_vec = wk_pool.map(worker_samp, range(len(x_vec)))
                                                    # Here we do not pass x_vec directly
                                                    # as it will invoke unnecissary 
                                                    # copy operations. Instead, we treat
                                                    # x_vec as global.
    x_vec = new_x_vec   # Since we use not-in-place rsvr_sample, we have to update x_vec

    print
    print "New x_vec is updated. Check the new x_vec."
    for i in range(len(x_vec)):
        print "x_vec[%d]: KWTable(%-12x)    size = %d" %(i, id(x_vec[i]), len(x_vec[i]))
    print
    print "---------- End of initialization ----------"

    # ---------- Recurrence: Transition and forecast ----------
    print
    print "---------- Start of recursion ----------"
    print "Make transition matrix..."
    m_mat = make_trans_matrix()
    print "Make U-vector..."
    u_vec = make_u_vec()

    while TS_CURR <= TS_END:
        # ---------- Forecast based on existing state vector ----------
        print 
        print "## Iteration of time %d ##" %(TS_CURR)
        fcast = forecast(1)     # Make one step forecast.

        # ---------- Transition ----------
        # First read in current time slot's records
        y = read_rec(TS_CURR) 
        print "Read time slot %d, size = %d" %(TS_CURR, len(y))
        
        # Then do the transition of current state vector
        print "Transition of x_vec"
        for i in range(len(x_vec)):
            print "x_vec[%d]: KWTable(%-12x)    size = %d" %(i, id(x_vec[i]), len(x_vec[i]))

        wk_pool = mp.Pool(N_WORKERS)
        new_x_vec = wk_pool.map(worker_row, range(len(x_vec)))
        x_vec = new_x_vec   # Since we use not-in-place rsvr_sample, we have to update x_vec

        print
        print "Transition of x_vec is complete. Check the new x_vec."
        for i in range(len(x_vec)):
            print "x_vec[%d]: KWTable(%-12x)    size = %d" %(i, id(x_vec[i]), len(x_vec[i]))
        print


        TS_CURR += INTERVAL
