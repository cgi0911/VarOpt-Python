#!/usr/bin/python

import pyvaropt as pv
import sys
import os
import time

RSVR_SIZE= 30000

#DATA_HOME= "/home/users/cgi0911/Data/Waikato_5/"
DATA_HOME= "/home/cgi0911/sg3000/PacketTraces/Waikato_5/"
DATA_DIR = os.path.join(DATA_HOME, "hourly_flowbin")

#RES_HOME = "/home/users/cgi0911/Results/Waikato_5/"
RES_HOME = "/home/cgi0911/Results/Waikato_5/"
RES_DIR  = os.path.join(RES_HOME, "temp/%s/" %(time.strftime("%Y%m%d-%H%M%S", time.localtime())))

INTERVAL = 3600         # Seconds in a time slot
TS_START = 1181088000   # Starting timestamp (in seconds)
TS_END   = TS_START + INTERVAL * 16 * 24    # Ending timestamp
FILETYPE = "flowbin"
PERIOD   = 24   # # of time slots in a period
M        = PERIOD + 1   # Dimension of transition matrix
R        = 1    # Forecast # of steps
ALPHA    = 0.800
BETA     = 0.100
GAMMA    = 0.200

SAMPLE_MODE = "abswt"

# ---------- Global variables and objects ----------
TS_CURR  = 0                # Current timestamp
x_vec    = []               # State vector. Make it globally accessible (shared dictionary)
y        = pv.KWTable()     # Current time slot's observation
m_mat    = []               # Transition matrix
u_vec    = []               # U-Vector for transition

RES_DIR_FCAST = os.path.join(RES_DIR, "fcast")
if not os.path.exists(RES_DIR_FCAST):   os.makedirs(RES_DIR_FCAST)




def make_trans_matrix():
    ret = []
    for i in range(M):
        ret.append([0.0] * (M))

    ret[0][0] = 1.0 - ALPHA
    ret[0][1] = 1.0 - ALPHA
    ret[0][2] = -1.0 * ALPHA
    ret[1][0] = -1.0 * BETA
    ret[1][1] = 1.0 - BETA
    ret[1][2] = -1.0 * BETA

    for i in range(2, M-1):
        ret[i][0]   = GAMMA/PERIOD
        ret[i][1]   = GAMMA/PERIOD
        ret[i][2]   = GAMMA/PERIOD
        ret[i][i+1] = 1.0

    ret[M-1][0] = GAMMA/PERIOD
    ret[M-1][1] = GAMMA/PERIOD
    ret[M-1][2] = GAMMA/PERIOD - 1.0
    for i in range(3,M):    ret[M-1][i] = -1.0

    return ret




def make_u_vec():
    ret = [0.0] * M
    ret[0] = ALPHA
    ret[1] = BETA
    for i in range(2, M):   ret[i] = -1.0 * GAMMA / PERIOD
    return ret




def read_rec(ts):
    return pv.KWTable(filetype=FILETYPE, fn=os.path.join(DATA_DIR, str(ts)+".rec"))




def forecast(r):
    # r is the # of steps to look ahead.
    if r < 1:
        print "Cannot make forecast! Wrong number of steps = %d." %(r)
        return None

    st_time = time.time()
    ret = pv.KWTable()
    l = x_vec[0]
    b = x_vec[1]
    s = x_vec[r % PERIOD + 1]  # S_{t-w+1}
    ret.aggr_inplace(l, 1.0, 1.0)
    ret.aggr_inplace(b, 1.0, float(r))
    ret.aggr_inplace(s, 1.0, 1.0)
    if SAMPLE_MODE == "abswt":
        ret = ret.rsvr_samp_abswt(RSVR_SIZE)
    else:
        ret = ret.rsvr_samp_split(RSVR_SIZE)
    el_time = time.time() - st_time
    print "Making %d-step forecast. Time stamp = %d. Elapsed time = %f" %(R, TS_CURR + INTERVAL * (R-1), el_time)
    return ret




def samp_x_vec():
    for i in range(len(x_vec)):
        st_time = time.time()
        print "Sampling x_vec[%d]: KWTable(%-12x). Current size = %d" %(i, id(x_vec[i]), len(x_vec[i])),
        if SAMPLE_MODE == "abswt":
            res = x_vec[i].rsvr_samp_abswt(RSVR_SIZE)
        else:
            res = x_vec[i].rsvr_samp_split(RSVR_SIZE)
        el_time = time.time() - st_time
        print "   Sampled size = %d    Elapsed time = %f" %(len(res), el_time)
        x_vec[i] = res
    return




def transition():
    ret = [pv.KWTable() for _ in range(M)]              # Must first return a new x_vec, then overwrite the x_vec
                                                        # Every element must be initialized individually!!
    
    print "SUM of x_vec:",
    for i in range(len(x_vec)):     print "%e" %(x_vec[i].get_sum()),
    print

    print "ABSSUM of x_vec:",
    for i in range(len(x_vec)):     print "%e" %(x_vec[i].get_abssum()),
    print



    for i in range(M):
        row_vec = m_mat[i]
        u       = u_vec[i]

        for j in range(len(row_vec)):
            if row_vec[j] == 0.0:   continue    # No need to do 0-coeff aggregation
            print "%.2e * %.4f +" %(x_vec[j].get_sum(), row_vec[j]),
            ret[i].aggr_inplace(x_vec[j], 1.0, row_vec[j])

        print "(%.2e * %.4f)" %(y.get_sum(), u),
        if not u == 0.0:    ret[i].aggr_inplace(y, 1.0, u)
        print " -> sum = %.2e, abssum = %.2e," %(ret[i].get_sum(), ret[i].get_abssum()),
        if SAMPLE_MODE == "abswt":
            ret[i] = ret[i].rsvr_samp_abswt(RSVR_SIZE)
        else:
            ret[i] = ret[i].rsvr_samp_split(RSVR_SIZE)
        print "thr = %.2e, posthr = %.2e, negthr = %.2e" %(ret[i].thresh, ret[i].posthresh, ret[i].negthresh)

    return ret




if __name__ == "__main__":
    l0 = pv.KWTable()
    b0 = pv.KWTable()
    s0_list = []
    #s0_list = [pv.KWTable() for _ in range(PERIOD)]     # create (w-1) empty tables
                                                        # note that s0_list[0] is not used
                                                        # indices 1..(PERIOD-1) are used

    print
    print "-" * 80
    print "Using Python interpreter:", sys.executable
    print "Data source folder:", DATA_DIR
    print "Result folder:", RES_DIR
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
    st_time_init = time.time()
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
            #s0_list[i+1].aggr_inplace(y, 1.0, 1.0)
            s0_list.append(pv.KWTable())
            s0_list[-1].aggr_inplace(y, 1.0, 1.0) 
            print "   s0_list[%d] += y" %(i),
        TS_CURR += INTERVAL
        el_time = time.time() - st_time
        print "   Elapsed time =", el_time


    for i in range(len(s0_list)):     # 0...(PERIOD-1)
        st_time = time.time()
        s0_list[i].aggr_inplace(l0, 1.0, -1.0)
        print "s0_list[%d] -= l0" %(i),
        el_time = time.time() - st_time
        print "   Elapsed time =", el_time

    # ---------- Form state vector ----------
    print
    print "--------- Form initial state vector and sample each element ----------"
    x_vec = [l0, b0] + s0_list


    # ---------- Sample each element ----------
    for i in range(len(x_vec)):
        print "x_vec[%d]: KWTable(%-12x)    size = %d" %(i, id(x_vec[i]), len(x_vec[i]))

    samp_x_vec()   # Sample each element in x_vec

    del l0, b0, s0_list         # Remove unused KWTables to release memory

    print
    print "New x_vec is updated. Check the new x_vec."
    for i in range(len(x_vec)):
        print "x_vec[%d]: KWTable(%-12x)    size = %d" %(i, id(x_vec[i]), len(x_vec[i]))
    print
    el_time_init = time.time() - st_time_init
    print "---------- End of initialization. Elapsed time = %f ----------" %(el_time_init)

    # ---------- Recurrence: Transition and forecast ----------
    print
    print "---------- Start of recursion ----------"
    print "Make transition matrix..."
    m_mat = make_trans_matrix()
    print "Make U-vector..."
    u_vec = make_u_vec()

    while TS_CURR < TS_END:
        st_time_recur = time.time()
        # ---------- Forecast based on existing state vector ----------
        print
        print "## Iteration of time %d ##" %(TS_CURR)
        fcast = forecast(1)     # Make one step forecast.
        fcast_fn = os.path.join(RES_DIR_FCAST, "%d.rec" %(TS_CURR + (R-1)*INTERVAL))
        fcast.to_flowbin(fn=fcast_fn)
        print "fcast abssum = %e" %(fcast.get_abssum())

        # ---------- Transition ----------
        # First read in current time slot's records
        y = read_rec(TS_CURR)
        print "Read time slot %d, size = %d" %(TS_CURR, len(y))

        # Then do the transition of current state vector
        print "Transition of x_vec..."
        new_x_vec = transition()
        x_vec = new_x_vec

        el_time_recur = time.time() - st_time_recur
        print "Transition of x_vec is complete. Elapsed time = %f" %(el_time_recur)

        TS_CURR += INTERVAL
