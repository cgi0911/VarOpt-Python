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
PERIOD   = 24    # # of time slots in a period
R        = 1    # Forecast # of steps
ALPHA    = 0.2
BETA     = 0.2
GAMMA    = 0.2
N_WORKERS= 4

# ---------- Global variables and objects ----------
TS_CURR  = 0                # Current timestamp
x_vec    = mp.Manager().dict()  # State vector. Make it globally accessible (shared dictionary)
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
    print "Making %d-step forecast. Time stamp = %d. Elapsed time = %f" %(R, TS_CURR + INTERVAL * (R-1), el_time)
    return ret




def worker_samp(wkid, task_queue, res_dict):
    while True:
        st_time = time.time()
        i = task_queue.get()
        if i == -1:     break       # Break when stop token received
        print "Worker #%d is working on x_vec[%d]: KWTable(%-12x). Current size = %d" %(wkid, i, id(x_vec[i]), len(x_vec[i])),
        res = x_vec[i].rsvr_sample(RSVR_SIZE, in_place=False)
        el_time = time.time() - st_time
        print "   Sampled size = %d    Elapsed time = %f" %(len(res), el_time)
        res_dict[i] = res
    return




def worker_row(wkid, task_queue, res_dict):
    while True:
        st_time = time.time()

        i = task_queue.get()
        if i == -1:     break       # Break when stop token received
        print "Worker #%d is working on transition of row #%d." %(wkid, i),

        res = pv.KWTable()

        row_vec = m_mat[i]
        u       = u_vec[i]

        for j in range(len(row_vec)):
            if row_vec[j] == 0.0:   continue    # No need to do 0-coeff aggregation
            res.aggr_inplace(x_vec[j], 1.0, row_vec[j])

        if not u == 0.0:    res.aggr_inplace(y, 1.0, u)

        ori_size = len(res)
        res.rsvr_sample(RSVR_SIZE, in_place=True)
        print "   Size = %d -> %d" %(ori_size, len(res)),
        el_time = time.time() - st_time
        print "   Elapsed time = %f" %(el_time)
        res_dict[i] = res
    return


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

    task_queue = mp.Queue()
    res_dict = mp.Manager().dict()
    wk_pool = []

    for i in range(len(x_vec)):     task_queue.put(i)
    for i in range(N_WORKERS):      task_queue.put(-1)

    for i in range(N_WORKERS):
        p = mp.Process(target=worker_samp, args=(i, task_queue, res_dict))
        wk_pool.append(p)
        p.start()

    for p in wk_pool:               p.join()

    for i in range(len(x_vec)):
        x_vec[i] = res_dict[i]

    del l0
    del b0
    del s0_list         # Remove unused KWTables to release memory

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
        st_time_recur = time.time()
        # ---------- Forecast based on existing state vector ----------
        print
        print "## Iteration of time %d ##" %(TS_CURR)
        fcast = forecast(1)     # Make one step forecast.
        print "Forcast abs_sum =", fcast.get_abssum()

        # ---------- Transition ----------
        # First read in current time slot's records
        y = read_rec(TS_CURR)
        print "Read time slot %d, size = %d" %(TS_CURR, len(y))

        # Then do the transition of current state vector
        print "Transition of x_vec..."
        #for i in range(len(x_vec)):
        #    print "x_vec[%d]: KWTable(%-12x)    size = %d" %(i, id(x_vec[i]), len(x_vec[i]))

        task_queue = mp.Queue()
        res_dict = mp.Manager().dict()
        wk_pool = []

        for i in range(len(x_vec)):     task_queue.put(i)
        for i in range(N_WORKERS):      task_queue.put(-1)

        for i in range(N_WORKERS):
            p = mp.Process(target=worker_row, args=(i, task_queue, res_dict))
            wk_pool.append(p)
            p.start()

        for p in wk_pool:               p.join()

        for i in range(len(x_vec)):     x_vec[i] = res_dict[i]

        el_time_recur = time.time() - st_time_recur
        print "Transition of x_vec is complete. Elapsed time = %f" %(el_time_recur)
        for i in range(len(x_vec)):
            print "x_vec[%d]: KWTable(%-12x)    size = %d" %(i, id(x_vec[i]), len(x_vec[i]))
        print


        TS_CURR += INTERVAL
