#!/usr/bin/python

import pyvaropt as pv
import sys
import os
import time
import multiprocessing as mp

RSVR_SIZE= 30000
DATA_DIR = "/home/users/cgi0911/Data/Waikato_5/hourly_flowbin/"
RES_DIR  = "/home/users/cgi0911/Results/Waikato_5/hourly_flowbin/%s/" %(time.strftime("%Y%m%d-%H%M%S", time.localtime()))
INTERVAL = 3600
TS_CURR  = 1181088000   # Current time stamp
FILETYPE = "flowbin"
PERIOD   = 24
ALPHA    = 0.2
BETA     = 0.2
GAMMA    = 0.2
N_WORKERS= 4

x_vec    = []           # Make it globally accessible

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




def read_rec(ts):
    return pv.KWTable(filetype=FILETYPE, fn=os.path.join(DATA_DIR, str(ts)+".rec"))




def worker_samp(i):
    # Grab a KWTable from queue and sample it down to size
    print "Worker is working on x_vec[%d]: KWTable(%x). Current size = %d" %(i, id(x_vec[i]), len(x_vec[i]))
    ret = x_vec[i].rsvr_sample(RSVR_SIZE, in_place=False)
        # To do multiprocessing, I have no choice but turn off in_place option. But I still get large speedup.
    return ret




if __name__ == "__main__":
    l0 = pv.KWTable()
    b0 = pv.KWTable()
    s0_list = [pv.KWTable() for i in range(PERIOD)]     # create (w-1) empty tables
                                                        # note that s0_list[0] is not used
                                                        # indices 1..(PERIOD-1) are used
    
    print "Using Python interpreter:", sys.executable

    # ---------- First training period ----------
    print "First training period"
    for i in range(PERIOD):
        st_time = time.time()
        y = read_rec(TS_CURR)
        b0.aggr_inplace(y, 1.0, -1.0)
        print "Read time slot", TS_CURR, "   b0 -= y",
        TS_CURR += INTERVAL
        el_time = time.time() - st_time
        print "   Elapsed time =", el_time
        
    # ---------- Second training period ----------
    print
    print "Second training period"
    for i in range(PERIOD):
        st_time = time.time()
        y = read_rec(TS_CURR)
        l0.aggr_inplace(y, 1.0, 1.0)
        b0.aggr_inplace(y, 1.0, 1.0)
        print "Read time slot", TS_CURR, "   l0 += y",
        print "   b0 += y",
        if i > 0:   # (PERIOD-1)..1
            s0_list[PERIOD-i].aggr_inplace(y, 1.0, 1.0)
            print "   s0_list[%d] += y" %(i),
        TS_CURR += INTERVAL
        el_time = time.time() - st_time
        print "   Elapsed time =", el_time


    for i in range(1, PERIOD):
        st_time = time.time()
        s0_list[i].aggr_inplace(l0, 1.0, -1.0)
        print "s0_list[%d] -= l0" %(i),
        el_time = time.time() - st_time
        print "   Elapsed time =", el_time
 
    # ---------- Form state vector ----------
    print
    print "Form initial state vector and VarOpt sample"

    x_vec = [l0, b0] + s0_list[1:]

    for i in range(len(x_vec)):
        print "x_vec[%d]: KWTable(%x)    size = %d" %(i, id(x_vec[i]), len(x_vec[i]))

    """
    for i in range(len(x_vec)):
        st_time = time.time()
        print "Sampling x_vec[%d]   size =" %(i), len(x_vec[i]), "->",
        x_vec[i].rsvr_sample(RSVR_SIZE, in_place=True)
        el_time = time.time() - st_time
        print len(x_vec[i]), "   Elapsed time =", el_time
    """
    
    """
    wk_pool = []

    x_idx_queue = mp.Queue()

    print
    for i in range(len(x_vec)):
        x_idx_queue.put(i)

    for i in range(N_WORKERS):
        wk_name = "WK#%d" %(i)
        p = mp.Process(target=worker_samp, args=(x_idx_queue, wk_name))
        wk_pool.append(p)
        x_idx_queue.put(-1)

    for wk in wk_pool:  wk.start()
    for wk in wk_pool:  wk.join()
    """
    wk_pool = mp.Pool(N_WORKERS)
    new_x_vec = wk_pool.map(worker_samp, range(len(x_vec)))
    x_vec = new_x_vec

    print
    print "Double check!"
    for i in range(len(x_vec)):
        print "new_x_vec[%d]: KWTable(%x)    size = %d" %(i, id(x_vec[i]), len(x_vec[i]))


    # ---------- Start transition and forecast ----------


