#!/usr/bin/python

import time
import os
from pyvaropt import KWTable

if __name__ == "__main__":
    tot_size = 0
    fns = sorted(os.listdir("/home/cgi0911/sg3000/Waikato_V_traces/hourly_flowbin/"))
    tables = []

    st = time.time()
    for i in range (20):
        t = KWTable(filetype="flowbin", fn="/home/cgi0911/sg3000/Waikato_V_traces/hourly_flowbin/" + fns[i])
        tables.append(t)
        tot_size += len(t)

    print "Reading time =", (time.time() - st)
    print "Original size =", tot_size

    ret = KWTable()

    st = time.time()
    for t in tables:
        ret.aggr_inplace(rhs=t)

    print "Aggregation time =", (time.time() - st)
    print "Size after aggregation =", len(ret)

    st = time.time()
    ret.rsvr_sample(100000)
    print "Sample time =", (time.time() - st)
    print "Size after sampling =", len(ret)
