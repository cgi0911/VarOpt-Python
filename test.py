#!/usr/bin/python

import sys
import varopt as vo
from varopt import KWTable
import time
import argparse as arp

parser = arp.ArgumentParser(description="Command line arguments.")
parser.add_argument('-fn1', type=str, default="rec1.rec")
parser.add_argument('-fn2', type=str, default="rec2.rec")
args = parser.parse_args()


#@profile
def mainfunc():
    KWTable.SHOW_DESTRUCTION = True
    st_time = time.time()
    t1      = vo.KWTable(fn=args.fn1, filetype="flowbin")
    ed_time = time.time()
    print "t1"
    print "Loading file to table. Elapsed time =", ed_time - st_time
    print len(t1)
    print "Sum = %e" %(t1.getsum())
    print "ABS sum = %e" %(t1.getabssum())
    print

    for i in range(5):    
        st_time = time.time()
        t3      = t1.rsvr_sample(300000)
        ed_time = time.time()
        print "t3"
        print "Do reservoir sampling. Elapsed time =", ed_time - st_time
        print len(t3)
        print "Sum = %e" %(t3.getsum())
        print "ABS sum = %e" %(t3.getabssum())
        print




if __name__ == "__main__":
    mainfunc()
