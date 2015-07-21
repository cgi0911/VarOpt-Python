#!/usr/bin/python

import sys
import pyvaropt as vo
from pyvaropt import KWTable
import time
import argparse as arp
import pprint

parser = arp.ArgumentParser(description="Command line arguments.")
parser.add_argument('-fn1', type=str, default="data/rec1.rec")
parser.add_argument('-fn2', type=str, default="data/rec2.rec")
parser.add_argument('-fn3', type=str, default="data/prefixes.txt")
args = parser.parse_args()
pp = pprint.PrettyPrinter(indent=2)


#@profile
def mainfunc():
    t1      = vo.KWTable(fn=args.fn1, filetype="flowbin")
    t2      = vo.KWTable(fn=args.fn2, filetype="flowbin")
    t3      = t1.aggr(t2, 1.0, -1.0)

    st_time = time.time()
    t4      = t3.rsvr_sample(300000)
    print time.time() - st_time
    st_time = time.time()
    t5      = t1.rsvr_sample(300000)
    print time.time() - st_time
    st_time = time.time()
    t6      = t2.rsvr_sample(300000)
    print time.time() - st_time

    print len(t3), len(t4)

    ptable1 = vo.PrefixQueryTable(fn=args.fn3)
    ptable2 = vo.PrefixQueryTable(in_table=ptable1)
    ptable3 = vo.PrefixQueryTable(in_table=ptable1)
    ptable4 = vo.PrefixQueryTable(in_table=ptable1)

    ptable3.query(t3)
    ptable4.query(t4)

    pp.pprint(ptable3.get_data())
    print
    pp.pprint(ptable4.get_data())




if __name__ == "__main__":
    vo.KWTable.SHOW_DESTRUCTION = True
    mainfunc()
