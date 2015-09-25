#!/usr/bin/python

import pyvaropt as pv
import os

DATA_DIR    = "/home/users/cgi0911/Results/Waikato_5/fcast"
FN_PREF     = "./test_data/prefix_all.txt"
QUERY_DIR   = "/home/users/cgi0911/Results/Waikato_5/temp/queried"

if not os.path.exists(QUERY_DIR):   os.makedirs(QUERY_DIR)

data_fns = sorted(os.listdir(DATA_DIR))
q_table  = pv.PrefixQueryTable(fn=FN_PREF)

for data_fn in data_fns:
    print data_fn
    query_fn = data_fn.replace(".rec", ".txt")
    kw_table = pv.KWTable(filetype="flowbin", fn=os.path.join(DATA_DIR, data_fn))
    q_table.query(kw_table)
    q_table.to_txt(os.path.join(QUERY_DIR, query_fn))
    q_table.reset()
