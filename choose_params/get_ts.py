#!/usr/bin/python

import pyvaropt as pv
import os

MYDIR = "/home/users/cgi0911/Data/Waikato_5/hourly_flowbin"

fns = os.listdir(MYDIR)

myfile = open("./ts.txt", "w")

for fn in fns:
    fn_path = os.path.join(MYDIR, fn)
    mytable = pv.KWTable(filetype="flowbin", fn=fn_path)
    mysum = mytable.get_sum()
    myfile.write("%f\n" %(mysum))

myfile.close()
