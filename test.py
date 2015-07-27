#!/usr/bin/python

import sys
import pyvaropt as pvo
from pyvaropt import KWTable
import time
import argparse as arp
import pprint

parser = arp.ArgumentParser(description="Command line arguments.")
parser.add_argument('-fn1', type=str, default="/home/cgi0911/sg3000/Waikato_V_traces/hourly_flowbin/")
#parser.add_argument('-fn2', type=str, default="data/rec2.rec")
#parser.add_argument('-fn3', type=str, default="data/prefixes.txt")
args = parser.parse_args()
pp = pprint.PrettyPrinter(indent=2)


#@profile
def mainfunc():
    myPredict = pvo.Prediction(folder=args.fn1, k=300000)
    print myPredict.folder
    print myPredict.fn_list




if __name__ == "__main__":
    pvo.KWTable.SHOW_DESTRUCTION = True
    mainfunc()
