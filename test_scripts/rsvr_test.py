#!/usr/bin/python

import pyvaropt as pv

print "Loading a"
a = pv.KWTable(filetype="flowbin", fn="/home/users/cgi0911/Data/Waikato_5/hourly_flowbin/1184932800.rec")
print "Loading b"
b = pv.KWTable(filetype="flowbin", fn="/home/users/cgi0911/Data/Waikato_5/hourly_flowbin/1184936400.rec")

asum = a.get_sum()
bsum = b.get_sum()
print "asum = %e" %(asum)
print "bsum = %e" %(bsum)

a_coeff = 0.5
b_coeff = -0.5

c = a.aggr(b, a_coeff, b_coeff)
csum = c.get_sum()
cabssum = c.get_abssum()
print "c = %f * a + %f * b" %(a_coeff, b_coeff)
print "%f * asum + %f * bsum = %e" %(a_coeff, b_coeff, a_coeff*asum+b_coeff*bsum)
print "csum = %e" %(csum)
print "cabssum = %e" %(cabssum)

print
print "d = c.rsvr_sample(30000, in_place=False)"
d = c.rsvr_sample(30000, in_place=False)
dsum = d.get_sum()
dabssum = d.get_abssum()
print "dsum = %e" %(dsum)
print "dabssum = %e" %(dabssum)
print "len(c) = %d" %(len(c))
print "len(d) = %d" %(len(d))

print
print "c.rsvr_sample(10000, in_place=True)"
c.rsvr_sample(10000, in_place=True)
csum = c.get_sum()
cabssum = c.get_abssum()
print "csum = %e" %(csum)
print "cabssum = %e" %(cabssum)
print "len(c) = %d" %(len(c))
