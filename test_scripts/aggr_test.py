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

a_coeff = 0.2
b_coeff = 0.8

c = a.aggr(b, a_coeff, b_coeff)
csum = c.get_sum()
print "c = %f * a + %f * b" %(a_coeff, b_coeff)
print "%f * asum + %f * bsum = %e" %(a_coeff, b_coeff, a_coeff*asum+b_coeff*bsum)
print "csum = %e" %(csum)

print
print "a aggregate in place with b."
a.aggr_inplace(b, a_coeff, b_coeff)
print "a = %f * a + %f * b" %(a_coeff, b_coeff)
print "%f * asum + %f * bsum = %e" %(a_coeff, b_coeff, a_coeff*asum+b_coeff*bsum)
aasum = a.get_sum()
print "aasum = %e" %(aasum)
