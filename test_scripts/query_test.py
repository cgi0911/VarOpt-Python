#!/usr/bin/python

import pyvaropt as pv

a = pv.KWTable(filetype="flowbin", fn="./test_data/rec1.rec")
b = pv.KWTable(filetype="flowbin", fn="./test_data/rec2.rec")
c = pv.PrefixQueryTable("./test_data/prefixes.txt")
d = pv.PrefixQueryTable("./test_data/prefixes.txt")
c.query(a)
d.query(b)

e = c.wtsum(d, 1.0, 2.0)

c.to_txt("c.txt")
d.to_txt("d.txt")
e.to_txt("e.txt")
