#!/usr/bin/python

import sys
import time
import struct
import netaddr as na
import collections as cl
import random as rd
import heapq
import math
import copy
import os
import multiprocessing as mp




class KWTable:
    """The Keyed Weight Table class.
    The main body is a dictionary in which
    - keys are strings (not restricted to a fixed length)
    - values are signed double number.
    I intended to keep the KWTable class general so that it can deal with various types of data formats.
    """
    # ---- Class initializations, variables and parameters ----
    id_counter = 0
    SHOW_DESTRUCTION = False
    SHOW_RSVR_SAMPLE = False
    rd.seed( int(time.time()) )

    def __init__(self, **kwargs):
        """Constructor of a KWTable object.
        It can read from different file formats:
        - "flowbin": 5-tuple records in binary format.
                     Each record is 24 bytes long. First 16 bytes is 5-tuple -
                     (uint32, uint32, uint16, uint16, uint8) padded to 16 bytes, corresponding to
                     (src_ip, dst_ip, src_port, dst_port, protocol).
                     The following 8 bytes is a signed double precision number.
        - "flowtxt": 5-tuple records in text format.
                     Each line is a record. Each record consists of six substrings:
                     src_ip, dst_ip, src_port, dst_port, protocol, weight.
                     Substrings are separated with each other by whitespaces.
                     src_ip and dst_ip are in IPv4 quadrant expressions.
                     port and protocol numbers are general unsigned integer expressed in strings.
                     weight is signed precision number expressed in string.
        - "trace":   Use libtrace. I'm working on it.

        If working 5-tuple-based filetypes:
        - The keys of the table shall be 5-tuples packed into 13-byte strings.
          4 bytes each for src_ip and dst_ip; 2 bytes each for src_port and dst_port; 1 byte for protocol. No padding.
        - The values shall be signed double precision numbers.
        """
        # ---- Member attributes and initializations ----
        self.table  = {}                        # Main body of table
        self.thresh = 0.0                       # Will be assigned value after any rsvr_sampling
        self.posthresh = 0.0                    # Positive component threshold
        self.negthresh = 0.0                    # Negative component threshold

        # ---- Read keyword arguments ----
        if len(kwargs) == 0:    return          # Simply return an empty KWTable
        fn          = kwargs.get("fn")
        filetype    = kwargs.get("filetype")
        in_table    = kwargs.get("in_table")

        # ---- Initialize from file ----
        if fn!=None and filetype!=None:
            if filetype == "flowbin":
                self.read_from_flowbin(fn)
            elif filetype == "flowtxt":
                self.read_from_flowtxt(fn)
            else:
                # traces not yet supported
                sys.stderr.write("Unrecognized file type. Please check your keyword arguments: %s\n" %(str(kwargs)))
                return

        # ---- Initialize from another KWTable ----
        elif in_table!=None:
            self.table = in_table.table.copy()  # Use the built-in copy method of defaultdict


        # ---- Show error message if no way to initialize from file or another KWTable ----
        else:
            sys.stderr.write("Cannot do initilization from file or KWTable. Please check your keyword arguments: %s\n" %(str(kwargs)))
            return




    def __del__(self):
        """Destructor of the KWTable class.
        """
        if KWTable.SHOW_DESTRUCTION == True:
            print "Destroy KWTable(%x)" %(id(self))




    def __repr__(self):
        ret = ""
        count = 0
        for key in self.table:
            if count >= 1000:
                ret += "Printing only 1000 items."
                break
            else:
                src_ip_int, dst_ip_int = struct.unpack("II", key)
                wt = self.table[key]
                ret += "%4d:  %s %s -> %.10f\n" %(count+1, na.IPAddress(src_ip_int), na.IPAddress(dst_ip_int), wt)
                count += 1
        return ret




    def get_sum(self):
        """Simply get the sum of all weights.
        """
        return sum( self.table.values() )




    def get_abssum(self):
        """Get the sum of all absolute weights.
        """
        return sum( [abs(v) for v in self.table.values()] )




    def read_from_flowbin(self, fn):
        """Read and add the information of a flow binary file into the KWTable.
        """
        inFile = open(fn, "rb")
        while True:
            key_str = inFile.read(16)   # Remove the padding zeroes.
            wt_str  = inFile.read(8)
            if len(key_str) != 16:  break
            #key     = key_str[:13]
            key     = key_str[:8]       # Taking only src and dst IP
            wt      = struct.unpack("d", wt_str)[0]
            if not key in self.table:
                self.table[key] = wt
            else:
                self.table[key] += wt

        inFile.close()
        return




    def read_from_flowtxt(self, fn):
        """Read and add the information of a flow text file into the KWTable.
        """
        inFile = open(fn, "r")
        for line in inFile:
            fields = line.rstrip().split(' ')
            src_ip_int = int(na.IPAddress(fields[0]))
            dst_ip_int = int(na.IPAddress(fields[1]))
            #src_port   = int(fields[2])
            #dst_port   = int(fields[3])
            #prot       = int(fields[4])
            wt         = float(fields[5])
            #key        = struct.pack("IIHHB", src_ip_int, dst_ip_int, src_port, dst_port, prot)
            key        = struct.pack("II", src_ip_int, dst_ip_int)  # Taking only src and dst IPs
            if not key in self.table:
                self.table[key] = wt
            else:
                self.table[key] += wt

        inFile.close()
        return




    def to_flowtxt(self, **kwargs):
        """Write the table to assigned buffer (stdout, stderr or file) in flow text format.
        """
        linenum    = kwargs.get("linenum")
        if linenum!=None:   linenum  = int(linenum)
        else:               linenum  = float('inf')

        buf_name    = kwargs.get("buf")
        buf         = None
        if buf_name == None:        buf = sys.stdout
        if buf_name == "stdout":    buf = sys.stdout
        elif buf_name == "stderr":  buf = sys.stderr
        else:                       buf = open(buf_name, "w")

        count = 0
        for key, wt in self.table.iteritems():
            #src_ip_int, dst_ip_int, src_port, dst_port, prot = struct.unpack("IIHHB", key)
            #buf.write("%s %s %d %d %d %.10f\n"  %(na.IPAddress(src_ip_int), na.IPAddress(dst_ip_int), src_port, dst_port, prot, wt))
            src_ip_int, dst_ip_int = struct.unpack("II", key)       # Taking only src and dst IPs
            buf.write("%s %s %.10f\n"  %(na.IPAddress(src_ip_int), na.IPAddress(dst_ip_int), wt))
            count += 1
            if count >= linenum:    break

        if buf!=sys.stdout and buf!=sys.stderr:     # stdout and stderr should not be closed
            buf.close()




    def to_flowbin(self, fn="temp.rec"):
        """Write the table to assigned file in flow binary format.
        """
        outFile = open(fn, "wb")

        for key, wt in self.table.iteritems():
            #rec = struct.pack("13sBBBd", key, 0, 0, 0, wt)
            rec = struct.pack("8s8sd", key, "", wt)
            outFile.write(rec)

        outFile.close()




    def __len__(self):
        return len(self.table)




    def add_item(self, key, wt):
        """Add an item to table. If the key exists, increment the weight by the given wt.
        """
        if not key in self.table:
            self.table[key] = wt
        else:
            self.table[key] += wt

        if self.table[key] == 0.0:  del self.table[key] # Sanity check




    def scaleby(self, coeff):
        """Every weight is scaled by a double precision coefficient.
        """
        ret = KWTable(in_table=self)
        for key in ret.table:
            ret.table[key] *= coeff

        return ret




    def scale_inplace(self, coeff):
        """Do scalar scaling in place.
        """
        for key in self.table:
            self.table[key] *= float(coeff)
        return




    def aggr(self, rhs=None, lcoeff=1.0, rcoeff=1.0):
        """Aggregate two KWTables by given coefficients. Return a new KWTable.
        ret = lcoeff * lhs + rcoeff * rhs
        """
        self.thresh = 0.0               # Reset threshold
        ret = KWTable(in_table=self)


        if lcoeff != 1.0:
            ret.scale_inplace(lcoeff)

        if rhs != None and rcoeff != 0.0:
            for key, wt in rhs.table.iteritems():
                if not key in ret.table:
                    if wt == 0.0:   continue    # Sanity check.
                    ret.table[key] = rcoeff * wt
                else:
                    ret.table[key] += rcoeff * wt
                    if ret.table[key] == 0.0:   del ret.table[key]  # Delete the zero item.

        return ret




    def aggr_inplace(self, rhs=None, lcoeff=1.0, rcoeff=1.0):
        """
        """
        self.thresh = 0.0               # Reset threshold
        for key in self.table:
            self.table[key] *= float(lcoeff)

        for key in rhs.table:
            if key in self.table:   self.table[key] += float(rcoeff) * rhs.table[key]
            else:                   self.table[key] = float(rcoeff) * rhs.table[key]
        
        del_keys = [k for k, v in self.table.iteritems() if v == 0.0]
        for k in del_keys: del self.table[k]




    def rsvr_samp_abswt(self, k):
        """Use reservoir sampling to draw a k-entry sample out of the original KWTable.
        - Shall return a new KWTable.
        - Uses threshold sampling.
        - Uses Horvitz-Thompson estimator to adjust the weights.
        - Please refer to:
          [1] E. Cohen et al., "Composable, Scalable, and Accurate Weight Summarization
              of Unaggregated Data Sets", in VLDB '09.
        """
        # Note that we are working on a KWTable object,
        # of which the table is already uniquely keyed.
        if KWTable.SHOW_RSVR_SAMPLE:    print "VarOpt sample of %s down to size %d." %(hex(id(self)), k)

        if len(self) <= k:              return self # Table size less than or equal to k.
        else:                           return self.rsvr_samp_general(k)




    def rsvr_samp_split(self, k):
        """Split the table into positive and negative components and samp separatedly.
        """
        if KWTable.SHOW_RSVR_SAMPLE:    print "VarOpt sample of %s down to size %d." %(hex(id(self)), k)

        if len(self) <= k:              return self # Table size less than or equal to k.
        else:
            pos = KWTable()
            neg = KWTable()
            # Split table into pos and neg components
            dict_iter = self.table.iteritems()
            for key, val in dict_iter:
                if val > 0:   pos.add_item(key, val)
                if val < 0:   neg.add_item(key, val)
                else:       continue    # Sanity check. Remove 0-weight items.
            
            # Decide rsvr size for pos and neg components
            k_pos = int( k * float(len(pos)) / float(len(pos) + len(neg)) )
            k_neg = k - k_pos
            
            # Sample separately
            pos = pos.rsvr_samp_general(k_pos)
            neg = neg.rsvr_samp_general(k_neg)

            # Aggregate the table
            ret = KWTable()
            ret.aggr_inplace(pos, 1.0, 1.0)
            ret.aggr_inplace(neg, 1.0, 1.0)
            ret.posthresh = pos.thresh
            ret.negthresh = neg.thresh * -1.0
            return ret




    def rsvr_samp_general(self, k):
        """Whether the input is positive, mixed or negative-signed, we use abs(wt) as 
        the sampling weight.
        """
        if len(self) <= k:              return self # Table size less than or equal to k.
        if k == 0:                      return KWTable()    # Empty sanity check

        dict_iter   = self.table.iteritems()    # Dictionary iterator
        l_heap      = []                        # A min-heap for large (weight over threshold) items.
        thresh      = 0.0                       # Threshold
        temp_thresh = 0.0                       # Temporary threshold
        t_list      = []                        # A list for items with abswt = threshold
        x_dict      = {}                        # A dict for pending small items.
        small_sum   = 0.0                       # Sum of all small items.
        rand_num    = 0.0                       # A uniform random number between 0 and 1
        del_keys    = {}                        # Memorize the keys to be deleted (for in-place mode)

        while len(l_heap) < k:                      # Push k items to heap.
            key, wt = dict_iter.next()
            if wt == 0.0:   continue                # Ignore zero items
            heapq.heappush(l_heap, (abs(wt), key))  # Note that each element in the heap is a
                                                    # 2-tuple (abs(wt), key). It is because tuples
                                                    # compared in lexicographical order (compare first, then second)
                                                    # Here we use abs(wt) as our heap weight.

        while True:     # Add an item, drop an item
            try:
                key, wt     = dict_iter.next()
                if wt == 0.0:   continue        # Sanity check. Skip zero-weight item.
                x_dict      = {}
                small_sum   = thresh * len(t_list)

                if abs(wt) > thresh:                # Add new item to the l_heap
                    heapq.heappush(l_heap, (abs(wt), key))
                else:                               # Add new item to the x_dict
                    x_dict[key] = abs(wt)
                    small_sum   += abs(wt)

                while len(l_heap) > 0 and small_sum > (len(t_list) + len(x_dict) -1) * l_heap[0][0]:
                                                            # Make sure to check if l_heap is non-empty...
                    abswt, key  = heapq.heappop(l_heap)     # Move small items and add to x_dict.
                    x_dict[key] = abswt
                    small_sum   += abswt

                temp_thresh = small_sum / (len(t_list) + len(x_dict) - 1)   # Future threshold
                rand_num    = rd.uniform(0, 1)

                if rand_num < len(t_list) * (1 - thresh/temp_thresh):   # Drop an item from t_list
                    d = int(rand_num / (1 - thresh/temp_thresh))
                    dkey = t_list[d][0]
                    t_list[d] = t_list[-1]      # Swap the item to be deleted with list tail.
                                                # Swap-then-delete makes deletion quicker.
                    del t_list[-1]
                else:                                                   # Drop an item from x_dict
                    rand_num -= len(t_list) * (1 - thresh/temp_thresh)
                    x_iter   = x_dict.iteritems()
                    while rand_num > 0:
                        dkey, abswt = x_iter.next()
                        rand_num -= (1 - abswt / temp_thresh)
                    del x_dict[dkey]
                thresh = temp_thresh
                for key, abswt in x_dict.iteritems():   t_list.append((key, abswt))


            except StopIteration:
                break

        # Add items to ret
        ret = KWTable()
        ret.thresh = thresh
        for abswt, key in l_heap:   ret.table[key] = math.copysign(abswt, self.table[key])
        for key, abswt in t_list:   ret.table[key] = math.copysign(thresh, self.table[key])
        return ret




class PrefixQueryTable:
    """The prefix query table. In this table, each row is an IPv4 prefix either in
    the source IP or destination IP. Prefix query table is typically initialized from
    a text file, but can also add/delete items in an ad-hoc manner.
    """
    def __init__(self, fn=None, in_table=None):
        """
        """
        self.table = cl.OrderedDict()                       # OrderedDict will remember item insert order
        self.checkdict = {}

        if fn != None:
            self.from_txt_file(fn)

        elif in_table != None:
            self.table = copy.copy(in_table.table)
            self.checkdict = copy.copy(in_table.checkdict)






    def from_txt_file(self, fn):
        """Parse a text file into prefix query table. Each row is in the format:
        <s or d>,<IP segment in CIDR expression>,<default value>
        An example:
        s,123.45.67.0/24,1000.0

        This function can be repeatedly call to add entries from various files.
        All entries in the table shall be unique.
        """
        inFile = open(fn, "r")
        for line in inFile:
            fields = line.rstrip().split(',')           # Parse the line
            try:
                src_dst     = "d" if fields[0] == "d" else "s"  # src_dst default to s
            except:
                src_dst     = "s"

            try:
                ip_str          = fields[1]
                ip_nw           = na.IPNetwork(ip_str)      # IP network in netaddr format
                ip_prefixlen    = ip_nw.prefixlen           # Prefix length
                ip              = ip_nw.cidr.ip             # CIDR IP in netaddr format
                ip_int          = int(ip)                   # CIDR IP in integer
            except:                                         # Cannot convert the string into IP network
                print "Failed to parse into IP network expression!"
                ip_prefixlen    = 0
                ip_int          = 0

            try:
                wt              = fields[2]
            except:
                wt              = 0.0                       # default to 0.0

            #Add key to self.table and corresponding check item to self.checkdict
            key = (ip_int, ip_prefixlen, src_dst)       # Key is a tuple of (ip in iteger, prefix length, src/dst string)

            if key in self.table:               self.table[key] += 0.0      # Already existing item. Skip processing. Keep counter value.
            else:                               self.table[key] = 0.0       # New item. Add an entry. Counter set to 0.0.

            checkitem = (src_dst, ip_prefixlen)     # A combination of src_dst and prefix length
            if checkitem in self.checkdict:     self.checkdict[checkitem] += 1  # Already existing item.
            else:                               self.checkdict[checkitem] = 1   # New item. Add an entry. Counter set to 1.




    def __del__(self):
        # Use default destructors.
        pass




    def __len__(self):
        return len(self.table)




    def __repr__(self):
        ret = ""
        for key in self.table:
            ip_nw = na.IPNetwork(na.IPAddress(key[0]))
            ip_nw.prefixlen = key[1]
            ip_nw = ip_nw.cidr                          # Remove redundant bits
            src_dst = key[2]
            ret += "%-4s %-20s %-.10f\n" %(src_dst, str(ip_nw), self.table[key])
        return ret



    def wtsum(self, rhs, lcoeff, rcoeff):
        """ The two prefix query tables must have identical key sets!
        """
        if not self.table.keys() == rhs.table.keys() or not self.checkdict.keys() == rhs.checkdict.keys():
            print "Cannot do weighted sum on two prefix query tables. Not perfectly aligned."
            return None

        ret = PrefixQueryTable()
        for key in self.table:      ret.table[key] = self.table[key] * lcoeff + rhs.table[key] * rcoeff
        for key in self.checkdict:  ret.checkdict[key] = self.checkdict[key]

        return ret




    def to_txt(self, fn):
        """
        """
        myFile = open(fn, "w")
        for k in self.table:
            ip_nw = na.IPNetwork(na.IPAddress(k[0]))
            ip_nw.prefixlen = k[1]
            ip_nw = ip_nw.cidr                          # Remove redundant bits
            src_dst = k[2]
            myFile.write("%s,%s,%.10f\n" %(src_dst, str(ip_nw), self.table[k]))
        myFile.close()




    def get_data(self):
        ret = []
        for k in self.table.keys():                     # Key order matters!
            ip_nw = na.IPNetwork(na.IPAddress(k[0]))
            ip_nw.prefixlen = k[1]
            ip_nw = ip_nw.cidr                          # Remove redundant bits
            src_dst = k[2]
            ret.append((src_dst, str(ip_nw), self.table[k]))     # Convert to string
        return ret




    def get_keys(self):
        ret = []
        for k in self.table.keys():                     # Key order matters!
            ip_nw = na.IPNetwork(na.IPAddress(k[0]))
            ip_nw.prefixlen = k[1]
            ip_nw = ip_nw.cidr                          # Remove redundant bits
            src_dst = k[2]
            ret.append(str(ip_nw) + ' ' + src_dst)                      # Convert to string
        return ret




    def get_values(self):
        ret = []
        for k in self.table.keys():                     # Key order matters!
            ip_nw = na.IPNetwork(na.IPAddress(k[0]))
            ip_nw.prefixlen = k[1]
            ip_nw = ip_nw.cidr                          # Remove redundant bits
            ret.append(self.table[k])     # Convert to string
        return ret




    def query(self, kwtable):
        """
        """
        self.reset()
        for kwkey, wt in kwtable.table.iteritems():
            sip, dip = struct.unpack("II", kwkey[:8])
            for checkitem in self.checkdict:
                src_dst = checkitem[0]
                if src_dst == "src":        ipkey = sip
                else:                       ipkey = dip
                prefixlen = checkitem[1]
                ipkey = ipkey >> (32 - prefixlen) << (32 - prefixlen)
                prefixkey = (ipkey, prefixlen, src_dst)

                if prefixkey in self.table:
                    self.table[prefixkey] += wt




    def reset(self):
        """
        """
        for key in self.table:
            self.table[key] = 0.0




if __name__ == "__main__":
    pass
