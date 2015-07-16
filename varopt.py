#!/usr/bin/python

import sys
import time
import struct
import netaddr as na 
import collections as cl
import random as rd
import heapq
import math




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
        self.table  = cl.defaultdict(float)     # Main body of table
        self.tid    = KWTable.id_counter
        KWTable.id_counter  += 1

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
                sys.stderr.write("Unrecognized file type. Please check your keyword arguments: %s" %(str(kwargs)))
                return
        
        # ---- Initialize from another KWTable ----
        elif in_table!=None:
            self.table = in_table.table.copy()


        # ---- Show error message if no way to initialize from file or another KWTable ----
        else:
            sys.stderr.write("Cannot do initilization from file or KWTable. Please check your keyword arguments: %s" %(str(kwargs)))
            return




    def __del__(self):
        """Destructor of the KWTable class.
        """
        if KWTable.SHOW_DESTRUCTION == True:
            print "Deconstructing KWTable #%d at address: %s" %(self.tid, hex(id(self)))



    
    def getsum(self):
        """Simply get the sum of all weights.
        """
        return sum( self.table.values() )




    def getabssum(self):
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
            key     = key_str[:13]
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
            src_ip_int = int(na.IPAddress(field[0]))
            dst_ip_int = int(na.IPAddress(field[1]))
            src_port   = int(field[2])
            dst_port   = int(field[3])
            prot       = int(field[4])
            wt         = float(field[5])
            key        = struct.pack("IIHHB", src_ip_int, dst_ip_int, src_port, dst_port, prot)
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
            src_ip_int, dst_ip_int, src_port, dst_port, prot = struct.unpack("IIHHB", key)
            buf.write("%s %s %d %d %d %.10f\n"  %(na.IPAddress(src_ip_int), na.IPAddress(dst_ip_int), src_port, dst_port, prot, wt))
            count += 1
            if count >= linenum:    break
        
        if buf!=sys.stdout and buf!=sys.stderr:     # stdout and stderr should not be closed
            buf.close()




    def to_flowbin(self, fn="temp.rec"):
        """Write the table to assigned file in flow binary format.
        """
        outFile = open(fn, "wb")

        for key, wt in self.table.iteritems():
            rec = struct.pack("13sBBBd", key, 0, 0, 0, wt)
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




    def __add__(lhs, rhs):
        """Aggregate two KWTables. If there are overlapped keys in both KWTables, add their weight.
        Shall return a new KWTable.
        """
        ret = KWTable()
        for key, wt in lhs.table.iteritems():
            ret.add_item(key, wt)
        for key, wt in rhs.table.iteritems():
            ret.add_item(key, wt)

        return ret




    def __iadd__(lhs, rhs):
        """Aggregate rhs to lhs in place.
        """
        for key, wt in rhs.table.iteritems():
            #lhs.table[key] = lhs.table.get(key, 0.0) + wt
            if not key in lhs.table:
                lhs.table[key] = wt
            else:
                lhs.table[key] += wt
        return lhs




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
            self.table[key] *= coeff

        return
        



    @classmethod
    def aggr(cls, lhs, rhs, lcoeff, rcoeff):
        """Aggregate two KWTables by given coefficients. Return a new KWTable.
        ret = lcoeff * lhs + rcoeff * rhs
        """
        ret = KWTable(in_table=lhs)
        ret.scale_inplace(lcoeff)
        for key, wt in rhs.table.iteritems():
            if not key in ret.table:
                ret.table[key] = rcoeff * wt
            else:
                ret.table[key] += rcoeff * wt

        return ret



    #@profile
    def rsvr_sample(self, k):
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

        if len(self) <= k:  # Table size less than or equal to k.
            return KWTable(in_table=self)

        else:
            ret         = KWTable()
            dict_iter   = self.table.iteritems()    # Dictionary iterator
            l_heap      = []                        # A min-heap for large (weight over threshold) items.
            thresh      = 0.0                       # Threshold
            temp_thresh = 0.0                       # Temporary threshold
            t_list      = []                        # A list for items with abswt = threshold
            x_dict      = cl.defaultdict(float)     # A dict for interim items.
            small_sum   = 0.0                       # Sum of all small items.
            rand_num    = 0.0                       # A uniform random number between 0 and 1


            for i in range(k):                      # Push k items to heap.
                key, wt = dict_iter.next()
                heapq.heappush(l_heap, (abs(wt), key))  # Note that each element in the heap is a 
                                                        # 2-tuple (abs(wt), key). It is because tuples
                                                        # compared in lexicographical order (compare first, then second)
                                                        # Here we use abs(wt) as our heap weight.

            while True:     # Add an item, drop an item
                try:
                    key, wt     = dict_iter.next()
                    x_dict      = cl.defaultdict(float)
                    small_sum   = thresh * len(t_list)

                    if abs(wt) > thresh:    heapq.heappush(l_heap, (abs(wt), key))
                    else:
                        x_dict[key] = abs(wt)
                        small_sum   += abs(wt)

                    while small_sum >= (len(t_list) + len(x_dict) -1) * l_heap[0][0]:
                        abswt, key  = heapq.heappop(l_heap)
                        x_dict[key] = abswt
                        small_sum   += abswt

                    temp_thresh = small_sum / (len(t_list) + len(x_dict) - 1)
                    rand_num    = rd.uniform(0, 1)
                    
                    if rand_num < len(t_list) * (1 - thresh/temp_thresh):   # Drop an item from t_list
                        d = int(rand_num / (1 - thresh/temp_thresh))
                        t_list[d] = t_list[-1]
                        del t_list[-1]
                    else:                                                   # Drop an item from x_list
                        rand_num -= len(t_list) * (1 - thresh/temp_thresh)
                        x_iter   = x_dict.iteritems() 
                        while rand_num > 0:
                            key, abswt = x_iter.next()
                            rand_num -= (1 - abswt / temp_thresh)
                        del x_dict[key]

                    thresh = temp_thresh
                    for key, abswt in x_dict.iteritems():   t_list.append((key, abswt))


                except StopIteration:
                    break



            
            for abswt, key in l_heap:
                ret.table[key] = math.copysign(abswt, self.table[key])

            for key, abswt in t_list:
                ret.table[key] = math.copysign(thresh, self.table[key])

            return ret
                        



if __name__ == "__main__":
    t1      = KWTable(fn="rec3.rec", filetype="flowbin")
    t2      = t1.rsvr_sample(10000)
