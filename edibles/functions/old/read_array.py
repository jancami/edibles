import numpy as np

def read_array(filename, dtype, separator=','):
     """ Read a file with an arbitrary number of columns.
         The type of data in each column is arbitrary
         It will be cast to the given dtype at runtime
     """
     cast = np.cast
     data = [[] for dummy in xrange(len(dtype))]
     for line in open(filename, 'r'):
         fields = line.strip().split(separator)
         for i, number in enumerate(fields):
             data[i].append(number)
     for i in xrange(len(dtype)):
         data[i] = cast[dtype[i]](data[i])
     return np.rec.array(data, dtype=dtype)

