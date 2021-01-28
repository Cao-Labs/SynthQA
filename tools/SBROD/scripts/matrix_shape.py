#!/usr/bin/python
#sbrod

import sys
from scipy.io import loadmat


if len(sys.argv) != 2:
  print('Usage: {} <matrix path>'.format(sys.argv[0]))
  exit(1)

mat_path = sys.argv[1]

x = loadmat(mat_path)['SpMat_{}'.format(mat_path.split('.')[-2])]
print("shape: {}".format(x.shape))
print("sum: {}".format(x.sum()))
print("nonzero dimensions: {}".format(len(x.nonzero()[0])))
