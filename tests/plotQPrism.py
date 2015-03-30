#! /usr/bin/env python
# -*- coding: utf-8 -*-
#

import sys

import pylab

if __name__ == '__main__':
  file = sys.argv[1]
  out0 = []
  out1 = []
  out2 = []
  out3 = []
  out4 = []
  out5 = []
  execfile(file)

  pylab.plot(out2, label='prismLeft', linewidth=1)
  pylab.plot(out5, label='prismRight', linewidth=1)
  pylab.legend()
  pylab.show()
