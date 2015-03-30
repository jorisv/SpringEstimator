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
  execfile(file)



  pylab.plot(out0, label='positon error', linewidth=1)
  pylab.plot(out1, label='orientation error', linewidth=1)
  pylab.plot(out2, label='target error', linewidth=1)
  # out4 is only fill if there is the prismatic error set
  if len(out4) > 0:
    pylab.plot(out3, label='prism error', linewidth=1)
    pylab.plot(out4, label='velocity error', linewidth=1)
  else:
    pylab.plot(out3, label='velocity error', linewidth=1)
  pylab.legend()
  pylab.show()
