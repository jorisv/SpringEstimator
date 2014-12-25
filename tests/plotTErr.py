#! /usr/bin/env python
# -*- coding: utf-8 -*-
#

import sys

import pylab

if __name__ == '__main__':
  file = sys.argv[1]
  t1Err = []
  t2Err = []
  t3Err = []
  t4Err = []
  execfile(file)

  pylab.plot(t1Err, label='positon error')
  pylab.plot(t2Err, label='orientation error')
  pylab.plot(t3Err, label='target error')
  pylab.plot(t4Err, label='velocity error')
  pylab.legend()
  pylab.show()
