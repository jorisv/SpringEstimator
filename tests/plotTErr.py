#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2015-2016 CNRS-UM LIRMM, CNRS-AIST JRL
#
# This file is part of SpringEstimator.
#
# Tasks is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Tasks is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Tasks.  If not, see <http://www.gnu.org/licenses/>.

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
