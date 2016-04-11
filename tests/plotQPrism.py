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
  out5 = []
  execfile(file)

  pylab.plot(out2, label='prismLeft', linewidth=1)
  pylab.plot(out5, label='prismRight', linewidth=1)
  pylab.legend()
  pylab.show()
