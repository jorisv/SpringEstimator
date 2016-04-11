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

set(SOURCES SpringEstimator.cpp)
set(HEADERS SpringEstimator.h)

add_library(SpringEstimator SHARED ${SOURCES} ${HEADERS})
PKG_CONFIG_USE_DEPENDENCY(SpringEstimator SpaceVecAlg)
PKG_CONFIG_USE_DEPENDENCY(SpringEstimator RBDyn)

set(INSTALL_PATH include/SpringEstimator)
install(FILES ${HEADERS} DESTINATION ${INSTALL_PATH})
install(TARGETS SpringEstimator DESTINATION "lib")
