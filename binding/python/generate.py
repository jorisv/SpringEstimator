# Copyright 2015-2016 CNRS-UM LIRMM, CNRS-AIST JRL
#
# This file is part of SpringEstimator.
#
# SpringEstimator is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SpringEstimator is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with SpringEstimator.  If not, see <http://www.gnu.org/licenses/>.

from pybindgen import *
import sys



def import_eigen3_types(mod):
  mod.add_class('Matrix3d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')
  mod.add_class('VectorXd', foreign_cpp_namespace='Eigen', import_from_module='eigen3')



def import_rbd_types(mod):
  mod.add_class('MultiBody', foreign_cpp_namespace='rbd', import_from_module='rbdyn')



if __name__ == '__main__':
  if len(sys.argv) < 2:
    sys.exit(1)

  se = Module('_spring_estimator', cpp_namespace='::spring_estimator')
  se.add_include('<SpringEstimator.h>')

  dom_ex = se.add_exception('std::domain_error', foreign_cpp_namespace=' ',
                            message_rvalue='%(EXC)s.what()')
  run_ex = se.add_exception('std::runtime_error', foreign_cpp_namespace=' ',
                            message_rvalue='%(EXC)s.what()')
  out_ex = se.add_exception('std::out_of_range', foreign_cpp_namespace=' ',
                            message_rvalue='%(EXC)s.what()')

  # import Eigen3 and sva types
  import_eigen3_types(se)
  import_rbd_types(se)

  # build list type
  se.add_container('std::vector<rbd::MultiBody>', 'rbd::MultiBody', 'vector')

  # build JointTarget structur
  jointTarget = se.add_struct('JointTarget')
  jointTarget.add_constructor([])
  jointTarget.add_constructor([param('int', 'armIndex'),
                               param('const std::string&', 'jointName'),
                               param('double', 'target')])
  jointTarget.add_copy_constructor()
  jointTarget.add_instance_attribute('armIndex', 'int')
  jointTarget.add_instance_attribute('jointName', 'std::string')
  jointTarget.add_instance_attribute('target', 'double')

  se.add_container('std::vector<spring_estimator::JointTarget>',
                   'spring_estimator::JointTarget', 'vector')

  springEst = se.add_class('SpringEstimator')
  springEst.add_constructor([])
  springEst.add_copy_constructor()

  springEst.add_method('initArms', None,
                       [param('const std::vector<rbd::MultiBody>&', 'arms'),
                        param('const std::vector<spring_estimator::JointTarget>&',
                              'jointsTarget')])
  springEst.add_method('updateArms', None,
                       [param('const std::vector<rbd::MultiBody>&', 'arms')])

  springEst.add_method('target', None,
                       [param('const Eigen::Matrix3d&', 'target')])
  springEst.add_method('target', retval('Eigen::Matrix3d'),
                        [], is_const=True)

  springEst.add_method('q', None,
                       [param('const Eigen::VectorXd&', 'q')])
  springEst.add_method('q', retval('Eigen::VectorXd'),
                       [], is_const=True)

  springEst.add_method('qd', retval('Eigen::VectorXd'),
                       [], is_const=True)

  springEst.add_method('taskError', retval('Eigen::VectorXd'),
                       [param('int', 'taskIndex')], is_const=True)

  def doubleSetter(name):
    springEst.add_method(name, retval('double'), [], is_const=True)
    springEst.add_method(name, None, [param('double', 'threshold')])

  map(doubleSetter, ['leastSquareRelTol', 'leastSquareMinTol',
                     'projectorRelTol', 'projectorMinTol'])

  springEst.add_method('update', retval('double'),
                       [param('double', 'timeStep'),
                        param('int', 'nrIter')])


  with open(sys.argv[1], 'w') as f:
    se.generate(f)

