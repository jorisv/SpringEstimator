// This file is part of SpringEstimator.
//
// SpringEstimator is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// SpringEstimator is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with SpringEstimator.  If not, see <http://www.gnu.org/licenses/>.

// check memory allocation in some method
#define EIGEN_RUNTIME_NO_MALLOC

// includes
// std
#include <array>
#include <fstream>
#include <iostream>

// boost
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SpringEstimatorTest
#include <boost/test/unit_test.hpp>

// RBDyn
#include <RBDyn/FK.h>

// SpringEstimator
#include "SpringEstimator.h"

// Arms
#include "TwoDofArm.h"
#include "ThreeDofArm.h"


template<int ArraySize>
void toPython(const std::vector<std::array<double, ArraySize>>& data,
              std::size_t dim, const std::string& name, std::ofstream& out)
{
  out << name << "= [";
  for(std::size_t i = 0; i < data.size(); ++i)
  {
    out << data[i][dim] << ",";
  }
  out << "]" << std::endl;
}


template<int ArraySize>
void toPython(const std::vector<std::array<double, ArraySize>>& data,
              const std::string& filename)
{
  std::ofstream out(filename);
  for(int i = 0; i < ArraySize; ++i)
  {
    std::stringstream str;
    str << "out" << i;
    toPython<ArraySize>(data, i, str.str(), out);
  }
}


void testEndEffectors(const rbd::MultiBodyConfig& mbc1, const rbd::MultiBodyConfig& mbc2,
  int endEffectorIndex, const Eigen::Matrix3d& target, double tol=1e-6)
{
  BOOST_CHECK_SMALL((mbc1.bodyPosW[endEffectorIndex].translation() -
                     mbc2.bodyPosW[endEffectorIndex].translation()).norm(), tol);
  BOOST_CHECK_SMALL(sva::rotationError(mbc1.bodyPosW[endEffectorIndex].rotation(),
                    mbc2.bodyPosW[endEffectorIndex].rotation(), 1e-7).norm(), tol);
  BOOST_CHECK_SMALL(sva::rotationError(mbc1.bodyPosW[endEffectorIndex].rotation(),
                    target, 1e-7).norm(), tol);
  BOOST_CHECK_SMALL(sva::rotationError(mbc2.bodyPosW[endEffectorIndex].rotation(),
                    target, 1e-7).norm(), tol);
}


BOOST_AUTO_TEST_CASE(SpringEstimatorTest)
{
  using namespace Eigen;
  using namespace sva;
  using namespace rbd;

  rbd::MultiBody mb;
  rbd::MultiBodyConfig mbc;
  rbd::MultiBodyGraph mbg;

  std::tie(mb, mbc, mbg) = makeTwoDofArm();

  PTransformd X_0_l(Vector3d(-0.1, 0., 0.));
  PTransformd X_0_r(Vector3d(0.1, 0., 0.));

  MultiBody mbLeft(mbg.makeMultiBody(0, true, X_0_l));
  MultiBodyConfig mbcLeft(mbLeft);
  MultiBody mbRight(mbg.makeMultiBody(0, true, X_0_r));
  MultiBodyConfig mbcRight(mbRight);

  spring_estimator::SpringEstimator est;
  est.initArms({mbLeft, mbRight});

  PTransformd X_0_acc(Vector3d(0., 0., 0.5));
  PTransformd X_l_acc = X_0_acc*X_0_l.inv();
  PTransformd X_r_acc = X_0_acc*X_0_r.inv();

  mbLeft.transform(3, X_l_acc);
  mbRight.transform(3, X_r_acc);
  est.updateArms({mbLeft, mbRight});

  // check getter and setter
  est.leastSquareMinTol(1e-6);
  BOOST_CHECK_EQUAL(est.leastSquareMinTol(), 1e-6);
  est.leastSquareRelTol(2e-6);
  BOOST_CHECK_EQUAL(est.leastSquareRelTol(), 2e-6);
  est.projectorMinTol(3e-6);
  BOOST_CHECK_EQUAL(est.projectorMinTol(), 3e-6);
  est.projectorRelTol(4e-6);
  BOOST_CHECK_EQUAL(est.projectorRelTol(), 4e-6);

  // set the threshold
  // we just change the projectorRelTol threshold
  // to a lower value. This must avoid to block
  // lower priority task.
  est.leastSquareMinTol(1e-8);
  est.leastSquareRelTol(1e-8);
  est.projectorMinTol(1e-8);
  est.projectorRelTol(1e-2);

  // the estimator is stuck during 2000 iteration if the threshold is to high
  // est.projectorRelTol(1e-8);

  VectorXd q(4);
  q << 1., 1., 1., 1.;

  est.q(q);
  est.target(X_0_acc.rotation());
  est.target(sva::RotX(0.1));

  const std::size_t nrIter = 6000;
  std::vector<std::array<double, 4>> taskErrorLog(nrIter);

  internal::set_is_malloc_allowed(false);
  for(std::size_t i = 0; i < nrIter; ++i)
  {
    est.update(0.005, 1);
    taskErrorLog[i] = std::array<double, 4>{
      {est.taskError(0).norm(), est.taskError(1).norm(),
       est.taskError(2).norm(), est.taskError(3).norm()}};
  }
  internal::set_is_malloc_allowed(true);

  toPython<4>(taskErrorLog, "2arms.py");

  rbd::vectorToParam(est.q().segment(0, 2), mbcLeft.q);
  rbd::vectorToParam(est.q().segment(2, 2), mbcRight.q);

  rbd::forwardKinematics(mbLeft, mbcLeft);
  rbd::forwardKinematics(mbRight, mbcRight);

  testEndEffectors(mbcLeft, mbcRight, 3, est.target());

  // test the estimator with one arm.
  est.initArms({mbLeft});
  q.resize(2);
  q << 1., 1.;
  est.q(q);
  est.target(sva::RotX(0.1));

  internal::set_is_malloc_allowed(false);
  est.update(0.1, 200);
  internal::set_is_malloc_allowed(true);

  rbd::vectorToParam(est.q().segment(0, 2), mbcLeft.q);
  rbd::forwardKinematics(mbLeft, mbcLeft);
  BOOST_CHECK_SMALL(sva::rotationError(mbcLeft.bodyPosW[3].rotation(),
                    est.target(), 1e-7).norm(), 1e-6);
}


BOOST_AUTO_TEST_CASE(SpringEstimatorTestPrism)
{
  using namespace Eigen;
  using namespace sva;
  using namespace rbd;

  rbd::MultiBody mb;
  rbd::MultiBodyConfig mbc;
  rbd::MultiBodyGraph mbg;

  std::tie(mb, mbc, mbg) = makeThreeDofArm();

  PTransformd X_0_l(Vector3d(-0.1, 0., 0.));
  PTransformd X_0_r(Vector3d(0.1, 0., 0.));

  MultiBody mbLeft(mbg.makeMultiBody(0, true, X_0_l));
  MultiBodyConfig mbcLeft(mbLeft);
  MultiBody mbRight(mbg.makeMultiBody(0, true, X_0_r));
  MultiBodyConfig mbcRight(mbRight);

  spring_estimator::SpringEstimator est;
  est.initArms({mbLeft, mbRight});

  PTransformd X_0_acc(Vector3d(0., 0., 0.5));
  PTransformd X_l_acc = X_0_acc*X_0_l.inv();
  PTransformd X_r_acc = X_0_acc*X_0_r.inv();

  mbLeft.transform(4, X_l_acc);
  mbRight.transform(4, X_r_acc);
  est.updateArms({mbLeft, mbRight});

  // set the threshold
  // see explanation in prévious test
  est.leastSquareMinTol(1e-8);
  est.leastSquareRelTol(1e-8);
  est.projectorMinTol(1e-8);
  est.projectorRelTol(1e-2);

  VectorXd q(6);
  q << 1., 1., 0., 1., 1., 0.;


  // first we test X rotation target without joint target objective
  est.q(q);
  est.target(sva::RotX(0.1));

  const std::size_t nrIter = 6000;
  std::vector<std::array<double, 4>> taskErrorLog(nrIter);
  std::vector<std::array<double, 5>> taskErrorPrismLog(nrIter);
  std::vector<std::array<double, 6>> qLog(nrIter);

  internal::set_is_malloc_allowed(false);
  for(std::size_t i = 0; i < nrIter; ++i)
  {
    est.update(0.05, 1);
    taskErrorLog[i] = std::array<double, 4>{
      {est.taskError(0).norm(), est.taskError(1).norm(),
       est.taskError(2).norm(), est.taskError(3).norm()}};
    Eigen::Map<Eigen::Matrix<double, 6, 1>>(qLog[i].data()) = est.q();
  }
  internal::set_is_malloc_allowed(true);

  toPython<4>(taskErrorLog, "2arms3dof_x_err.py");
  toPython<6>(qLog, "2arms3dof_x_q.py");

  rbd::vectorToParam(est.q().segment(0, 3), mbcLeft.q);
  rbd::vectorToParam(est.q().segment(3, 3), mbcRight.q);
  rbd::forwardKinematics(mbLeft, mbcLeft);
  rbd::forwardKinematics(mbRight, mbcRight);

  testEndEffectors(mbcLeft, mbcRight, 4, est.target());
  // Z prisms are not at 0
  BOOST_CHECK_GE(std::abs(est.q()(2)), 0.001);
  BOOST_CHECK_GE(std::abs(est.q()(5)), 0.001);


  // same test but with rotation around Y axis
  est.initArms({mbLeft, mbRight});
  est.q(q);
  est.target(sva::RotY(0.1));

  internal::set_is_malloc_allowed(false);
  for(std::size_t i = 0; i < nrIter; ++i)
  {
    est.update(0.05, 1);
    taskErrorLog[i] = std::array<double, 4>{
      {est.taskError(0).norm(), est.taskError(1).norm(),
       est.taskError(2).norm(), est.taskError(3).norm()}};
    Eigen::Map<Eigen::Matrix<double, 6, 1>>(qLog[i].data()) = est.q();
  }
  internal::set_is_malloc_allowed(true);

  toPython<4>(taskErrorLog, "2arms3dof_y_err.py");
  toPython<6>(qLog, "2arms3dof_y_q.py");

  rbd::vectorToParam(est.q().segment(0, 3), mbcLeft.q);
  rbd::vectorToParam(est.q().segment(3, 3), mbcRight.q);
  rbd::forwardKinematics(mbLeft, mbcLeft);
  rbd::forwardKinematics(mbRight, mbcRight);

  testEndEffectors(mbcLeft, mbcRight, 4, est.target(), 1e-3);
  // Z prisms are not equals
  BOOST_CHECK_GE(std::abs(std::abs(est.q()(2)) - std::abs(est.q()(5))), 0.001);


  // Test X rot with prisms joint targeted to 0
  est.initArms({mbLeft, mbRight}, {{0,2,0.},{1,2,0.}});
  est.q(q);
  est.target(X_0_acc.rotation());
  est.target(sva::RotX(0.1));

  internal::set_is_malloc_allowed(false);
  for(std::size_t i = 0; i < nrIter; ++i)
  {
    est.update(0.005, 1);
    taskErrorPrismLog[i] = std::array<double, 5>{
      {est.taskError(0).norm(), est.taskError(1).norm(),
       est.taskError(2).norm(), est.taskError(3).norm(),
       est.taskError(4).norm()}};
    Eigen::Map<Eigen::Matrix<double, 6, 1>>(qLog[i].data()) = est.q();
  }
  internal::set_is_malloc_allowed(true);

  toPython<5>(taskErrorPrismLog, "2arms3dof_x_fixed_err.py");
  toPython<6>(qLog, "2arms3dof_x_fixed_q.py");

  rbd::vectorToParam(est.q().segment(0, 3), mbcLeft.q);
  rbd::vectorToParam(est.q().segment(3, 3), mbcRight.q);
  rbd::forwardKinematics(mbLeft, mbcLeft);
  rbd::forwardKinematics(mbRight, mbcRight);

  testEndEffectors(mbcLeft, mbcRight, 4, est.target());
  // prisms Z must be near zero
  BOOST_CHECK_LE(std::abs(est.q()(2)), 1e-4);
  BOOST_CHECK_LE(std::abs(est.q()(5)), 1e-4);


  // Test Y rot with prisms joint targeted to 0
  est.initArms({mbLeft, mbRight}, {{0,2,0.},{1,2,0.}});
  est.q(q);
  est.target(X_0_acc.rotation());
  est.target(sva::RotY(0.1));

  internal::set_is_malloc_allowed(false);
  for(std::size_t i = 0; i < nrIter; ++i)
  {
    est.update(0.005, 1);
    taskErrorPrismLog[i] = std::array<double, 5>{
      {est.taskError(0).norm(), est.taskError(1).norm(),
       est.taskError(2).norm(), est.taskError(3).norm(),
       est.taskError(4).norm()}};
    Eigen::Map<Eigen::Matrix<double, 6, 1>>(qLog[i].data()) = est.q();
  }
  internal::set_is_malloc_allowed(true);

  toPython<5>(taskErrorPrismLog, "2arms3dof_y_fixed_err.py");
  toPython<6>(qLog, "2arms3dof_y_fixed_q.py");

  rbd::vectorToParam(est.q().segment(0, 3), mbcLeft.q);
  rbd::vectorToParam(est.q().segment(3, 3), mbcRight.q);
  rbd::forwardKinematics(mbLeft, mbcLeft);
  rbd::forwardKinematics(mbRight, mbcRight);

  testEndEffectors(mbcLeft, mbcRight, 4, est.target(), 1e-3);
  // Z prisms are equals
  BOOST_CHECK_LE(std::abs(est.q()(2)) - std::abs(est.q()(5)), 1e-4);


  // test the estimator with one arm.
  est.initArms({mbLeft}, {{0,2,0.}});
  q.resize(3);
  q << 1., 1., 0.;
  est.q(q);
  est.target(sva::RotX(0.1));

  internal::set_is_malloc_allowed(false);
  est.update(0.1, 200);
  internal::set_is_malloc_allowed(true);

  rbd::vectorToParam(est.q().segment(0, 3), mbcLeft.q);
  rbd::forwardKinematics(mbLeft, mbcLeft);
  BOOST_CHECK_SMALL(sva::rotationError(mbcLeft.bodyPosW[3].rotation(),
                    est.target(), 1e-7).norm(), 1e-6);
}
