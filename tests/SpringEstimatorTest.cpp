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


void toPython(const std::vector<std::array<double, 4>>& data,
              std::size_t dim, const std::string& name, std::ofstream& out)
{
  out << name << "= [";
  for(std::size_t i = 0; i < data.size(); ++i)
  {
    out << data[i][dim] << ",";
  }
  out << "]" << std::endl;
}


void toPython(const std::vector<std::array<double, 4>>& data,
              const std::string& filename)
{
  std::ofstream out(filename);
  toPython(data, 0, "t1Err", out);
  toPython(data, 1, "t2Err", out);
  toPython(data, 2, "t3Err", out);
  toPython(data, 3, "t4Err", out);
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
  internal::set_is_malloc_allowed(false);

  const std::size_t nrIter = 6000;
  std::vector<std::array<double, 4>> taskErrorLog(nrIter);
  for(std::size_t i = 0; i < nrIter; ++i)
  {
    est.update(0.005, 1);
    taskErrorLog[i] = std::array<double, 4>{
      {est.taskError(0).norm(), est.taskError(1).norm(),
       est.taskError(2).norm(), est.taskError(3).norm()}};
  }
  toPython(taskErrorLog, "2arms.py");

  internal::set_is_malloc_allowed(true);

  rbd::vectorToParam(est.q().segment(0, 2), mbcLeft.q);
  rbd::vectorToParam(est.q().segment(2, 2), mbcRight.q);

  rbd::forwardKinematics(mbLeft, mbcLeft);
  rbd::forwardKinematics(mbRight, mbcRight);

  BOOST_CHECK_SMALL((mbcLeft.bodyPosW[3].translation() -
                     mbcRight.bodyPosW[3].translation()).norm(), 1e-6);
  BOOST_CHECK_SMALL(sva::rotationError(mbcLeft.bodyPosW[3].rotation(),
                    mbcRight.bodyPosW[3].rotation(), 1e-7).norm(), 1e-6);
  BOOST_CHECK_SMALL(sva::rotationError(mbcLeft.bodyPosW[3].rotation(),
                    est.target(), 1e-7).norm(), 1e-6);
  BOOST_CHECK_SMALL(sva::rotationError(mbcRight.bodyPosW[3].rotation(),
                    est.target(), 1e-7).norm(), 1e-6);

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

  VectorXd q(6);
  q << 1., 1., 0., 1., 1., 0.;

  est.q(q);
  est.target(X_0_acc.rotation());
  est.target(sva::RotX(0.1));
  internal::set_is_malloc_allowed(false);

  const std::size_t nrIter = 6000;
  std::vector<std::array<double, 4>> taskErrorLog(nrIter);
  for(std::size_t i = 0; i < nrIter; ++i)
  {
    est.update(0.005, 1);
    taskErrorLog[i] = std::array<double, 4>{
      {est.taskError(0).norm(), est.taskError(1).norm(),
       est.taskError(2).norm(), est.taskError(3).norm()}};
  }
  toPython(taskErrorLog, "2arms3dof.py");

  internal::set_is_malloc_allowed(true);

  rbd::vectorToParam(est.q().segment(0, 3), mbcLeft.q);
  rbd::vectorToParam(est.q().segment(3, 3), mbcRight.q);

  rbd::forwardKinematics(mbLeft, mbcLeft);
  rbd::forwardKinematics(mbRight, mbcRight);

  BOOST_CHECK_SMALL((mbcLeft.bodyPosW[4].translation() -
                     mbcRight.bodyPosW[4].translation()).norm(), 1e-6);
  BOOST_CHECK_SMALL(sva::rotationError(mbcLeft.bodyPosW[4].rotation(),
                    mbcRight.bodyPosW[4].rotation(), 1e-7).norm(), 1e-6);
  BOOST_CHECK_SMALL(sva::rotationError(mbcLeft.bodyPosW[4].rotation(),
                    est.target(), 1e-7).norm(), 1e-6);
  BOOST_CHECK_SMALL(sva::rotationError(mbcRight.bodyPosW[4].rotation(),
                    est.target(), 1e-7).norm(), 1e-6);
}
