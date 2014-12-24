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
  VectorXd q(4);
  q << 1., 1., 1., 1.;

  est.q(q);
  est.target(X_0_acc.rotation());
  est.target(sva::RotX(0.1));
  internal::set_is_malloc_allowed(false);
  est.update(0.005, 6000);
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
