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

#pragma once

// includes
// std
#include <tuple>

// RBDyn
#include <RBDyn/Body.h>
#include <RBDyn/Joint.h>
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>
#include <RBDyn/MultiBodyGraph.h>


/// @return An simple rotoide XY  and prismatic Z arm with Z as up axis.
std::tuple<rbd::MultiBody, rbd::MultiBodyConfig, rbd::MultiBodyGraph>
makeThreeDofArm()
{
  using namespace Eigen;
  using namespace sva;
  using namespace rbd;

  MultiBodyGraph mbg;

  double mass = 1.;
  Matrix3d inertia = Matrix3d::Identity();
  Vector3d h = Vector3d::Zero();

  RBInertiad rbi(mass, h, inertia);

  Body b0(rbi, 0, "b0");
  Body b1(rbi, 1, "b1");
  Body b2(rbi, 2, "b2");
  Body b3(rbi, 3, "b3");
  Body b4(rbi, 4, "b4");

  mbg.addBody(b0);
  mbg.addBody(b1);
  mbg.addBody(b2);
  mbg.addBody(b3);
  mbg.addBody(b4);

  Joint j0(Joint::RevX, true, 0, "j0");
  Joint j1(Joint::RevY, true, 1, "j1");
  Joint j2(Joint::PrismZ, true, 2, "j2");
  Joint j3(Joint::Fixed, true, 3, "j3");

  mbg.addJoint(j0);
  mbg.addJoint(j1);
  mbg.addJoint(j2);
  mbg.addJoint(j3);


  PTransformd I(sva::PTransformd::Identity());
  PTransformd to(Vector3d(0., 0., 0.5));


  mbg.linkBodies(0, I, 1, I, 0);
  mbg.linkBodies(1, I, 2, I, 1);
  mbg.linkBodies(2, I, 3, I, 2);
  mbg.linkBodies(3, to, 4, I, 3);

  MultiBody mb = mbg.makeMultiBody(0, true);

  MultiBodyConfig mbc(mb);
  mbc.zero(mb);

  return std::make_tuple(mb, mbc, mbg);
}
