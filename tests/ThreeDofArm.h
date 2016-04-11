// Copyright 2015-2016 CNRS-UM LIRMM, CNRS-AIST JRL
//
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

  Body b0(rbi, "b0");
  Body b1(rbi, "b1");
  Body b2(rbi, "b2");
  Body b3(rbi, "b3");
  Body b4(rbi, "b4");

  mbg.addBody(b0);
  mbg.addBody(b1);
  mbg.addBody(b2);
  mbg.addBody(b3);
  mbg.addBody(b4);

  Joint j0(Joint::RevX, true, "j0");
  Joint j1(Joint::RevY, true, "j1");
  Joint j2(Joint::PrismZ, true, "j2");
  Joint j3(Joint::Fixed, true, "j3");

  mbg.addJoint(j0);
  mbg.addJoint(j1);
  mbg.addJoint(j2);
  mbg.addJoint(j3);


  PTransformd I(sva::PTransformd::Identity());
  PTransformd to(Vector3d(0., 0., 0.5));


  mbg.linkBodies("b0", I, "b1", I, "j0");
  mbg.linkBodies("b1", I, "b2", I, "j1");
  mbg.linkBodies("b2", I, "b3", I, "j2");
  mbg.linkBodies("b3", to, "b4", I, "j3");

  MultiBody mb = mbg.makeMultiBody("b0", true);

  MultiBodyConfig mbc(mb);
  mbc.zero(mb);

  return std::make_tuple(mb, mbc, mbg);
}
