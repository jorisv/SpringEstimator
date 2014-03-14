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
#include <vector>

// RBDyn
#include <RBDyn/Jacobian.h>
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>

namespace spring_estimator
{

class SpringEstimator
{
public:
  void initArms(const std::vector<rbd::MultiBody>& arms);
  void updateArms(const std::vector<rbd::MultiBody>& arms);

  void target(const Eigen::Matrix3d& t);
  const Eigen::Matrix3d& target() const;

  void q(const Eigen::VectorXd& q);
  const Eigen::VectorXd& q() const;

  const Eigen::VectorXd& qd() const;

  double update(double timeStep, int nrIter);

private:
  void updateArmsData();
  double update1Arm(double timeStep, int nrIter);
  double updateNArm(double timeStep, int nrIter);

private:
  struct ArmData
  {
    rbd::MultiBody mb;
    rbd::MultiBodyConfig mbc;
    rbd::Jacobian jac;
    Eigen::MatrixXd jacMat;
    int endEffectorIndex;
  };

  struct Stage1Data
  {
    Eigen::VectorXd err;
    Eigen::MatrixXd jac;
    Eigen::JacobiSVD<Eigen::MatrixXd> jacSvd;
    Eigen::VectorXd svdSingular;
    Eigen::MatrixXd preResult;
    Eigen::MatrixXd jacPseudoInv;
    Eigen::VectorXd qd;
  };

  struct Stage2Data
  {
    Eigen::Vector3d err;
    Eigen::MatrixXd jac;
    Eigen::MatrixXd projector;
    Eigen::MatrixXd projectorJac;
    Eigen::JacobiSVD<Eigen::MatrixXd> projectorJacSvd;
    Eigen::VectorXd svdSingular;
    Eigen::MatrixXd preResult;
    Eigen::MatrixXd projectorJacPseudoInv;
  };

private:
  std::vector<ArmData> arms_;
  Eigen::Matrix3d target_;

  Eigen::VectorXd q_;
  Eigen::VectorXd qd_;

  Stage1Data s1data;
  Stage2Data s2data;
};


} // spring estimator
