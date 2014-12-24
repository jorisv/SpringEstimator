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

  struct TaskData
  {
    TaskData(int dim, int dof);

    // input
    Eigen::VectorXd err;
    Eigen::MatrixXd jac;

    // result
    Eigen::MatrixXd pseudoInv;
    Eigen::VectorXd qd;

    // buffer
    Eigen::VectorXd errMinusPrev;
    Eigen::MatrixXd projectorJac;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd;
    Eigen::VectorXd svdSingular;
    Eigen::MatrixXd prePseudoInv;
  };

  struct ProjectorData
  {
    ProjectorData(int dim, int dof);

    // input
    Eigen::MatrixXd jacA;

    // output
    Eigen::MatrixXd projector;

    // buffer
    Eigen::JacobiSVD<Eigen::MatrixXd> svd;
    Eigen::VectorXd svdSingular;
    Eigen::MatrixXd preProjector;
  };

private:
  static void solveT1(TaskData& task1);
  static void solveTN(const TaskData& taskPrev, const ProjectorData& projPrev,
                      TaskData& taskN);
  static void projector(const TaskData& taskPrev,
                        const ProjectorData& projPrev,
                        ProjectorData& proj);

private:
  std::vector<ArmData> arms_;
  Eigen::Matrix3d target_;

  Eigen::VectorXd q_;
  Eigen::VectorXd qd_;

  std::vector<TaskData> tasks_;
  std::vector<ProjectorData> projs_;
};


} // spring estimator
