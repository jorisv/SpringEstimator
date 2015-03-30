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

// associated header
#include "SpringEstimator.h"

// includes
// RBDyn
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>

namespace spring_estimator
{

SpringEstimator::TaskData::TaskData(int dim, int dof):
  err(Eigen::VectorXd::Zero(dim)),
  jac(Eigen::MatrixXd::Zero(dim, dof)),
  pseudoInv(Eigen::MatrixXd::Zero(dof, dim)),
  qd(Eigen::VectorXd::Zero(dof)),
  errMinusPrev(Eigen::VectorXd::Zero(dim)),
  projectorJac(Eigen::MatrixXd::Zero(dim, dof)),
  svd(dim, dof, Eigen::ComputeThinU | Eigen::ComputeThinV),
  svdSingular(Eigen::VectorXd::Zero(std::min(dim, dof))),
  prePseudoInv(Eigen::MatrixXd::Zero(dof, svdSingular.rows()))
{}


SpringEstimator::ProjectorData::ProjectorData(int dim, int dof):
  jacA(Eigen::MatrixXd::Zero(dim, dof)),
  projector(Eigen::MatrixXd::Zero(dof, dof)),
  svd(dim, dof, Eigen::ComputeFullU | Eigen::ComputeFullV),
  svdSingular(Eigen::VectorXd::Zero(dof)),
  preProjector(Eigen::MatrixXd::Zero(dof, dof))
{}


void SpringEstimator::initArms(const std::vector<rbd::MultiBody>& arms,
  const std::vector<JointTarget>& jointTarget)
{
  arms_.clear();
  jTargetData_.clear();
  tasks_.clear();
  projs_.clear();

  int dof = 0;
  std::vector<int> armQDBegin;
  for(const rbd::MultiBody& mb: arms)
  {
    armQDBegin.push_back(dof);
    dof += mb.nrDof();
    rbd::MultiBodyConfig mbc(mb);
    mbc.zero(mb);
    rbd::forwardKinematics(mb, mbc);
    rbd::forwardVelocity(mb, mbc);

    int endEffectorId = mb.body(mb.nrBodies() - 1).id();
    rbd::Jacobian jac(mb, endEffectorId);
    arms_.push_back({mb, mbc, jac, Eigen::MatrixXd(6, jac.dof()),
                     mb.bodyIndexById(endEffectorId)});
  }

  q_.setZero(dof, 1);
  qd_.setZero(dof, 1);

  Eigen::MatrixXd jointTargetJac(
    Eigen::MatrixXd::Zero(jointTarget.size(), dof));

  // file joint targets data and joint target jacobian
  for(std::size_t i = 0; i < jointTarget.size(); ++i)
  {
    const JointTarget& jt = jointTarget[i];
    ArmData& ad = arms_[jt.armIndex];
    /// TODO manage arm with dof != params
    /// TODO manage joint.dof > 1
    int jointIndex = ad.mb.jointIndexById(jt.jointId);
    int jointIndexInQD = armQDBegin[jt.armIndex] +
      ad.mb.jointPosInDof(jointIndex);

    jointTargetJac(i, jointIndexInQD) = 1.;

    jTargetData_.push_back({jointIndexInQD, jt.target});
  }

  // create the first dummy projector
  projs_.emplace_back(0, dof);

  int cumDim = 0;

  // only add the end effector null translation
  // and rotation error task if the is more than 1 arm
  if(arms.size() > 1)
  {
    // end effector translation task
    tasks_.emplace_back((arms.size() - 1)*3, dof);
    cumDim += int(tasks_.back().jac.rows());
    projs_.emplace_back(cumDim, dof);

    // end effector rotation task
    tasks_.emplace_back((arms.size() - 1)*3, dof);
    cumDim += int(tasks_.back().jac.rows());
    projs_.emplace_back(cumDim, dof);
  }

  // minimize rotation error to target
  tasks_.emplace_back(3, dof);
  cumDim += int(tasks_.back().jac.rows());
  projs_.emplace_back(cumDim, dof);

  // joints target (jacobian is constant but error is computed each iter)
  if(jTargetData_.size() > 0)
  {
    tasks_.emplace_back(jTargetData_.size(), dof);
    tasks_.back().jac = jointTargetJac;
    cumDim += int(tasks_.back().jac.rows());
    projs_.emplace_back(cumDim, dof);
  }

  // minimize joints velocity
  // the jacobian and error are constants
  tasks_.emplace_back(dof, dof);
  tasks_.back().err.setZero();
  tasks_.back().jac.setIdentity();
}


void SpringEstimator::updateArms(const std::vector<rbd::MultiBody>& arms)
{
  for(std::size_t i = 0; i < arms.size(); ++i)
  {
    arms_[i].mb = arms[i];
  }
}


void SpringEstimator::target(const Eigen::Matrix3d& t)
{
  target_ = t;
}


const Eigen::Matrix3d& SpringEstimator::target() const
{
  return target_;
}


void SpringEstimator::q(const Eigen::VectorXd& q)
{
  q_ = q;
}


const Eigen::VectorXd& SpringEstimator::q() const
{
  return q_;
}


const Eigen::VectorXd& SpringEstimator::qd() const
{
  return qd_;
}


const Eigen::VectorXd& SpringEstimator::taskError(std::size_t taskIndex) const
{
  return tasks_[taskIndex].err;
}


void SpringEstimator::leastSquareMinTol(double lstsqMinTol)
{
  lstsqMinTol_ = lstsqMinTol;
}


double SpringEstimator::leastSquareMinTol() const
{
  return lstsqMinTol_;
}


void SpringEstimator::leastSquareRelTol(double lstsqRelTol)
{
  lstsqRelTol_ = lstsqRelTol;
}


double SpringEstimator::leastSquareRelTol() const
{
  return lstsqRelTol_;
}


void SpringEstimator::projectorMinTol(double projMinTol)
{
  projMinTol_ = projMinTol;
}


double SpringEstimator::projectorMinTol() const
{
  return projMinTol_;
}


void SpringEstimator::projectorRelTol(double projRelTol)
{
  projRelTol_ = projRelTol;
}


double SpringEstimator::projectorRelTol() const
{
  return projRelTol_;
}


void pseudoInverse(const Eigen::MatrixXd& jac,
                   Eigen::JacobiSVD<Eigen::MatrixXd>& svd,
                   Eigen::VectorXd& svdSingular,
                   Eigen::MatrixXd& prePseudoInv,
                   Eigen::MatrixXd& result,
                   double epsilon=std::numeric_limits<double>::epsilon(),
                   double minTol=1e-8)
{
  svd.compute(jac, Eigen::ComputeThinU | Eigen::ComputeThinV);

  // singular values are sorted in decreasing order
  // so the first one is the max one
  double tolerance =
      epsilon*double(std::max(jac.cols(), jac.rows()))*
      std::abs(svd.singularValues()[0]);
  tolerance = std::max(tolerance, minTol);

  svdSingular = ((svd.singularValues().array().abs() > tolerance).
      select(svd.singularValues().array().inverse(), 0.));

  prePseudoInv.noalias() = svd.matrixV()*svdSingular.asDiagonal();
  result.noalias() = prePseudoInv*svd.matrixU().adjoint();
}


void projectorFromSvd(const Eigen::MatrixXd& jac,
                      Eigen::JacobiSVD<Eigen::MatrixXd>& svd,
                      Eigen::VectorXd& svdSingular,
                      Eigen::MatrixXd& preResult,
                      Eigen::MatrixXd& result,
                      double epsilon=std::numeric_limits<double>::epsilon(),
                      double minTol=1e-8)
{
  // we are force to compute the Full matrix because of
  // the nullspace matrix computation
  svd.compute(jac, Eigen::ComputeFullU | Eigen::ComputeFullV);

  double tolerance =
      epsilon*double(std::max(jac.cols(), jac.rows()))*
      std::abs(svd.singularValues()[0]);
  tolerance = std::max(tolerance, minTol);

  svdSingular.setOnes();
  for(int i = 0; i < svd.singularValues().rows(); ++i)
  {
    svdSingular[i] = svd.singularValues()[i] > tolerance ? 0. : 1.;
  }

  preResult.noalias() = svd.matrixV()*svdSingular.asDiagonal();
  result.noalias() = preResult*svd.matrixV().adjoint();
}


void projectorDummy(const Eigen::MatrixXd& pseudoInv,
                    const Eigen::MatrixXd& jac,
                    Eigen::MatrixXd& result)
{
  result.setIdentity();
  result.noalias() -= pseudoInv*jac;
}


double SpringEstimator::update(double timeStep, int nrIter)
{
  if(arms_.size() > 1)
  {
    return updateNArm(timeStep, nrIter);
  }
  else
  {
    return update1Arm(timeStep, nrIter);
  }
}


void SpringEstimator::updateArmsData()
{
  int dof = 0;
  for(std::size_t i = 0; i < arms_.size(); ++i)
  {
    rbd::vectorToParam(q_.segment(dof, arms_[i].jac.dof()), arms_[i].mbc.q);
    rbd::forwardKinematics(arms_[i].mb, arms_[i].mbc);
    arms_[i].jacMat = arms_[i].jac.jacobian(arms_[i].mb, arms_[i].mbc);
    dof += arms_[i].jac.dof();
  }
}


double SpringEstimator::update1Arm(double timeStep, int nrIter)
{
  for(int iter = 0; iter < nrIter; ++iter)
  {
    updateArmsData();

    // compute task 3 error and jacobian
    // minimize distance between the arm 0 end effector and the target
    const sva::PTransformd& arm0 =
        arms_[0].mbc.bodyPosW[arms_[0].endEffectorIndex];
    tasks_[0].err.noalias() = sva::rotationError(target_, arm0.rotation(), 1e-7);
    tasks_[0].jac.block(0, 0, 3, arms_[0].jac.dof()).noalias() =
        arms_[0].jacMat.block(0, 0, 3, arms_[0].jac.dof());

    for(std::size_t i = 0; i < jTargetData_.size(); ++i)
    {
      const JointTargetData& jtd = jTargetData_[i];
      tasks_[1].err(i) = q_(jtd.jointIndexInQ) - jtd.target;
    }

    // solve all the tasks
    // minimize the distance and minimum joint velocity
    solveT1(tasks_[0], lstsqMinTol_, lstsqRelTol_);
    for(std::size_t i = 1; i < tasks_.size(); ++i)
    {
      projector(tasks_[i - 1], projs_[i - 1], projs_[i], projMinTol_,
                projRelTol_);
      solveTN(tasks_[i - 1], projs_[i], tasks_[i], lstsqMinTol_, lstsqRelTol_);
    }

    qd_.noalias() = tasks_.back().qd;
    q_.noalias() -= qd_*timeStep;
  }
  return 0.;
}


double SpringEstimator::updateNArm(double timeStep, int nrIter)
{
  for(int iter = 0; iter < nrIter; ++iter)
  {
    updateArmsData();

    // compute task 1 and 2 error and jacobian
    // (end effector position/orientation error)
    int dof = 0;
    for(std::size_t i = 1; i < arms_.size(); ++i)
    {
      dof += arms_[i-1].jac.dof();
      const sva::PTransformd& arm0 =
          arms_[0].mbc.bodyPosW[arms_[0].endEffectorIndex];
      const sva::PTransformd& armi =
          arms_[i].mbc.bodyPosW[arms_[i].endEffectorIndex];

      tasks_[0].err.segment((i-1)*3, 3).noalias() =
          arm0.translation() - armi.translation();
      tasks_[1].err.segment((i-1)*3, 3).noalias() =
          sva::rotationError(armi.rotation(), arm0.rotation(), 1e-7);

      tasks_[0].jac.block((i-1)*3, 0, 3, arms_[0].jac.dof()).noalias() =
          arms_[0].jacMat.block(3, 0, 3, arms_[0].jac.dof());
      tasks_[1].jac.block((i-1)*3, 0, 3, arms_[0].jac.dof()).noalias() =
          arms_[0].jacMat.block(0, 0, 3, arms_[0].jac.dof());

      tasks_[0].jac.block((i-1)*3, dof, 3, arms_[i].jac.dof()).noalias() =
          -arms_[i].jacMat.block(3, 0, 3, arms_[i].jac.dof());
      tasks_[1].jac.block((i-1)*3, dof, 3, arms_[i].jac.dof()).noalias() =
          -arms_[i].jacMat.block(0, 0, 3, arms_[i].jac.dof());
    }

    // compute task 3 error and jacobian
    // minimize distance between the arm 0 end effector and the target
    const sva::PTransformd& arm0 =
        arms_[0].mbc.bodyPosW[arms_[0].endEffectorIndex];
    tasks_[2].err.noalias() = sva::rotationError(target_, arm0.rotation(), 1e-7);
    tasks_[2].jac.block(0, 0, 3, arms_[0].jac.dof()).noalias() =
        arms_[0].jacMat.block(0, 0, 3, arms_[0].jac.dof());

    for(std::size_t i = 0; i < jTargetData_.size(); ++i)
    {
      const JointTargetData& jtd = jTargetData_[i];
      tasks_[3].err(i) = q_(jtd.jointIndexInQ) - jtd.target;
    }

    // solve all the tasks
    solveT1(tasks_[0], lstsqMinTol_, lstsqRelTol_);
    for(std::size_t i = 1; i < tasks_.size(); ++i)
    {
      projector(tasks_[i - 1], projs_[i - 1], projs_[i], projMinTol_,
                projRelTol_);
      solveTN(tasks_[i - 1], projs_[i], tasks_[i], lstsqMinTol_, lstsqRelTol_);
    }

    qd_.noalias() = tasks_.back().qd;
    q_.noalias() -= qd_*timeStep;
  }

  return 0.;
}


void SpringEstimator::solveT1(TaskData& task1,
                              double lstsqMinTol, double lstsqRelTol)
{
  // compute the least square solution of task1
  pseudoInverse(task1.jac, task1.svd,
                task1.svdSingular, task1.prePseudoInv,
                task1.pseudoInv, lstsqRelTol, lstsqMinTol);
  task1.qd.noalias() = task1.pseudoInv*task1.err;
}


void SpringEstimator::solveTN(const TaskData& taskPrev,
                              const ProjectorData& projPrev, TaskData& taskN,
                              double lstsqMinTol, double lstsqRelTol)
{
  // compute the least square solution of taskN with the
  // taskN jacobian project into taskPrev nullspace
  taskN.projectorJac.noalias() = taskN.jac*projPrev.projector;

  taskN.errMinusPrev.noalias() = taskN.err;
  taskN.errMinusPrev.noalias() -= taskN.jac*taskPrev.qd;

  pseudoInverse(taskN.projectorJac, taskN.svd,
                taskN.svdSingular, taskN.prePseudoInv,
                taskN.pseudoInv, lstsqRelTol, lstsqMinTol);
  taskN.qd.noalias() = taskPrev.qd;
  taskN.qd.noalias() += taskN.pseudoInv*taskN.errMinusPrev;
}


void SpringEstimator::projector(const TaskData& taskPrev,
                                const ProjectorData& projPrev,
                                ProjectorData& proj,
                                double projMinTol, double projRelTol)
{
  // fill the jacobian with the last projector and the last task
  const auto dof = proj.jacA.cols();
  const auto projPrevR = projPrev.jacA.rows();
  proj.jacA.block(0, 0, projPrevR, dof) = projPrev.jacA;
  proj.jacA.block(projPrevR, 0, taskPrev.jac.rows(), dof) = taskPrev.jac;

  // compute the projector from ProjectorData
  projectorFromSvd(proj.jacA, proj.svd, proj.svdSingular, proj.preProjector,
                   proj.projector, projRelTol, projMinTol);
}

} // spring estimator
