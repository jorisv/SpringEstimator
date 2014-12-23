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


void SpringEstimator::initArms(const std::vector<rbd::MultiBody>& arms)
{
  arms_.clear();

  int dof = 0;
  for(const rbd::MultiBody& mb: arms)
  {
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

  s1data.err.setZero((arms.size() - 1)*6, 1);
  s1data.jac.setZero((arms.size() - 1)*6, dof);
  s1data.jacSvd =
      Eigen::JacobiSVD<Eigen::MatrixXd>(s1data.jac.rows(), s1data.jac.cols(),
                                        Eigen::ComputeThinU | Eigen::ComputeThinV);
  s1data.svdSingular.setZero(std::min(s1data.jac.rows(), s1data.jac.cols()), 1);
  s1data.preResult.setZero(s1data.jac.cols(), s1data.svdSingular.rows());
  s1data.jacPseudoInv.setZero(s1data.jac.cols(), s1data.jac.rows());
  s1data.qd.setZero(dof, 1);

  s2data.err.setZero();
  s2data.jac.setZero(3, dof);
  s2data.projector.setZero(dof, dof);
  s2data.projectorJac.setZero(3, dof);
  s2data.projectorJacSvd =
      Eigen::JacobiSVD<Eigen::MatrixXd>(s2data.projectorJac.rows(),
                                        s2data.projectorJac.cols(),
                                        Eigen::ComputeThinU | Eigen::ComputeThinV);
  s2data.svdSingular.setZero(std::min(s2data.projectorJac.rows(),
                                      s2data.projectorJac.cols()), 1);
  s2data.preResult.setZero(s2data.projectorJac.cols(), s2data.svdSingular.rows());
  s2data.projectorJacPseudoInv.setZero(s2data.projectorJac.cols(),
                                       s2data.projectorJac.rows());

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



double pseudoInverse(const Eigen::MatrixXd& jac,
                   Eigen::JacobiSVD<Eigen::MatrixXd>& svd,
                   Eigen::VectorXd& svdSingular,
                   Eigen::MatrixXd& preResult,
                   Eigen::MatrixXd& result,
                   double epsilon=std::numeric_limits<double>::epsilon())
{
  svd.compute(jac, Eigen::ComputeThinU | Eigen::ComputeThinV);

  double tolerance =
      epsilon*double(std::max(jac.cols(), jac.rows()))*
      svd.singularValues().array().abs().maxCoeff();

  svdSingular = ((svd.singularValues().array().abs() > tolerance).
      select(svd.singularValues().array().inverse(), 0.));

  preResult.noalias() = svd.matrixV()*svdSingular.asDiagonal();
  result.noalias() = preResult*svd.matrixU().adjoint();

  return tolerance;
}


void projectorFromSvd(Eigen::JacobiSVD<Eigen::MatrixXd>& svd,
                 Eigen::VectorXd& svdSingular,
                 Eigen::MatrixXd& preResult,
                 Eigen::MatrixXd& result,
                 double tolerance)
{
  for(int i = 0; i < svdSingular.rows(); ++i)
  {
    svdSingular[i] = svdSingular[i] > tolerance ? 0. : 1.;
  }

  preResult.noalias() = svd.matrixV()*svdSingular.asDiagonal();
  result.noalias() = preResult*svd.matrixV().adjoint();
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

    const sva::PTransformd& arm0 = arms_[0].mbc.bodyPosW[arms_[0].endEffectorIndex];
    s2data.err.noalias() = sva::rotationError(target_, arm0.rotation(), 1e-7);
    s2data.jac.block(0, 0, 3, arms_[0].jac.dof()).noalias() =
        arms_[0].jacMat.block(0, 0, 3, arms_[0].jac.dof());

    pseudoInverse(s2data.jac, s2data.projectorJacSvd, s2data.svdSingular,
                  s2data.preResult, s2data.projectorJacPseudoInv, 1e-8);
    qd_.noalias() = s2data.projectorJacPseudoInv*s2data.err;

    q_.noalias() -= qd_*timeStep;
  }
  return 0.;
}


double SpringEstimator::updateNArm(double timeStep, int nrIter)
{
  for(int iter = 0; iter < nrIter; ++iter)
  {
    updateArmsData();

    // stage 1 end effector position/orientation error
    int dof = 0;
    for(std::size_t i = 1; i < arms_.size(); ++i)
    {
      dof += arms_[i-1].jac.dof();
      const sva::PTransformd& arm0 = arms_[0].mbc.bodyPosW[arms_[0].endEffectorIndex];
      const sva::PTransformd& armi = arms_[i].mbc.bodyPosW[arms_[i].endEffectorIndex];

      s1data.err.segment((i-1)*6 + 0, 3).noalias() =
          sva::rotationError(armi.rotation(), arm0.rotation(), 1e-7);
      s1data.err.segment((i-1)*6 + 3, 3).noalias() =
          arm0.translation() - armi.translation();

      s1data.jac.block((i-1)*6 + 0, 0, 3, arms_[0].jac.dof()).noalias() =
          arms_[0].jacMat.block(0, 0, 3, arms_[0].jac.dof());
      s1data.jac.block((i-1)*6 + 3, 0, 3, arms_[0].jac.dof()).noalias() =
          arms_[0].jacMat.block(3, 0, 3, arms_[0].jac.dof());

      s1data.jac.block((i-1)*6 + 0, dof, 3, arms_[i].jac.dof()).noalias() =
          -arms_[i].jacMat.block(0, 0, 3, arms_[i].jac.dof());
      s1data.jac.block((i-1)*6 + 3, dof, 3, arms_[i].jac.dof()).noalias() =
          -arms_[i].jacMat.block(3, 0, 3, arms_[i].jac.dof());
    }
    double tolerance = pseudoInverse(s1data.jac, s1data.jacSvd,
                                     s1data.svdSingular, s1data.preResult,
                                     s1data.jacPseudoInv, 1e-8);
    s1data.qd.noalias() = s1data.jacPseudoInv*s1data.err;

    // stage 2 arm0 orientation error with target
    const sva::PTransformd& arm0 = arms_[0].mbc.bodyPosW[arms_[0].endEffectorIndex];
    s2data.err.noalias() = sva::rotationError(target_, arm0.rotation(), 1e-7);
    s2data.jac.block(0, 0, 3, arms_[0].jac.dof()).noalias() =
        arms_[0].jacMat.block(0, 0, 3, arms_[0].jac.dof());

    // compute the projector into stage 1 jacobian nullspace
    projectorFromSvd(s1data.jacSvd, s1data.svdSingular, s1data.preResult,
                     s2data.projector, tolerance);
    s2data.projectorJac.noalias() = s2data.jac*s2data.projector;

    Eigen::Vector3d errS2MinusS1 = s2data.err;
    errS2MinusS1.noalias() -= s2data.jac*s1data.qd;

    pseudoInverse(s2data.projectorJac, s2data.projectorJacSvd,
                  s2data.svdSingular, s2data.preResult,
                  s2data.projectorJacPseudoInv, 1e-8);

    qd_.noalias() = s1data.qd;
    qd_.noalias() += s2data.projectorJacPseudoInv*errS2MinusS1;

    q_.noalias() -= qd_*timeStep;
  }

  return 0.;
}

} // spring estimator
