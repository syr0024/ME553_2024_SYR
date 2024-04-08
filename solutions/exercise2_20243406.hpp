#pragma once

#include <Eigen/Core>
#include <iostream>

class ArticulatedSystem {
public:
  ArticulatedSystem() {
    joint_pos_w_.setZero();
    joint_ori_w_.setIdentity();
  };
  
  void reset () {
    joint_pos_w_.setZero();
    joint_ori_w_.setIdentity();
  }
  
  Eigen::Matrix3d rpy2Rot (const Eigen::Vector3d& rpy) {
    Eigen::Matrix3d Rx, Ry, Rz;
    
    Rx << 1, 0, 0,
      0, cos(rpy[0]), -sin(rpy[0]),
      0, sin(rpy[0]), cos(rpy[0]);
    Ry << cos(rpy[1]), 0, sin(rpy[1]),
      0, 1, 0,
      -sin(rpy[1]), 0, cos(rpy[1]);
    Rz << cos(rpy[2]), -sin(rpy[2]), 0,
      sin(rpy[2]), cos(rpy[2]), 0,
      0, 0, 1;
    
    return Rx*Ry*Rz;
  }
  Eigen::Matrix3d quat2Rot (const Eigen::Vector4d& q) {
    Eigen::Matrix3d rot;
    rot << 2*(q(0)*q(0) + q(1)*q(1)) - 1,
      2*(q(1)*q(2) - q(0)*q(3)),
      2*(q(1)*q(3) + q(0)*q(2)),
      2*(q(1)*q(2) + q(0)*q(3)),
      2*(q(0)*q(0) + q(2)*q(2)) - 1,
      2*(q(2)*q(3) - q(0)*q(1)),
      2*(q(1)*q(3) - q(0)*q(2)),
      2*(q(2)*q(3) + q(0)*q(1)),
      2*(q(0)*q(0) + q(3)*q(3)) - 1;
    return rot;
  }
  
  Eigen::MatrixXd skewSymmetric (const Eigen::Vector3d& a) {
    Eigen::Matrix3d ax;
    ax << 0, -a[2], a[1],
          a[2], 0, -a[0],
          -a[1], a[0], 0;
    return ax;
  }
  
  Eigen::Vector3d forwardKinematics (const Eigen::Vector3d &parent_pos, const Eigen::Vector3d& parent_ori) {
    
    joint_pos_w_ = parent_pos + rpy2Rot(parent_ori)*joint_pos_w_;
    joint_ori_w_ = rpy2Rot(parent_ori) * joint_ori_w_;
    
    return joint_pos_w_;
  }
  
  Eigen::Vector3d joint_pos_w_;
  Eigen::Matrix3d joint_ori_w_;
  
};

class RobotParameters {
public:
  
  RobotParameters() {};
  
  void reset () {
    origin_ori.setZero(joint_num*3 +3);
    origin_pos.setZero(joint_num*3 +3);
    r_W_LH_HAA.setZero();
    p_W_LH_HAA.setZero();
    r_W_LH_HFE.setZero();
    p_W_LH_HFE.setZero();
    r_W_LH_KFE.setZero();
    p_W_LH_HFE.setZero();
  }
  
  void step (const Eigen::VectorXd& gc) {
    
    // number of dynamic joint
    Eigen::VectorXd theta; theta.setZero(dynamic_joint_num);
    
    // add base pose and orientation
    origin_pos.head(3) = gc.head(3);
    
    // add revolute joint difference
    for(int i = 0; i < dynamic_joint_num; i++){
      theta[i] = gc(13+i);
    }
    
    origin_pos.tail(joint_num*3) << -0.2999, 0.104, 0.0, // base_LH_HAA
      0,0,0, // LH_HAA (revolute)
      0,0,0, // LH_HIP_LH_hip_fixed
      -0.0599, 0.08381, 0.0, // LH_hip_fixed_LH_HFE
      0,0,0, // LH_HFE (revolute)
      0,0,0, // LH_THIGH_LH_thigh_fixed
      -0.0, 0.1003, -0.285, // LH_thigh_fixed_LH_KFE
      0,0,0, // LH_KFE (revolute)
      0,0,0, // LH_shank_LH_shank_fixed
      -0.08795, 0.01305, -0.33797; // LH_shank_fixed_LH_FOOT
    
    // Write down from the top joint origin
    origin_ori.tail(joint_num*3) << -2.61799387799, 0, -3.14159265359,
      -theta[0], 0, 0,
      -2.61799387799,0,-3.14159265359,
      0, 0, 1.57079632679,
      theta[1],0,0,
      0, 0, -1.57079632679,
      0, 0, 1.57079632679,
      theta[2], 0, 0,
      0, 0, -1.57079632679,
      0, 0, 0;
    
    Eigen::Matrix3d ROB; ROB = articulatedSystem_.quat2Rot(gc.block<4,1>(3, 0));
    
    articulatedSystem_.reset();
    for (int i = 3*3; i > 3; i -= 3){
      articulatedSystem_.forwardKinematics(origin_pos.block<3,1>(i-3, 0), origin_ori.block<3,1>(i-3, 0));
    }
    articulatedSystem_.joint_pos_w_ = gc.head(3) + ROB*articulatedSystem_.joint_pos_w_;
    r_W_LH_HAA = articulatedSystem_.joint_pos_w_;
    p_W_LH_HAA = - (ROB * articulatedSystem_.joint_ori_w_).block<3,1>(0,0);
    
    articulatedSystem_.reset();
    for (int i = 6*3; i > 3; i -= 3){
      articulatedSystem_.forwardKinematics(origin_pos.block<3,1>(i-3, 0), origin_ori.block<3,1>(i-3, 0));
    }
    articulatedSystem_.joint_pos_w_ = gc.head(3) + ROB*articulatedSystem_.joint_pos_w_;
    r_W_LH_HFE = articulatedSystem_.joint_pos_w_;
    p_W_LH_HFE = (ROB * articulatedSystem_.joint_ori_w_).block<3,1>(0,0);
    
    articulatedSystem_.reset();
    for (int i = 9*3; i > 3; i -= 3){
      articulatedSystem_.forwardKinematics(origin_pos.block<3,1>(i-3, 0), origin_ori.block<3,1>(i-3, 0));
    }
    articulatedSystem_.joint_pos_w_ = gc.head(3) + ROB*articulatedSystem_.joint_pos_w_;
    r_W_LH_KFE = articulatedSystem_.joint_pos_w_;
    p_W_LH_KFE = (ROB * articulatedSystem_.joint_ori_w_).block<3,1>(0,0);
    
    articulatedSystem_.reset();
    for (int i = joint_num*3 +3; i > 3; i -= 3){
      articulatedSystem_.forwardKinematics(origin_pos.block<3,1>(i-3, 0), origin_ori.block<3,1>(i-3, 0));
    }
    articulatedSystem_.joint_pos_w_ = gc.head(3) + ROB*articulatedSystem_.joint_pos_w_;
    r_W_e = articulatedSystem_.joint_pos_w_;
  }
  
  ArticulatedSystem articulatedSystem_;
  /// variables
  Eigen::VectorXd origin_ori;
  Eigen::VectorXd origin_pos;
  Eigen::Vector3d r_W_LH_HAA;
  Eigen::Vector3d p_W_LH_HAA;
  Eigen::Vector3d r_W_LH_HFE;
  Eigen::Vector3d p_W_LH_HFE;
  Eigen::Vector3d r_W_LH_KFE;
  Eigen::Vector3d p_W_LH_KFE;
  Eigen::Vector3d r_W_e;
  int joint_num = 10;
  int dynamic_joint_num = 3;
  
};

/// do not change the name of the method
inline Eigen::Vector3d getFootLinearVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
  
  ArticulatedSystem articulatedSystem_;
  RobotParameters robotParameters_;
  robotParameters_.reset();
  robotParameters_.step(gc);
  
  Eigen::MatrixXd positional_jacobain; positional_jacobain.setZero(3, 6+robotParameters_.dynamic_joint_num);
  positional_jacobain.block<3,3>(0,0) = Eigen::Matrix3d::Identity();
  positional_jacobain.block<3,3>(0,3) = - articulatedSystem_.skewSymmetric((robotParameters_.r_W_e - gc.head(3)));
  positional_jacobain.block<3,1>(0,6) = robotParameters_.p_W_LH_HAA.cross(robotParameters_.r_W_e - robotParameters_.r_W_LH_HAA);
  positional_jacobain.block<3,1>(0,7) = robotParameters_. p_W_LH_HFE.cross(robotParameters_.r_W_e - robotParameters_.r_W_LH_HFE);
  positional_jacobain.block<3,1>(0,8) = robotParameters_.p_W_LH_KFE.cross(robotParameters_.r_W_e - robotParameters_.r_W_LH_KFE);
  
  Eigen::VectorXd joint_gv; joint_gv.setZero(9);
  joint_gv << gv.head(6), gv.block<3,1>(12,0);
  
  return positional_jacobain * joint_gv;
}

/// do not change the name of the method
inline Eigen::Vector3d getFootAngularVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
  
  ArticulatedSystem articulatedSystem_;
  RobotParameters robotParameters_;
  robotParameters_.reset();
  robotParameters_.step(gc);
  
  Eigen::MatrixXd velocity_jacobian; velocity_jacobian.setZero(3, 6+robotParameters_.dynamic_joint_num);
  velocity_jacobian.block<3,3>(0,0) = Eigen::Matrix3d::Zero();
  velocity_jacobian.block<3,3>(0,3) = Eigen::Matrix3d::Identity();
  velocity_jacobian.block<3,1>(0,6) = robotParameters_.p_W_LH_HAA;
  velocity_jacobian.block<3,1>(0,7) = robotParameters_.p_W_LH_HFE;
  velocity_jacobian.block<3,1>(0,8) = robotParameters_.p_W_LH_KFE;
  
  Eigen::VectorXd joint_gv; joint_gv.setZero(9);
  joint_gv << gv.head(6), gv.block<3,1>(12,0);
  
  return velocity_jacobian*joint_gv;
}

inline Eigen::Vector3d getEndEffectorPosition (const Eigen::VectorXd& gc) {
  ArticulatedSystem articulatedSystem_;
  RobotParameters robotParameters_;
  robotParameters_.reset();
  robotParameters_.step(gc);
  
  return robotParameters_.r_W_e;
}