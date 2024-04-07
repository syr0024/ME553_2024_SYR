//
// Created by Jemin Hwangbo on 2022/03/17.
//

#ifndef ME553_2024_SOLUTIONS_EXERCISE1_20243406_HPP_
#define ME553_2024_SOLUTIONS_EXERCISE1_20243406_HPP_

#include <Eigen/Core>
#include <iostream>

class ArticulatedSystem {
public:
  ArticulatedSystem() {
    joint_pos_w_.setZero();
    joint_ori_w_.setIdentity();
  };
  
  Eigen::Matrix3d rpy2Rot(const Eigen::Vector3d& rpy) {
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
  Eigen::Matrix3d quat2Rot(const Eigen::Vector4d& q) {
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
  Eigen::Vector3d quat2Rpy(const Eigen::Vector4d& q) {
    double roll, pitch, yaw;
    double sinr_cosp = 2 * (q[0] * q[1] + q[2] * q[3]);
    double cosr_cosp = 1 - 2 * (q[1] * q[1] + q[2] * q[2]);
    roll = std::atan2(sinr_cosp, cosr_cosp);
    
    double sinp = 2 * (q[0] * q[2] - q[3] * q[1]);
    if (std::abs(sinp) >= 1)
      pitch = std::copysign(M_PI / 2, sinp); // use 90 degrees if out of range
    else
      pitch = std::asin(sinp);
    
    double siny_cosp = 2 * (q[0] * q[3] + q[1] * q[2]);
    double cosy_cosp = 1 - 2 * (q[2] * q[2] + q[3] * q[3]);
    yaw = std::atan2(siny_cosp, cosy_cosp);
    
    return Eigen::Vector3d(roll, pitch, yaw);
  }
  
  Eigen::Vector3d forwardKinematics(const Eigen::Vector3d &parent_pos, const Eigen::Vector3d& parent_ori) {
    
    joint_pos_w_ = parent_pos + rpy2Rot(parent_ori)*joint_pos_w_;
    joint_ori_w_ = rpy2Rot(parent_ori) * joint_ori_w_;
    
    return joint_pos_w_;
  }
  
  Eigen::VectorXd gc_;
  Eigen::Vector3d joint_pos_w_;
  Eigen::Matrix3d joint_ori_w_;
  
};

/// do not change the name of the method
inline Eigen::Vector3d getEndEffectorPosition (const Eigen::VectorXd& gc) {
  
  ArticulatedSystem articulatedSystem_;
  
  // initialized result
  Eigen::Vector3d posWe; posWe.setZero();
  
  // number of joint
  int joint_num = 10;
  Eigen::VectorXd origin_ori; origin_ori.setZero(joint_num*3 +3);
  Eigen::VectorXd origin_pos; origin_pos.setZero(joint_num*3 +3);
  
  // number of dynamic joint
  int dynamic_joint_num = 3;
  Eigen::VectorXd theta; theta.setZero(dynamic_joint_num);
  
  // add base pose and orientation
  origin_pos.head(3) = gc.head(3);
  origin_ori.head(3) << articulatedSystem_.quat2Rpy(gc.block<4,1>(3, 0));
  
  // add revolute joint difference
  for(int i = 0; i < dynamic_joint_num; i++){
    theta[i] = gc(13+i);
  }
  
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
  
  origin_pos.tail(joint_num*3) << -0.2999, 0.104, 0.0,
                                    0,0,0,
                                    0,0,0,
                                    -0.0599, 0.08381, 0.0,
                                    0,0,0,
                                    0,0,0,
                                    -0.0, 0.1003, -0.285,
                                    0,0,0,
                                    0,0,0,
                                    -0.08795, 0.01305, -0.33797;
  
  // calculate end effector position
  for (int i = joint_num*3 +3; i > 0; i -= 3){
    articulatedSystem_.forwardKinematics(origin_pos.block<3,1>(i-3, 0), origin_ori.block<3,1>(i-3, 0));
  }
  
  return articulatedSystem_.joint_pos_w_; /// replace this
}

#endif // ME553_2022_SOLUTIONS_EXERCISE1_STUDENTID_HPP_