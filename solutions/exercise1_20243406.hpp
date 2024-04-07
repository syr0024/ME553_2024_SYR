//
// Created by Jemin Hwangbo on 2022/03/17.
//

#ifndef ME553_2024_SOLUTIONS_EXERCISE1_20243406_HPP_
#define ME553_2024_SOLUTIONS_EXERCISE1_20243406_HPP_

#include <Eigen/Core>
#include <iostream>

Eigen::Matrix3d rpy_to_rot(double r, double p, double y) {
  Eigen::Matrix3d Rx, Ry, Rz;
  
  Rx << 1, 0, 0,
        0, cos(r), -sin(r),
        0, sin(r), cos(r);
  Ry << cos(p), 0, sin(p),
        0, 1, 0,
        -sin(p), 0, cos(p);
  Rz << cos(y), -sin(y), 0,
        sin(y), cos(y), 0,
        0, 0, 1;
  
  return Rx*Ry*Rz;
}

Eigen::Matrix3d quat_to_rot(Eigen::Vector4d q){
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

Eigen::Vector3d quat_to_rpy(const Eigen::Vector4d& q) {
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

/// do not change the name of the method
inline Eigen::Vector3d getEndEffectorPosition (const Eigen::VectorXd& gc) {
  
  // initialized result
  Eigen::Vector3d posWe; posWe.setZero();
  
  // number of joint
  int joint_num = 10;
  Eigen::VectorXd origin_rot; origin_rot.setZero((joint_num+1)*3);
  Eigen::VectorXd origin_pos; origin_pos.setZero((joint_num+1)*3);
  
  // number of dynamic joint
  int revolute_joint_num = 3;
  Eigen::VectorXd theta; theta.setZero(revolute_joint_num);
  
  // add base pose and orientation
  origin_pos.head(3) = gc.head(3);
  origin_rot.head(3) << quat_to_rpy(gc.block<4,1>(3, 0));
  
  // add revolute joint difference
  for(int i = 0; i < revolute_joint_num; i++){
    theta[i] = gc(13+i);
  }
  
  // Write down from the top joint origin
  origin_rot.tail(joint_num*3) << -2.61799387799, 0, -3.14159265359,
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
  for (int i = (joint_num+1)*3; i > 0; i -= 3){
    posWe = origin_pos.block<3,1>(i-3, 0) + rpy_to_rot(origin_rot[i-3], origin_rot[i-2], origin_rot[i-1])*posWe;
  }
  
  return posWe; /// replace this
}

#endif // ME553_2022_SOLUTIONS_EXERCISE1_STUDENTID_HPP_