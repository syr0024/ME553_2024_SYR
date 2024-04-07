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
  return Rz*Ry*Rx;
}

/// do not change the name of the method
inline Eigen::Vector3d getEndEffectorPosition (const Eigen::VectorXd& gc) {
  
  Eigen::Vector3d base_pos = gc.head(3);
  Eigen::Vector4d base_quat = gc.block<4,1>(3, 0);
  Eigen::Matrix3d base_rot;
  base_rot << 2*(base_quat(0)*base_quat(0) + base_quat(1)*base_quat(1)) - 1, 2*(base_quat(1)*base_quat(2) - base_quat(0)*base_quat(3)), 2*(base_quat(1)*base_quat(3) + base_quat(0)*base_quat(2)),
    2*(base_quat(1)*base_quat(2) + base_quat(0)*base_quat(3)), 2*(base_quat(0)*base_quat(0) + base_quat(2)*base_quat(2)) - 1, 2*(base_quat(2)*base_quat(3) - base_quat(0)*base_quat(1)),
    2*(base_quat(1)*base_quat(3) - base_quat(0)*base_quat(2)), 2*(base_quat(2)*base_quat(3) + base_quat(0)*base_quat(1)), 2*(base_quat(0)*base_quat(0) + base_quat(3)*base_quat(3)) - 1;
  
  double theta1 = gc(13); double theta2 = gc(14); double theta3 = gc(15);
  
  Eigen::Matrix3d rot_base2HAA, rot_HAA2HIP, rot_HIP2hip_fixed, rot_hip_fixed2HFE, rot_HFE2THIGH, rot_THIGH2thigh_fixed;
  Eigen::Matrix3d rot_thigh_fixed2KFE, rot_KFE2SHANK, rot_SHANK2shank_fixed, rot_shank_fixed2FOOT;
  Eigen::Vector3d pos_base2HAA, pos_HAA2HIP, pos_HIP2hip_fixed, pos_hip_fixed2HFE, pos_HFE2THIGH, pos_THIGH2thigh_fixed;
  Eigen::Vector3d pos_thigh_fixed2KFE, pos_KFE2SHANK, pos_SHANK2shank_fixed, pos_shank_fixed2FOOT;

  rot_base2HAA = rpy_to_rot(-2.61799387799,0,-3.14159265359);
  pos_base2HAA << -0.2999, 0.104, 0.0;
  
  rot_HAA2HIP = rpy_to_rot(-theta1,0,0);
  pos_HAA2HIP << 0,0,0;
  
  rot_HIP2hip_fixed = rpy_to_rot(-2.61799387799,0,-3.14159265359);
  pos_HIP2hip_fixed << 0,0,0;
  
  rot_hip_fixed2HFE = rpy_to_rot(0, 0, 1.57079632679);
  pos_hip_fixed2HFE << -0.0599, 0.08381, 0.0;
  
  rot_HFE2THIGH = rpy_to_rot(theta2,0,0);
  pos_HFE2THIGH << 0,0,0;
  
  rot_THIGH2thigh_fixed = rpy_to_rot(0, 0, -1.57079632679);
  pos_THIGH2thigh_fixed << 0,0,0;
  
  rot_thigh_fixed2KFE = rpy_to_rot(0, 0, 1.57079632679);
  pos_thigh_fixed2KFE << -0.0, 0.1003, -0.285;
  
  rot_KFE2SHANK = rpy_to_rot(theta3, 0, 0);
  pos_KFE2SHANK << 0,0,0;
  
  rot_SHANK2shank_fixed = rpy_to_rot(0, 0, -1.57079632679);
  pos_SHANK2shank_fixed << 0,0,0;
  
  rot_shank_fixed2FOOT = rpy_to_rot(0,0,0);
  pos_shank_fixed2FOOT << -0.08795, 0.01305, -0.33797;
  
  Eigen::Vector3d posWe = base_pos + base_rot * (pos_base2HAA + rot_base2HAA *
    (pos_HAA2HIP + rot_HAA2HIP * (pos_HIP2hip_fixed + rot_HIP2hip_fixed *
    (pos_hip_fixed2HFE + rot_hip_fixed2HFE * (pos_HFE2THIGH + rot_HFE2THIGH *
    (pos_THIGH2thigh_fixed + rot_THIGH2thigh_fixed * (pos_thigh_fixed2KFE + rot_thigh_fixed2KFE *
    (pos_KFE2SHANK + rot_KFE2SHANK * (pos_SHANK2shank_fixed + rot_SHANK2shank_fixed *
    (pos_shank_fixed2FOOT))))))))));
  
  return posWe; /// replace this
}

#endif // ME553_2022_SOLUTIONS_EXERCISE1_STUDENTID_HPP_