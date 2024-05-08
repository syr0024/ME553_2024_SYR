#pragma once

Eigen::Matrix3d hat (Eigen::Vector3d v) {
  Eigen::Matrix3d S;
  S <<  0, -v.z(),  v.y(),
        v.z(),  0, -v.x(),
        -v.y(),  v.x(), 0;
  return S;
}

void getComposite (const Eigen::Vector3d& r1_com, const Eigen::Vector3d& r2_com,
                    const double& m1, const double& m2,
                    const Eigen::Matrix3d& I1, const Eigen::Matrix3d& I2,
                    Eigen::Vector3d& r_com, double& mass, Eigen::Matrix3d& I_new) {
  r_com = (m1*r1_com + m2*r2_com)/(m1+m2);
  Eigen::Vector3d r1 = r1_com - r_com;
  Eigen::Vector3d r2 = r2_com - r_com;

  mass = m1 + m2;
  
  I_new = I1 + I2 - m1*hat(r1)*hat(r1) - m2*hat(r2)*hat(r2);
}

Eigen::Matrix3d rpy2rot (const double& r, const double& p, const double& y) {
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

Eigen::Matrix3d rpy2rot (const Eigen::Vector3d& rpy) {
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
  
  return Rz*Ry*Rx;
}

Eigen::Matrix3d quat2rot (const Eigen::Vector4d& q) {
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

Eigen::Vector3d getLinkCom (const Eigen::Vector3d &parent_pos, const Eigen::Vector3d& parent_ori,
                                 const Eigen::Vector3d &joint_pos, const Eigen::Vector3d& joint_ori,
                                 const Eigen::Vector3d &com_pos) {

  return parent_pos + rpy2rot(parent_ori)*(joint_pos + rpy2rot(joint_ori)*(com_pos));
}

void getBaseLinkCom (const Eigen::Vector3d &parent_pos, const Eigen::Vector4d& parent_ori,
                            const Eigen::Vector3d &joint_pos, const Eigen::Vector3d& joint_ori,
                            const Eigen::Vector3d &com_pos, const Eigen::Vector3d& com_ori,
                            Eigen::Vector3d& pos_w, Eigen::Matrix3d& ori_w) {
  
  pos_w = parent_pos + quat2rot(parent_ori)*(joint_pos + rpy2rot(joint_ori)*(com_pos));
  ori_w = quat2rot(parent_ori) * rpy2rot(joint_ori) * rpy2rot(com_ori);
}



//class Joint{
//public:
//  Joint();
//
//  enum class Type{
//    FIXED,
//    FLOATING,
//    PRISMATIC,
//    REVOLUTE,
//  }type=Type::REVOLUTE;
//
//  enum class JointName {
//    BASE,
//    FR_HIP,
//    FR_THIGH,
//    FR_CALF,
//    FL_HIP,
//    FL_THIGH,
//    FL_CALF,
//    RR_HIP,
//    RR_THIGH,
//    RR_CALF,
//    RL_HIP,
//    RL_THIGH,
//    RL_CALF,
//    ELSE,
//  }jointName=JointName::ELSE;
//
//  Eigen::Vector3d get_world_pos();
//  Eigen::Vector3d get_world_rot();
//
//};

class Body{
public:
  Body();
  Eigen::Matrix3d getBaseCompositeInertia(const Eigen::VectorXd& gc,
                                          const Eigen::VectorXd& joint_ori,
                                          const Eigen::VectorXd& joint_pos,
                                          const Eigen::VectorXd& com_ori,
                                          const Eigen::VectorXd& com_pos,
                                          const Eigen::VectorXd& mass,
                                          const Eigen::VectorXd& inertia) ;
  Eigen::Matrix3d getInertia (const Eigen::VectorXd& inertia_vec);
  
public:
  Eigen::Vector3d base_com;
  double base_mass = 0;
  Eigen::Matrix3d base_inertia;
};

Body::Body() {};

Eigen::Matrix3d Body::getBaseCompositeInertia  (const Eigen::VectorXd& gc,
                                                const Eigen::VectorXd& joint_ori,
                                                const Eigen::VectorXd& joint_pos,
                                                const Eigen::VectorXd& com_ori,
                                                const Eigen::VectorXd& com_pos,
                                                const Eigen::VectorXd& mass,
                                                const Eigen::VectorXd& inertia) {
  Eigen::Vector3d gc_pos = gc.head(3);
  Eigen::Matrix3d gc_rot = quat2rot(gc.segment(3,4));
  base_com.setZero(); base_inertia.setZero();
  
  for(int i=0; i<mass.size(); i++){
    Eigen::Vector3d com_pos_w; Eigen::Matrix3d com_rot_w; Eigen::Matrix3d inertia_w; Eigen::Vector3d base_com_new;
    com_pos_w = gc_pos + gc_rot*(joint_pos.segment(i*3, 3) + rpy2rot(joint_ori.segment(i*3, 3))*com_pos.segment(i*3, 3));
    com_rot_w = gc_rot * rpy2rot(joint_ori.segment(i*3, 3)) *  rpy2rot(com_ori.segment(i*3, 3));
    inertia_w = com_rot_w * getInertia(inertia.segment(i*6, 6)) * com_rot_w.transpose();
    
    base_com_new = (base_mass*base_com + mass(i)*com_pos_w) / (base_mass+mass(i));
    Eigen::Vector3d r1 = base_com - base_com_new;
    Eigen::Vector3d r2 = com_pos_w - base_com_new;
    
    base_inertia = base_inertia + inertia_w - base_mass*hat(r1)*hat(r1) - mass(i)*hat(r2)*hat(r2);

    base_mass = base_mass + mass(i);
    
    base_com = base_com_new;
  }
};

Eigen::Matrix3d Body::getInertia (const Eigen::VectorXd& inertia_vec) {
  Eigen::Matrix3d inertia;
  inertia << inertia_vec[0], inertia_vec[1], inertia_vec[2],
              inertia_vec[1], inertia_vec[3], inertia_vec[4],
              inertia_vec[2], inertia_vec[4], inertia_vec[5];
  return inertia;
};

/// do not change the name of the method
inline Eigen::MatrixXd getMassMatrix (const Eigen::VectorXd& gc) {

  // Step 1. get center of mass of each body
  // << Base >>
  Body body;
  
  Eigen::VectorXd base_joint_ori, base_joint_pos, base_com_ori, base_com_pos, base_mass, base_inertia;
  int base_child_num = 11;
  base_joint_ori.setZero(base_child_num*3); base_joint_pos.setZero(base_child_num*3);
  base_com_ori.setZero(base_child_num*3); base_com_pos.setZero(base_child_num*3);
  base_mass.setZero(base_child_num); base_inertia.setZero(base_child_num*6);
  base_joint_ori << 0,0,0, 0,0,0, 0,0,3.14159265359, 0,0,0, 0,0,0, 0.0,0.0,-1.57079632679, 0,0,0,
                    2.61799387799,0,0,  -2.61799387799,0,0,  -2.61799387799,0,-3.14159265359, 2.61799387799,0,-3.14159265359;
  
  base_joint_pos << 0,0,0, 0.4145,0,0, -0.4145,0,0, 0,0,0, 0.343,0.0,-0.07, -0.364,0.0,0.0735+0.0687, 0.0,0.0,0.0,
                    0.2999,0.104,0.0,  0.2999,-0.104,0.0,  -0.2999,0.104,0.0, -0.2999,-0.104,0.0;
  
  base_com_ori << 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
                  0,0,0,  0,0,0,  0,0,0, 0,0,0;
  
  base_com_pos << -0.018,-0.002,0.024,  0.042,-0.001,0.004, 0.042,-0.001,0.004, -0.00067,-0.00023,-0.03362, -0.003,0.0,0.005, -0.012,0.001,-0.008, 0.116,0.0,0.0758,
                  -0.063,7e-05,0.00046,  -0.063,7e-05,0.00046,  -0.063,7e-05,0.00046, -0.063,7e-05,0.00046;
                  
  base_mass << 6.222, 0.73, 0.73, 5.53425, 0.065, 0.695, 0.142, 2.04, 2.04, 2.04, 2.04;
  
  base_inertia << 0.017938806, 0.00387963, 0.001500772, 0.370887745, 6.8963e-05, 0.372497653,
    0.005238611, 1.7609e-05, 7.2167e-05, 0.002643098, 1.9548e-05, 0.004325938,
    0.005238611, 1.7609e-05, 7.2167e-05, 0.002643098, 1.9548e-05, 0.004325938,
    0.00749474794, 0.00016686282, 7.82763e-05, 0.0722338913, 1.42902e-06, 0.07482717535,
    0.00063283, 0.0, 3.45e-07, 0.00110971, 0.0, 0.00171883,
    0.000846765, 6.9565e-05, 0.00027111, 0.001367583, 5.8984e-05, 0.001363673,
    0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
    0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827,
    0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827,
    0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827,
    0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827;
  
  body.getBaseCompositeInertia(gc, base_joint_ori, base_joint_pos, base_com_ori, base_com_pos, base_mass, base_inertia);
  
 
  std::cout << "\n base_com ... \n" << body.base_com << std::endl;
  std::cout << "\n base_mass ... \n" << body.base_mass << std::endl;
  std::cout << "\n base_inertia ... \n" << body.base_inertia << std::endl;
  
  // LF_HAA
  
  
  return Eigen::MatrixXd::Ones(18,18);
}