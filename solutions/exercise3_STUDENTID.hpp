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



class Joint{
public:
  Joint();
  void set_gc(const Eigen::VectorXd& gc);
  void jointKinematics();
  
  
  
  Eigen::VectorXd gc_;
  std::vector<Eigen::Vector3d> joint_pos_w_;
  std::vector<Eigen::Matrix3d> joint_rot_w_;
  
private:
  static void forwardKinematics(const Eigen::Vector3d& parent_pos_w, const Eigen::Vector3d& parent_ori_w,
                         const Eigen::Vector3d& joint_pos, const Eigen::Matrix3d& joint_ori,
                         Eigen::Vector3d& joint_pos_w, Eigen::Matrix3d& joint_ori_w);
  int revolute_joint_num_, joint_num_;
  std::vector<Eigen::VectorXd> joint_ori_, joint_pos_;
};

Joint::Joint() {
  revolute_joint_num_ = 12;
  joint_num_ = 3; // 전체 joint 개수 (fixed 포함)
  
  // lf [0]
  Eigen::VectorXd lf_joint_ori; lf_joint_ori.setZero((3*3+1)*3);
  Eigen::VectorXd lf_joint_pos; lf_joint_pos.setZero((3*3+1)*3);
  lf_joint_ori << 2.61799387799,0,0, 0,0,0, -2.61799387799,0,0,
                  0,0,1.57079632679, 0,0,0, 0,0,-1.57079632679,
                  0,0,1.57079632679, 0,0,0, 0,0,-1.57079632679,
                  0,0,0;
  lf_joint_pos << 0.2999,0.104,0.0, 0,0,0, 0,0,0,
                  0.0599,0.08381,0.0, 0,0,0, 0,0,0,
                  0.0,0.1003,-0.285, 0,0,0, 0,0,0,
                  0.08795,0.01305,-0.33797;
  joint_ori_.push_back(lf_joint_ori);
  joint_pos_.push_back(lf_joint_pos);
  
  // rf [1]
  Eigen::VectorXd rf_joint_ori; rf_joint_ori.setZero((3*3+1)*3);
  Eigen::VectorXd rf_joint_pos; rf_joint_pos.setZero((3*3+1)*3);
  rf_joint_ori << -2.61799387799,0,0, 0,0,0, 2.61799387799,0,0,
    0,0,-1.57079632679, 0,0,0, 0,0,1.57079632679,
    0,0,-1.57079632679, 0,0,0, 0,0,1.57079632679,
    0,0,0;
  rf_joint_pos << 0.2999,-0.104,0, 0,0,0, 0,0,0,
    0.0599,-0.08381,0, 0,0,0, 0,0,0,
    0,-0.1003,-0.285, 0,0,0, 0,0,0,
    0.08795,-0.01305,-0.33797;
  joint_ori_.push_back(rf_joint_ori);
  joint_pos_.push_back(rf_joint_pos);
  
  // lh [2]
  Eigen::VectorXd lh_joint_ori; lh_joint_ori.setZero((3*3+1)*3);
  Eigen::VectorXd lh_joint_pos; lh_joint_pos.setZero((3*3+1)*3);
  lh_joint_ori << -2.61799387799,0,-3.14159265359, 0,0,0, -2.61799387799,0,-3.14159265359,
    0,0,1.57079632679, 0,0,0, 0,0,-1.57079632679,
    0,0,1.57079632679, 0,0,0, 0,0,-1.57079632679,
    0,0,0;
  lh_joint_pos << -0.2999,0.104,0, 0,0,0, 0,0,0,
    -0.0599,0.08381,0, 0,0,0, 0,0,0,
    0,0.1003,-0.285, 0,0,0, 0,0,0,
    -0.08795,0.01305,-0.33797;
  joint_ori_.push_back(lh_joint_ori);
  joint_pos_.push_back(lh_joint_pos);
  
  // rh [3]
  Eigen::VectorXd rh_joint_ori; rh_joint_ori.setZero((3*3+1)*3);
  Eigen::VectorXd rh_joint_pos; rh_joint_pos.setZero((3*3+1)*3);
  rh_joint_ori << 2.61799387799,0,-3.14159265359, 0,0,0, 2.61799387799,0,-3.14159265359,
    0,0,-1.57079632679, 0,0,0, 0,0,1.57079632679,
    0,0,-1.57079632679, 0,0,0, 0,0,1.57079632679,
    0,0,0;
  rh_joint_pos << -0.2999,-0.104,0, 0,0,0, 0,0,0,
    -0.0599,-0.08381,0, 0,0,0, 0,0,0,
    0,-0.1003,-0.285, 0,0,0, 0,0,0,
    -0.08795,-0.01305,-0.33797;
  joint_ori_.push_back(rh_joint_ori);
  joint_pos_.push_back(rh_joint_pos);
  
};

void Joint::set_gc(const Eigen::VectorXd& gc) {
  gc_ = gc;
  joint_ori_[0].segment(3,1) = gc_.segment(7,1);
  joint_ori_[0].segment(12,1) = gc_.segment(8,1);
  joint_ori_[0].segment(21,1) = gc_.segment(9,1);
  joint_ori_[1].segment(3,1) = gc_.segment(10,1);
  joint_ori_[1].segment(12,1) = -gc_.segment(11,1);
  joint_ori_[1].segment(21,1) = -gc_.segment(12,1);
  joint_ori_[2].segment(3,1) = -gc_.segment(13,1);
  joint_ori_[2].segment(12,1) = gc_.segment(14,1);
  joint_ori_[2].segment(21,1) = gc_.segment(15,1);
  joint_ori_[3].segment(3,1) = -gc_.segment(16,1);
  joint_ori_[3].segment(12,1) = -gc_.segment(17,1);
  joint_ori_[3].segment(21,1) = -gc_.segment(18,1);
};

void Joint::jointKinematics() {
  for (int i=0; i < joint_pos_.size(); i++) {
    for (int j=0; j < joint_pos_[i].size(); j+=3) { // j: tree 전체의 joint 개수
      Eigen::Vector3d joint_pos_w; Eigen::Matrix3d joint_rot_w;
      joint_pos_w.setZero(); joint_rot_w.setIdentity();
      for (int k = j + 3; k>=3; k-=3) {
        forwardKinematics(joint_pos_[i].segment(k-3,3), joint_ori_[i].segment(k-3,3),
                          joint_pos_w, joint_rot_w,
                          joint_pos_w, joint_rot_w);
      }
      joint_pos_w = gc_.head(3) + quat2rot(gc_.segment(3,4))*joint_pos_w;
      joint_rot_w = quat2rot(gc_.segment(3,4)) * joint_rot_w;
      joint_pos_w_.push_back(joint_pos_w);
      joint_rot_w_.push_back(joint_rot_w);
    }
  }
}

void Joint::forwardKinematics(const Eigen::Vector3d& parent_pos_w, const Eigen::Vector3d& parent_ori_w,
                              const Eigen::Vector3d& joint_pos, const Eigen::Matrix3d& joint_ori,
                              Eigen::Vector3d& joint_pos_w, Eigen::Matrix3d& joint_ori_w) {
  joint_pos_w = parent_pos_w + rpy2rot(parent_ori_w)*joint_pos;
  joint_ori_w = rpy2rot(parent_ori_w) * joint_ori;
}




class Body{
public:
  Body();
  void getBaseCompositeInertia(const Eigen::VectorXd& gc,
                                const Eigen::VectorXd& joint_ori,
                                const Eigen::VectorXd& joint_pos,
                                const Eigen::VectorXd& com_ori,
                                const Eigen::VectorXd& com_pos,
                                const Eigen::VectorXd& mass,
                                const Eigen::VectorXd& inertia) ;
  void getLegCompositeInertia(const std::vector<Eigen::Matrix3d>& joint_ori_w,
                              const std::vector<Eigen::Vector3d>& joint_pos_w,
                              const Eigen::VectorXd& com_ori,
                              const Eigen::VectorXd& com_pos,
                              const Eigen::VectorXd& mass,
                              const Eigen::VectorXd& inertia);
  Eigen::Matrix3d getInertia (const Eigen::VectorXd& inertia_vec);
  
public:
  Eigen::Vector3d base_com;
  double base_mass = 0;
  Eigen::Matrix3d base_inertia;
  std::vector<Eigen::Vector3d> leg_com_w;
  std::vector<Eigen::Matrix3d> leg_inertia_w;
  std::vector<double> leg_mass_;
};

Body::Body() {};

void Body::getBaseCompositeInertia (const Eigen::VectorXd& gc,
                                    const Eigen::VectorXd& joint_ori,
                                    const Eigen::VectorXd& joint_pos,
                                    const Eigen::VectorXd& com_ori,
                                    const Eigen::VectorXd& com_pos,
                                    const Eigen::VectorXd& mass,
                                    const Eigen::VectorXd& inertia) {
  Eigen::Vector3d gc_pos = gc.head(3);
  Eigen::Matrix3d gc_rot = quat2rot(gc.segment(3, 4));
  base_com.setZero();
  base_inertia.setZero();
  
  for (int i = 0; i < mass.size(); i++) {
    Eigen::Vector3d com_pos_w;
    Eigen::Matrix3d com_rot_w;
    Eigen::Matrix3d inertia_w;
    Eigen::Vector3d base_com_new;
    com_pos_w = gc_pos + gc_rot * (joint_pos.segment(i * 3, 3) +
                                   rpy2rot(joint_ori.segment(i * 3, 3)) * com_pos.segment(i * 3, 3));
    com_rot_w = gc_rot * rpy2rot(joint_ori.segment(i * 3, 3)) * rpy2rot(com_ori.segment(i * 3, 3));
    inertia_w = com_rot_w * getInertia(inertia.segment(i * 6, 6)) * com_rot_w.transpose();
    
    base_com_new = (base_mass * base_com + mass(i) * com_pos_w) / (base_mass + mass(i));
    Eigen::Vector3d r1 = base_com - base_com_new;
    Eigen::Vector3d r2 = com_pos_w - base_com_new;
    
    base_inertia = base_inertia + inertia_w - base_mass * hat(r1) * hat(r1) - mass(i) * hat(r2) * hat(r2);
    
    base_mass = base_mass + mass(i);
    
    base_com = base_com_new;
  };
};
  
void Body::getLegCompositeInertia(const std::vector<Eigen::Matrix3d>& joint_ori_w,
                                  const std::vector<Eigen::Vector3d>& joint_pos_w,
                                  const Eigen::VectorXd& com_ori,
                                  const Eigen::VectorXd& com_pos,
                                  const Eigen::VectorXd& mass,
                                  const Eigen::VectorXd& inertia) {
  int num = static_cast<int>(mass.size());
  for(int i=num-1; i>=0; i--) {
    Eigen::Vector3d com_pos_w, leg_com_new;
    Eigen::Matrix3d com_rot_w, inertia_w, leg_inertia_new;
    com_pos_w = joint_pos_w[i] + joint_ori_w[i] * com_pos.segment(i*3, 3);
    com_rot_w = joint_ori_w[i] * rpy2rot(com_ori.segment(i*3, 3));
    inertia_w = com_rot_w * getInertia(inertia.segment(i*6, 6)) * com_rot_w.transpose();
    
    if(i == num-1) {
      leg_com_w.push_back(com_pos_w);
      leg_inertia_w.push_back(inertia_w);
      leg_mass_.push_back(mass(i));
    }
    else {
      leg_com_new = (leg_mass_.back() * leg_com_w.back() + mass(i) * com_pos_w) / (leg_mass_.back() + mass(i));
      Eigen::Vector3d r1 = leg_com_w.back() - leg_com_new;
      Eigen::Vector3d r2 = com_pos_w - leg_com_new;
      
      leg_inertia_new = leg_inertia_w.back() + inertia_w - leg_mass_.back() * hat(r1) * hat(r1) - mass(i) * hat(r2) * hat(r2);
      leg_inertia_w.push_back(leg_inertia_new);
      leg_com_w.push_back(leg_com_new);
      leg_mass_.push_back(leg_mass_.back() + mass(i));
    }
    std::cout << leg_com_w.back().transpose() << std::endl;
    std::cout << leg_mass_.back() << std::endl;
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
  
  Body body;
  Joint joint;
  joint.set_gc(gc);
  joint.jointKinematics();
  
  // Step 1. get center of mass of each body
  // << Base >>
  
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
  
  // << LF >>
  Eigen::VectorXd lf_joint_ori, lf_joint_pos, lf_com_ori, lf_com_pos, lf_mass, lf_inertia;
  int lf_child_num = 9;
  lf_joint_ori.setZero(lf_child_num*3); lf_joint_pos.setZero(lf_child_num*3);
  lf_com_ori.setZero(lf_child_num*3); lf_com_pos.setZero(lf_child_num*3);
  lf_mass.setZero(lf_child_num); lf_inertia.setZero(lf_child_num*6);
  lf_com_ori << 0,0,0, 0,0,0, 0,0,0,
    0,0,0, 0,0,0, 0,0,0,
    0,0,0, 0,0,0, 0,0,0;
  lf_com_pos << 0,0,0, 0.048,0.008,-0.003, -0.063,7e-05,0.00046,
    0,0,0, 0.0,0.018,-0.169, -0.063,7e-05,0.00046,
    0,0,0, 0.03463,0.00688,0.00098, 0.00948,-0.00948,0.1468;
  lf_mass << 0.001, 0.74, 2.04, 0.001, 1.03, 2.04, 0.001, 0.33742, 0.25;
  lf_inertia << 0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001,
    0.001393106, 8.4012e-05, 2.3378e-05, 0.003798579, 7.1319e-05, 0.003897509,
    0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827,
    0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001,
    0.018644469, 5.2e-08, 1.0157e-05, 0.019312599, 0.002520077, 0.002838361,
    0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827,
    0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001,
    0.00032748005, 2.142561e-05, 1.33942e-05, 0.00110974122, 7.601e-08, 0.00089388521,
    0.00317174097, 2.63048e-06, 6.815581e-05, 0.00317174092, 6.815583e-05, 8.319196e-05;
  
  std::vector<Eigen::Matrix3d> lf_joint_rot_w(joint.joint_rot_w_.begin() + 1, joint.joint_rot_w_.begin() + lf_child_num + 1);
  std::vector<Eigen::Vector3d> lf_joint_pos_w(joint.joint_pos_w_.begin() + 1, joint.joint_pos_w_.begin() + lf_child_num + 1);
  body.getLegCompositeInertia(lf_joint_rot_w,lf_joint_pos_w,
                              lf_com_ori, lf_com_pos, lf_mass, lf_inertia);
  
  // << RF >>
  Eigen::VectorXd rf_joint_ori, rf_joint_pos, rf_com_ori, rf_com_pos, rf_mass, rf_inertia;
  int rf_child_num = 9;
  rf_joint_ori.setZero(rf_child_num*3); rf_joint_pos.setZero(rf_child_num*3);
  rf_com_ori.setZero(rf_child_num*3); rf_com_pos.setZero(rf_child_num*3);
  rf_mass.setZero(rf_child_num); rf_inertia.setZero(rf_child_num*6);
  rf_com_ori << 0,0,0, 0,0,0, 0,0,0,
    0,0,0, 0,0,0, 0,0,0,
    0,0,0, 0,0,0, 0,0,0;
  rf_com_pos << 0,0,0, 0.048,-0.008,-0.003, -0.063,7e-05,0.00046,
    0,0,0, 0.0,-0.018,-0.169, -0.063,7e-05,0.00046,
    0,0,0, 0.03463,-0.00688,0.00098, 0.00948,0.00948,0.1468;
  rf_mass << 0.001, 0.74, 2.04, 0.001, 1.03, 2.04, 0.001, 0.33742, 0.25;
  rf_inertia << 0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001,
    0.001393106, -8.4012e-05, 2.3378e-05, 0.003798579, -7.1319e-05, 0.003897509,
    0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827,
    0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001,
    0.018644469, -5.2e-08, 1.0157e-05, 0.019312599, -0.002520077, 0.002838361,
    0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827,
    0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001,
    0.00032748005, -2.142561e-05, 1.33942e-05, 0.00110974122, -7.601e-08, 0.00089388521,
    0.00317174097, -2.63048e-06, 6.815581e-05, 0.00317174092, -6.815583e-05, 8.319196e-05;
  
  std::vector<Eigen::Matrix3d> rf_joint_rot_w(joint.joint_rot_w_.begin() + 9 + 2, joint.joint_rot_w_.begin() + 9 + rf_child_num + 2);
  std::vector<Eigen::Vector3d> rf_joint_pos_w(joint.joint_pos_w_.begin() + 9 + 2, joint.joint_pos_w_.begin() + 9 + rf_child_num + 2);
  body.getLegCompositeInertia(rf_joint_rot_w,rf_joint_pos_w,
                              rf_com_ori, rf_com_pos, rf_mass, rf_inertia);
  
  // << LH >>
  Eigen::VectorXd lh_joint_ori, lh_joint_pos, lh_com_ori, lh_com_pos, lh_mass, lh_inertia;
  int lh_child_num = 9;
  lh_joint_ori.setZero(lh_child_num*3); lh_joint_pos.setZero(lh_child_num*3);
  lh_com_ori.setZero(lh_child_num*3); lh_com_pos.setZero(lh_child_num*3);
  lh_mass.setZero(lh_child_num); lh_inertia.setZero(lh_child_num*6);
  lh_com_ori << 0,0,0, 0,0,0, 0,0,0,
    0,0,0, 0,0,0, 0,0,0,
    0,0,0, 0,0,0, 0,0,0;
  lh_com_pos << 0,0,0, -0.048,0.008,-0.003, -0.063,7e-05,0.00046,
    0,0,0, 0.0,0.018,-0.169, -0.063,7e-05,0.00046,
    0,0,0, -0.03463,0.00688,0.00098, -0.00948,-0.00948,0.1468;
  lh_mass << 0.001, 0.74, 2.04, 0.001, 1.03, 2.04, 0.001, 0.33742, 0.25;
  lh_inertia << 0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001,
    0.001393106, -8.4012e-05, -2.3378e-05, 0.003798579, 7.1319e-05, 0.003897509,
    0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827,
    0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001,
    0.018644469, -5.2e-08, -1.0157e-05, 0.019312599, 0.002520077, 0.002838361,
    0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827,
    0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001,
    0.00032748005, -2.142561e-05, -1.33942e-05, 0.00110974122, 7.601e-08, 0.00089388521,
    0.00317174097, -2.63048e-06, -6.815581e-05, 0.00317174092, 6.815583e-05, 8.319196e-05;
  
  std::vector<Eigen::Matrix3d> lh_joint_rot_w(joint.joint_rot_w_.begin() + 9*2 + 3, joint.joint_rot_w_.begin() + 9*2 + lh_child_num + 3);
  std::vector<Eigen::Vector3d> lh_joint_pos_w(joint.joint_pos_w_.begin() + 9*2 + 3, joint.joint_pos_w_.begin() + 9*2 + lh_child_num + 3);
  body.getLegCompositeInertia(lh_joint_rot_w,lh_joint_pos_w,
                              lh_com_ori, lh_com_pos, lh_mass, lh_inertia);
  
  // << RH >>
  Eigen::VectorXd rh_joint_ori, rh_joint_pos, rh_com_ori, rh_com_pos, rh_mass, rh_inertia;
  int rh_child_num = 9;
  rh_joint_ori.setZero(rh_child_num*3); rh_joint_pos.setZero(rh_child_num*3);
  rh_com_ori.setZero(rh_child_num*3); rh_com_pos.setZero(rh_child_num*3);
  rh_mass.setZero(rh_child_num); rh_inertia.setZero(rh_child_num*6);
  rh_com_ori << 0,0,0, 0,0,0, 0,0,0,
    0,0,0, 0,0,0, 0,0,0,
    0,0,0, 0,0,0, 0,0,0;
  rh_com_pos << 0,0,0, -0.048,-0.008,-0.003, -0.063,7e-05,0.00046,
    0,0,0, 0.0,-0.018,-0.169, -0.063,7e-05,0.00046,
    0,0,0, -0.03463,-0.00688,0.00098, -0.00948,0.00948,0.1468;
  rh_mass << 0.001, 0.74, 2.04, 0.001, 1.03, 2.04, 0.001, 0.33742, 0.25;
  rh_inertia << 0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001,
    0.001393106, 8.4012e-05, -2.3378e-05, 0.003798579, -7.1319e-05, 0.003897509,
    0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827,
    0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001,
    0.018644469, 5.2e-08, -1.0157e-05, 0.019312599, -0.002520077, 0.002838361,
    0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827,
    0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001,
    0.00032748005, 2.142561e-05, -1.33942e-05, 0.00110974122, -7.601e-08, 0.00089388521,
    0.00317174097, 2.63048e-06, -6.815581e-05, 0.00317174092, -6.815583e-05, 8.319196e-05;
  
  std::vector<Eigen::Matrix3d> rh_joint_rot_w(joint.joint_rot_w_.begin() + 9*3 + 4, joint.joint_rot_w_.begin() + 9*3 + rh_child_num + 4);
  std::vector<Eigen::Vector3d> rh_joint_pos_w(joint.joint_pos_w_.begin() + 9*3 + 4, joint.joint_pos_w_.begin() + 9*3 + rh_child_num + 4);
  body.getLegCompositeInertia(rh_joint_rot_w,rh_joint_pos_w,
                              rh_com_ori, rh_com_pos, rh_mass, rh_inertia);
  
  
  return Eigen::MatrixXd::Ones(18,18);
}