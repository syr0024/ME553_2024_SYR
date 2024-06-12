//
// Created by jemin on 6/15/22.
//

#pragma once

Eigen::Matrix3d skew (Eigen::Vector3d v) {
  Eigen::Matrix3d S;
  S <<  0, -v.z(),  v.y(),
    v.z(),  0, -v.x(),
    -v.y(),  v.x(), 0;
  return S;
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
  
  return Rx*Ry*Rz;
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

class Joint{
public:
  Joint(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv);
  void set_gc(const Eigen::VectorXd& gc); // gc가 update 될 때마다 불러줘야함
  void set_gv(const Eigen::VectorXd& gv);
  void jointKinematics();
  void getWorldAngularVel();
  void getWorldLinearVel();
  void getSmatrix();
  void getWorldAcc();
  
  // change urdf ! set jacobian size
  std::vector<Eigen::MatrixXd> pos_jacobian_;
  std::vector<Eigen::MatrixXd> vel_jacobian_;
  std::vector<int> joint_leg_idx;
  
  Eigen::VectorXd gc_, gv_;
  Eigen::Vector3d gravity;
  
  std::vector<Eigen::Vector3d> joint_pos_w_; // world frame에서 본 joint position (fixed joint 포함)
  std::vector<Eigen::Matrix3d> joint_rot_w_; // world frame에서 본 joint orientation (3,3 matrix) (fixed joint 포함)
  std::vector<Eigen::Vector3d> angvel_w_;   // world frame에서 본 joint angular velocity (fixed joint 제외, base 포함)
  std::vector<Eigen::Vector3d> linvel_w_;   // world frame에서 본 joint linear velocity (fixed joint 제외, base 포함)
  std::vector<Eigen::MatrixXd> s_;  // motion subspace matrix (fixed joint 제외, base 포함)
  std::vector<Eigen::MatrixXd> s_dot_;  // time derivative of motion subspace matrix (fixed joint 제외, base 포함)
  std::vector<Eigen::Vector3d> a_w_, alpha_w_;  // joint acceleration vector expressed in world frame
  std::vector<Eigen::VectorXd> wrench_; // generalized acceleration vector expressed in world frame (fixed joint 제외, base 포함)
  std::vector<double> joint_sign_;

private:
  static void forwardKinematics(const Eigen::Vector3d& parent_pos_w, const Eigen::Vector3d& parent_ori_w,
                                const Eigen::Vector3d& joint_pos, const Eigen::Matrix3d& joint_ori,
                                Eigen::Vector3d& joint_pos_w, Eigen::Matrix3d& joint_ori_w);
  std::vector<Eigen::Vector3d> joint_ori_, joint_pos_;
  int leg_num_;
  int leg_joint_num_;
};

Joint::Joint(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
  
  /// change urdf !!
  // joint origin, position urdf에 있는 순서로 넣어주기
  // joint sign 넣어주기 (현재 코드는 x축 기준으로 dynamic joint만 고려되어 있음)
  joint_sign_ = {1, 1, 1}; // slider[1,0,0], rod_revolute[0,1,0], rod_revolute2[0,1,0] //TODO: 축이 다 안 동일함
  leg_joint_num_ = 3;
  
  // joint orientation
  Eigen::Vector3d world_to_sliderBar; world_to_sliderBar << 0,0,0;
  Eigen::Vector3d slider; slider << 0,0,0;
  Eigen::Vector3d rod_revolute; rod_revolute << 0,0,0;
  Eigen::Vector3d rod_revolute2; rod_revolute2 << 0,0,0;
  joint_ori_.push_back(world_to_sliderBar);
  joint_ori_.push_back(slider);
  joint_ori_.push_back(rod_revolute);
  joint_ori_.push_back(rod_revolute2);
  
  // joint position
  world_to_sliderBar << 0,0,5.0;
  slider << 0,0,0;
  rod_revolute << 0,0,0;
  rod_revolute2 << 0,0,0;
  joint_pos_.push_back(world_to_sliderBar);
  joint_pos_.push_back(slider);
  joint_pos_.push_back(rod_revolute);
  joint_pos_.push_back(rod_revolute2);
  
  set_gc(gc);
  set_gv(gv);
};

void Joint::set_gc(const Eigen::VectorXd& gc) {
  /// change urdf !
  // gc가 update 될 때마다 불러 줘야 함
  gc_ = gc;

  joint_pos_[1](0) = joint_sign_[0]*gc_(0);
  joint_ori_[2](1) = joint_sign_[1]*gc_(1);
  joint_ori_[3](1) = joint_sign_[2]*gc_(2);
};

void Joint::set_gv(const Eigen::VectorXd& gv) { gv_=gv; };

void Joint::jointKinematics() {
  // joint 값들 초기화
  joint_pos_w_.clear();
  joint_rot_w_.clear();
  
//  // tree 구조 일 경우: 각 joint의 world frame position, orientation 을 update
//  for (int i=0; i < joint_pos_.size(); i++) {
//      Eigen::Vector3d joint_pos_w; Eigen::Matrix3d joint_rot_w;
//      joint_pos_w.setZero(); joint_rot_w.setIdentity();
//      // joint pos, rot 구하는 과정
//      for (int j = i; j>=0; j--) {
//        forwardKinematics(joint_pos_[j], joint_ori_[j],
//                          joint_pos_w, joint_rot_w,
//                          joint_pos_w, joint_rot_w);
//      }
//      // 최종 world frame 기준 joint pos, rot (urdf에 있는 순서로 들어감)
//      joint_pos_w_.push_back(joint_pos_w);
//      joint_rot_w_.push_back(joint_rot_w);
//      std::cout << "pos: " << joint_pos_w.transpose() << std::endl;
//      std::cout << "ori:\n" << joint_rot_w << std::endl;
//  }

// tree 구조 아닐경우
  // world_to_sliderBar (fixed joint)
  joint_pos_w_.push_back(joint_pos_[0]);
  joint_rot_w_.push_back(rpy2rot(joint_ori_[0]));
  
  /// slider
  Eigen::Vector3d joint_pos_w; Eigen::Matrix3d joint_rot_w;
  joint_pos_w.setZero(); joint_rot_w.setIdentity();
  forwardKinematics(joint_pos_[1], joint_ori_[1],
                    joint_pos_w, joint_rot_w,
                    joint_pos_w, joint_rot_w);
  forwardKinematics(joint_pos_[0], joint_ori_[0],
                    joint_pos_w, joint_rot_w,
                    joint_pos_w, joint_rot_w);
  joint_pos_w_.push_back(joint_pos_w);
  joint_rot_w_.push_back(joint_rot_w);
//  std::cout << "pos: " << joint_pos_w.transpose() << std::endl;
//  std::cout << "ori:\n" << joint_rot_w << std::endl;

  /// rod_revolute
  joint_pos_w.setZero(); joint_rot_w.setIdentity();
  forwardKinematics(joint_pos_[2], joint_ori_[2],
                    joint_pos_w, joint_rot_w,
                    joint_pos_w, joint_rot_w);
  forwardKinematics(joint_pos_[1], joint_ori_[1],
                    joint_pos_w, joint_rot_w,
                    joint_pos_w, joint_rot_w);
  forwardKinematics(joint_pos_[0], joint_ori_[0],
                    joint_pos_w, joint_rot_w,
                    joint_pos_w, joint_rot_w);
  joint_pos_w_.push_back(joint_pos_w);
  joint_rot_w_.push_back(joint_rot_w);
//  std::cout << "pos: " << joint_pos_w.transpose() << std::endl;
//  std::cout << "ori:\n" << joint_rot_w << std::endl;
  
  /// rod_revolute2
  joint_pos_w.setZero(); joint_rot_w.setIdentity();
  forwardKinematics(joint_pos_[3], joint_ori_[3],
                    joint_pos_w, joint_rot_w,
                    joint_pos_w, joint_rot_w);
  forwardKinematics(joint_pos_[1], joint_ori_[1],
                    joint_pos_w, joint_rot_w,
                    joint_pos_w, joint_rot_w);
  forwardKinematics(joint_pos_[0], joint_ori_[0],
                    joint_pos_w, joint_rot_w,
                    joint_pos_w, joint_rot_w);
  joint_pos_w_.push_back(joint_pos_w);
  joint_rot_w_.push_back(joint_rot_w);
//  std::cout << "pos: " << joint_pos_w.transpose() << std::endl;
//  std::cout << "ori:\n" << joint_rot_w << std::endl;
  
  getWorldAngularVel();
  getWorldLinearVel();
  getSmatrix();
  getWorldAcc();
}

void Joint::getWorldAngularVel() {
  angvel_w_.clear();
  
  /// slider
  angvel_w_.push_back(Eigen::Vector3d::Zero());
//  std::cout << "angvel_w_: " << angvel_w_.back().transpose() << std::endl;
  
  Eigen::Vector3d joint_axis;
  /// rod_revolute
  joint_axis.setZero();
  joint_axis = joint_sign_[1] * joint_rot_w_[2].block<3,1>(0,1);
  angvel_w_.push_back( angvel_w_[0] + joint_axis * gv_[1] );
//  std::cout << "angvel_w_: " << angvel_w_.back().transpose() << std::endl;
  
  /// rod_revolute2
  joint_axis.setZero();
  joint_axis = joint_sign_[2] * joint_rot_w_[3].block<3,1>(0,1);
  angvel_w_.push_back( angvel_w_[0] + joint_axis * gv_[2] );
//  std::cout << "angvel_w_: " << angvel_w_.back().transpose() << std::endl;
  
}

void Joint::getWorldLinearVel() {
  linvel_w_.clear();
  
  Eigen::MatrixXd positional_jacobain; positional_jacobain.setZero(3, 3);
  /// slider
  positional_jacobain.block<3,1>(0,0) = joint_sign_[0] * joint_rot_w_[1].block<3,1>(0,0);
  linvel_w_.push_back(positional_jacobain * gv_);
//  std::cout << "linvel_w_: " << linvel_w_.back().transpose() << std::endl;
  /// rod_revolute
  positional_jacobain.block<3,1>(0,1) = (joint_sign_[1] * joint_rot_w_[2].block<3,1>(0,1)).cross(joint_pos_w_[2] - joint_pos_w_[1]);
  linvel_w_.push_back(positional_jacobain * gv_);
//  std::cout << "linvel_w_: " << linvel_w_.back().transpose() << std::endl;
  /// rod_revolute2
  positional_jacobain.block<3,1>(0,1) = (joint_sign_[1] * joint_rot_w_[2].block<3,1>(0,1)).cross(joint_pos_w_[2] - joint_pos_w_[1]);
  linvel_w_.push_back(positional_jacobain * gv_);
//  std::cout << "linvel_w_: " << linvel_w_.back().transpose() << std::endl;
}

void Joint::getSmatrix() {
  s_.clear();
  s_dot_.clear();
  
  Eigen::Vector3d joint_axis; 
  Eigen::VectorXd s, s_dot;
  /// silder
  joint_axis.setZero(); s.setZero(6); s_dot.setZero(6);
  joint_axis = joint_sign_[0] * joint_rot_w_[1].block<3,1>(0,0); // joint: axis xyz 고려해줘야함
  s.segment(0,3) = joint_axis;
  s_.push_back(s);
  s_dot.segment(3,3) = skew(angvel_w_[0]) * joint_axis;
  s_dot_.push_back(s_dot);
  /// rod_revolute
  joint_axis.setZero(); s.setZero(6); s_dot.setZero(6);
  joint_axis = joint_sign_[1] * joint_rot_w_[2].block<3,1>(0,1);
  s.segment(3,3) = joint_axis;
  s_.push_back(s);
  s_dot.segment(3,3) = skew(angvel_w_[1]) * joint_axis;
  s_dot_.push_back(s_dot);
  /// rod_revolute2
  joint_axis.setZero(); s.setZero(6); s_dot.setZero(6);
  joint_axis = joint_sign_[2] * joint_rot_w_[3].block<3,1>(0,1);
  s.segment(3,3) = joint_axis;
  s_.push_back(s);
  s_dot.segment(3,3) = skew(angvel_w_[2]) * joint_axis;
  s_dot_.push_back(s_dot);
}

void Joint::getWorldAcc() {
  a_w_.clear();
  alpha_w_.clear();
  wrench_.clear();
  
  Eigen::Vector3d a_w, alpha_w;
  Eigen::VectorXd wrench; wrench.setZero(6);
  /// slider
  a_w.setZero(); alpha_w.setZero();
  wrench << a_w, alpha_w;
  wrench += s_dot_[0] * gv_[0];
  wrench_.push_back(wrench);
  a_w << wrench.segment(0,3);
  alpha_w << wrench.segment(3,3);
  a_w_.push_back(a_w);
  alpha_w_.push_back(alpha_w);
//  std::cout << "a_w_: " << a_w_.back().transpose() << std::endl;
//  std::cout << "alpha_w_: " << alpha_w_.back().transpose() << std::endl;
  /// rod_revolute
  a_w.setZero(); alpha_w.setZero();
  a_w = a_w_[0] + skew(alpha_w_[0]) * (joint_pos_w_[1] - joint_pos_w_[0])
        + skew(angvel_w_[0])*skew(angvel_w_[0]) * (joint_pos_w_[1] - joint_pos_w_[0]);
  alpha_w = alpha_w_[0];
  wrench << a_w, alpha_w;
  wrench += s_dot_[1] * gv_[1];
  wrench_.push_back(wrench);
  a_w << wrench.segment(0,3);
  alpha_w << wrench.segment(3,3);
  a_w_.push_back(a_w);
  alpha_w_.push_back(alpha_w);
//  std::cout << "a_w_: " << a_w_.back().transpose() << std::endl;
//  std::cout << "alpha_w_: " << alpha_w_.back().transpose() << std::endl;
  /// rod_revolute2
  a_w.setZero(); alpha_w.setZero();
  a_w = a_w_[1] + skew(alpha_w_[1]) * (joint_pos_w_[2] - joint_pos_w_[0])
        + skew(angvel_w_[1])*skew(angvel_w_[1]) * (joint_pos_w_[2] - joint_pos_w_[0]);
  alpha_w = alpha_w_[0];
  wrench << a_w, alpha_w;
  wrench += s_dot_[2] * gv_[2];
  wrench_.push_back(wrench);
  a_w << wrench.segment(0,3);
  alpha_w << wrench.segment(3,3);
  a_w_.push_back(a_w);
  alpha_w_.push_back(alpha_w);
//  std::cout << "a_w_: " << a_w_.back().transpose() << std::endl;
//  std::cout << "alpha_w_: " << alpha_w_.back().transpose() << std::endl;
}

void Joint::forwardKinematics(const Eigen::Vector3d& parent_pos_w, const Eigen::Vector3d& parent_ori_w,
                              const Eigen::Vector3d& joint_pos, const Eigen::Matrix3d& joint_ori,
                              Eigen::Vector3d& joint_pos_w, Eigen::Matrix3d& joint_ori_w) {
  joint_pos_w = parent_pos_w + rpy2rot(parent_ori_w)*joint_pos;
  joint_ori_w = rpy2rot(parent_ori_w) * joint_ori;
}

class Body{
public:
  Body(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv);
  void reset(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv);
  Eigen::Matrix3d getInertia (const Eigen::VectorXd& inertia_vec);
  void  getComposite (const Eigen::Vector3d& r1_com, const Eigen::Vector3d& r2_com,
                      const double& m1, const double& m2,
                      const Eigen::Matrix3d& I1, const Eigen::Matrix3d& I2,
                      Eigen::Vector3d& r_com, double& mass, Eigen::Matrix3d& I_new);
  void revComposite (const Eigen::Vector3d& r_com, const Eigen::Vector3d& r2_com,
                     const double& m, const double& m2,
                     const Eigen::Matrix3d& I, const Eigen::Matrix3d& I2,
                     Eigen::Vector3d& r1_com, double& m1, Eigen::Matrix3d& I1);
  Eigen::MatrixXd SpatialInertiaMat (const double& mass, const Eigen::Matrix3d& I_j, const Eigen::Vector3d& r_jc);

public:
  Eigen::VectorXd gc_, gv_;
  
  // composite com, mass, inertia in world frame containing child (child 를 포함한 composite com, mass, inertia)
  std::vector<Eigen::Vector3d> com_COM_child_w_;
  std::vector<double> com_Mass_child_;
  std::vector<Eigen::Matrix3d> com_Inertia_child_w_;
  
  // composite com, mass, inertia in world frame containing child (각 body의 composite com, mass, inertia)
  std::vector<Eigen::Vector3d> com_COM_w_;
  std::vector<double> com_Mass_;
  std::vector<Eigen::Matrix3d> com_Inertia_w_;
  
  /// Calculated by RNE
  // spatial inertia matrix (각각의 body들만 고려한 spatial inertia matrix)
  std::vector<Eigen::MatrixXd> Mc_;
  // nonlinearity b only consider body (각각의 body만 고려해 구한 b)
  std::vector<Eigen::VectorXd> b_;
  /// Calculated by ABA
  // spatial inertia matrix of ABA
  std::vector<Eigen::MatrixXd> M_A_; // Articulated body inertia
  std::vector<Eigen::VectorXd> b_A_; // Articulated body bias
  
};

Body::Body(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) { gc_ = gc; gv_ = gv;};

void Body::reset(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
  gc_ = gc;
  gv_ = gv;
  /// 이전에 구해놓은 값들을 초기화, gc, gv가 바뀌었을 경우 다시 구해줘야함
  com_COM_child_w_.clear();
  com_Mass_child_.clear();
  com_Inertia_child_w_.clear();
  com_COM_w_.clear();
  com_Mass_.clear();
  com_Inertia_w_.clear();
  Mc_.clear();
  b_.clear();
}

void Body::getComposite (const Eigen::Vector3d& r1_com, const Eigen::Vector3d& r2_com,
                         const double& m1, const double& m2,
                         const Eigen::Matrix3d& I1, const Eigen::Matrix3d& I2,
                         Eigen::Vector3d& r_com, double& mass, Eigen::Matrix3d& I_new) {
  /**
   2개의 rigid body가 주어졌을 때, composite com, mass, inertia를 구해주는 함수
   **/
  
  Eigen::Vector3d r_com_new = (m1*r1_com + m2*r2_com)/(m1+m2);
  Eigen::Vector3d r1 = r1_com - r_com_new;
  Eigen::Vector3d r2 = r2_com - r_com_new;
  
  I_new = I1 + I2 - m1*skew(r1)*skew(r1) - m2*skew(r2)*skew(r2);
  mass = m1 + m2;
  r_com = r_com_new;
}

void Body::revComposite (const Eigen::Vector3d& r_com, const Eigen::Vector3d& r2_com,
                         const double& m, const double& m2,
                         const Eigen::Matrix3d& I, const Eigen::Matrix3d& I2,
                         Eigen::Vector3d& r1_com, double& m1, Eigen::Matrix3d& I1) {
  // reverse composite
  m1 = m - m2;
  r1_com = (m*r_com - m2*r2_com)/m1;
  Eigen::Vector3d r1 = r1_com - r_com;
  Eigen::Vector3d r2 = r2_com - r_com;
  I1 = I - I2 + m1*skew(r1)*skew(r1) + m2*skew(r2)*skew(r2);
}

Eigen::MatrixXd Body::SpatialInertiaMat (const double& mass, const Eigen::Matrix3d& I_j, const Eigen::Vector3d& r_jc) {
  
  Eigen::Matrix<double, 6, 6> Mc;
  Mc.block<3,3>(0,0) = mass*Eigen::Matrix3d::Identity();
  Mc.block<3,3>(3,0) = mass*skew(r_jc);
  Mc.block<3,3>(0,3) = - mass*skew(r_jc);
  Mc.block<3,3>(3,3) = I_j - mass*skew(r_jc)*skew(r_jc);
  
  return Mc;
}

Eigen::Matrix3d Body::getInertia (const Eigen::VectorXd& inertia_vec) {
  /// urdf 에서 6개의 variable로 주어진 inertia를 Matrix3d로 변환해주는 함수
  Eigen::Matrix3d inertia;
  inertia << inertia_vec[0], inertia_vec[1], inertia_vec[2],
    inertia_vec[1], inertia_vec[3], inertia_vec[4],
    inertia_vec[2], inertia_vec[4], inertia_vec[5];
  return inertia;
}

/// do not change the name of the method
inline Eigen::VectorXd getNonlinearities (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
  
  // gravity
  Eigen::Vector3d gravity;
  gravity << 0, 0, -9.81;
  
  Joint joint(gc, gv);
  joint.gravity << gravity;
  joint.jointKinematics();
  
  Body body(gc, gv);
  body.reset(gc, gv);
  
  int idx;
  
  /// CRBA ///
  /// Step 1. get center of mass of each body
  Eigen::Vector3d com_pos_w;
  Eigen::Matrix3d com_rot_w, inertia_w;
  Eigen::Vector3d com_new;
  double mass_new;
  Eigen::Matrix3d inertia_new;
  
  /// rod2
  Eigen::Vector3d rod2_com_ori; rod2_com_ori << 0,0,0;
  Eigen::Vector3d rod2_com_pos; rod2_com_pos << 0,0,0.5;
  Eigen::VectorXd rod2_inertia; rod2_inertia.setZero(6); rod2_inertia << 1.0,0.0,0.0,1.0,0.0,1.0;
  double rod2_mass; rod2_mass = 5;
  com_pos_w = joint.joint_pos_w_[3] + joint.joint_rot_w_[3] * rod2_com_pos;
  com_rot_w = joint.joint_rot_w_[3] * rpy2rot(rod2_com_ori);
  inertia_w = com_rot_w * body.getInertia(rod2_inertia) * com_rot_w.transpose();
  // child 제외한 body만의 COM, Mass, Inertia
  body.com_COM_w_.push_back(com_pos_w);
  body.com_Mass_.push_back(rod2_mass);
  body.com_Inertia_w_.push_back(inertia_w);
  // child 포함한 COM, Mass, Inertia : leaf body 이므로 추가만 해줌
  body.com_COM_child_w_.push_back(com_pos_w);
  body.com_Mass_child_.push_back(rod2_mass);
  body.com_Inertia_child_w_.push_back(inertia_w);
//  std::cout << "com_test: " << body.com_COM_child_w_.back().transpose() << std::endl;
//  std::cout << body.com_Mass_child_.back() << std::endl;
//  std::cout << body.com_Inertia_child_w_.back() << std::endl;
  
  /// rod
  Eigen::Vector3d rod_com_ori; rod_com_ori << 0,0,0;
  Eigen::Vector3d rod_com_pos; rod_com_pos << 0,0,0.5;
  Eigen::VectorXd rod_inertia; rod_inertia.setZero(6); rod_inertia << 1.0,0.0,0.0,1.0,0.0,1.0;
  double rod_mass; rod_mass = 5;
  com_pos_w = joint.joint_pos_w_[2] + joint.joint_rot_w_[2] * rod2_com_pos;
  com_rot_w = joint.joint_rot_w_[2] * rpy2rot(rod2_com_ori);
  inertia_w = com_rot_w * body.getInertia(rod2_inertia) * com_rot_w.transpose();
  // child 제외한 body만의 COM, Mass, Inertia
  body.com_COM_w_.push_back(com_pos_w);
  body.com_Mass_.push_back(rod_mass);
  body.com_Inertia_w_.push_back(inertia_w);
  // child 포함한 COM, Mass, Inertia : leaf body 이므로 추가만 해줌
  body.com_COM_child_w_.push_back(com_pos_w);
  body.com_Mass_child_.push_back(rod_mass);
  body.com_Inertia_child_w_.push_back(inertia_w);
//  std::cout << "com_test: "<< body.com_COM_child_w_.back().transpose() << std::endl;
//  std::cout << body.com_Mass_child_.back() << std::endl;
//  std::cout << body.com_Inertia_child_w_.back() << std::endl;
  
  /// slider
  Eigen::Vector3d slider_com_ori; slider_com_ori << 0,0,0;
  Eigen::Vector3d slider_com_pos; slider_com_pos << 0,0,0;
  Eigen::VectorXd slider_inertia; slider_inertia.setZero(6); slider_inertia << 2.0,0.0,0.0,1.0,0.0,2.0;
  double slider_mass; slider_mass = 2;
  com_pos_w = joint.joint_pos_w_[1] + joint.joint_rot_w_[1] * slider_com_pos;
  com_rot_w = joint.joint_rot_w_[1] * rpy2rot(slider_com_ori);
  inertia_w = com_rot_w * body.getInertia(slider_inertia) * com_rot_w.transpose();
  // child 제외한 body만의 COM, Mass, Inertia
  body.com_COM_w_.push_back(com_pos_w);
  body.com_Mass_.push_back(slider_mass);
  body.com_Inertia_w_.push_back(inertia_w);
  // child 포함한 COM, Mass, Inertia
  body.getComposite(body.com_COM_child_w_[0], com_pos_w, body.com_Mass_child_[0], slider_mass, body.com_Inertia_child_w_[0], inertia_w,
                    com_new, mass_new, inertia_new);
  body.getComposite(body.com_COM_child_w_[1], com_new, body.com_Mass_child_[1], mass_new, body.com_Inertia_child_w_[1], inertia_new,
                    com_new, mass_new, inertia_new);
  body.com_COM_child_w_.push_back(com_new);
  body.com_Mass_child_.push_back(mass_new);
  body.com_Inertia_child_w_.push_back(inertia_new);
//  std::cout << "com_test: " << body.com_COM_child_w_.back().transpose() << std::endl;
//  std::cout << body.com_Mass_child_.back() << std::endl;
//  std::cout << body.com_Inertia_child_w_.back() << std::endl;
  
  std::reverse(body.com_COM_w_.begin(), body.com_COM_w_.end());
  std::reverse(body.com_Mass_.begin(), body.com_Mass_.end());
  std::reverse(body.com_Inertia_w_.begin(), body.com_Inertia_w_.end());

  std::reverse(body.com_COM_child_w_.begin(), body.com_COM_child_w_.end());
  std::reverse(body.com_Mass_child_.begin(), body.com_Mass_child_.end());
  std::reverse(body.com_Inertia_child_w_.begin(), body.com_Inertia_child_w_.end());
  
  /// Step 2. get Mass matrix
  // com_COM_child_w_, com_Mass_child_, com_Inertia_child_w_: j 기준
  // joint_*_w_: fixed joint 포함
//  std::cout << "joint_pos_w_ size: " << joint.joint_pos_w_.size() << std::endl;
  // com_*_w_: fixed joint 제외
//  std::cout << "com_COM_child_w_ size: " << body.com_COM_child_w_.size() << std::endl;
  
  Eigen::Matrix<double, 3, 3> Mass_mat; Mass_mat.setZero();
  Eigen::Vector3d r_jc, r_ij;
  Eigen::Matrix<double, 6, 6> Mc;
  Eigen::Matrix<double,6,6> A;
  
  // [i,j]=[2,2]
  r_jc = body.com_COM_child_w_[2] - joint.joint_pos_w_[3];
  r_ij = joint.joint_pos_w_[3] - joint.joint_pos_w_[3];
  Mc = body.SpatialInertiaMat(body.com_Mass_child_[2], body.com_Inertia_child_w_[2], r_jc);
  A.setZero();
  A.block<3,3>(0,0) = Eigen::Matrix3d::Identity();
  A.block<3,3>(3,3) = Eigen::Matrix3d::Identity();
  A.block<3,3>(0,3) = - skew(r_ij);
  Mass_mat(2,2) = (joint.s_[2].transpose() * Mc * A * joint.s_[2])(0,0);
  
  // [i,j]=[1,2] = 0
  
  // [i,j]=[0,2]
  r_jc = body.com_COM_child_w_[2] - joint.joint_pos_w_[3];
  r_ij = joint.joint_pos_w_[3] - joint.joint_pos_w_[1];
  Mc = body.SpatialInertiaMat(body.com_Mass_child_[2], body.com_Inertia_child_w_[2], r_jc);
  A.setZero();
  A.block<3,3>(0,0) = Eigen::Matrix3d::Identity();
  A.block<3,3>(3,3) = Eigen::Matrix3d::Identity();
  A.block<3,3>(0,3) = - skew(r_ij);
  Mass_mat(0,2) = (joint.s_[2].transpose() * Mc * A * joint.s_[0])(0,0); // s_j.transpose() * Mc * A * s_i
  
  // [i,j]=[1,1]
  r_jc = body.com_COM_child_w_[1] - joint.joint_pos_w_[2];
  r_ij = joint.joint_pos_w_[2] - joint.joint_pos_w_[2];
  Mc = body.SpatialInertiaMat(body.com_Mass_child_[1], body.com_Inertia_child_w_[1], r_jc);
  A.setZero();
  A.block<3,3>(0,0) = Eigen::Matrix3d::Identity();
  A.block<3,3>(3,3) = Eigen::Matrix3d::Identity();
  A.block<3,3>(0,3) = - skew(r_ij);
  Mass_mat(1,1) = (joint.s_[1].transpose() * Mc * A * joint.s_[1])(0,0); // s_j.transpose() * Mc * A * s_i
  
  // [i,j]=[0,1]
  r_jc = body.com_COM_child_w_[1] - joint.joint_pos_w_[1];
  r_ij = joint.joint_pos_w_[1] - joint.joint_pos_w_[1];
  Mc = body.SpatialInertiaMat(body.com_Mass_child_[1], body.com_Inertia_child_w_[1], r_jc);
  A.setZero();
  A.block<3,3>(0,0) = Eigen::Matrix3d::Identity();
  A.block<3,3>(3,3) = Eigen::Matrix3d::Identity();
  A.block<3,3>(0,3) = - skew(r_ij);
  Mass_mat(0,1) = (joint.s_[1].transpose() * Mc * A * joint.s_[0])(0,0); // s_j.transpose() * Mc * A * s_i
  
  // [0,0]
  r_jc = body.com_COM_child_w_[0] - joint.joint_pos_w_[1];
  r_ij = joint.joint_pos_w_[1] - joint.joint_pos_w_[1];
  Mc = body.SpatialInertiaMat(body.com_Mass_child_[0], body.com_Inertia_child_w_[0], r_jc);
  A.setZero();
  A.block<3,3>(0,0) = Eigen::Matrix3d::Identity();
  A.block<3,3>(3,3) = Eigen::Matrix3d::Identity();
  A.block<3,3>(0,3) = - skew(r_ij);
  Mass_mat(0,0) = (joint.s_[0].transpose() * Mc * A * joint.s_[0])(0,0); // s_j.transpose() * Mc * A * s_i
  
  std::cout << "my Mass matrix is\n" << Mass_mat << std::endl;
  
  /// RNE ///
  std::vector<Eigen::Vector3d> F_w_;
  std::vector<Eigen::Vector3d> tau_w_;
  Eigen::VectorXd base_gen_F; base_gen_F.setZero(6);
  
  /// Step 1. Compute body force (leaves to the root)
  Eigen::Matrix<double,6,1> gen_F;
  
  /// rod2
  idx = 2;
  // get Spatial Inertia Matrix
  r_jc = body.com_COM_w_[idx] - joint.joint_pos_w_[idx+1];
  Mc = body.SpatialInertiaMat(body.com_Mass_[idx], body.com_Inertia_w_[idx], r_jc);
  // get Generalized Force
  gen_F.segment(0,3) = body.com_Mass_[idx]*skew(joint.angvel_w_[idx])*skew(joint.angvel_w_[idx])*r_jc;
  gen_F.segment(3,3) = skew(joint.angvel_w_[idx]) * (body.com_Inertia_w_[idx] - body.com_Mass_[idx]*skew(r_jc)*skew(r_jc)) * joint.angvel_w_[idx];
  // + force applied from the child : leaf body -> skip
  // + M*[a alpha]
  gen_F += Mc * joint.wrench_[idx];
  // gravity
  gen_F.segment(0,3) -= body.com_Mass_[idx]*gravity;
  gen_F.segment(3,3) -= body.com_COM_w_[idx].cross(body.com_Mass_[idx]*gravity);
  
  F_w_.push_back(gen_F.segment(0,3));
  tau_w_.push_back(gen_F.segment(3,3));
  
  /// rod1
  idx = 1;
  // get Spatial Inertia Matrix
  r_jc = body.com_COM_w_[idx] - joint.joint_pos_w_[idx+1];
  Mc = body.SpatialInertiaMat(body.com_Mass_[idx], body.com_Inertia_w_[idx], r_jc);
  // get Generalized Force
  gen_F.segment(0,3) = body.com_Mass_[idx]*skew(joint.angvel_w_[idx])*skew(joint.angvel_w_[idx])*r_jc;
  gen_F.segment(3,3) = skew(joint.angvel_w_[idx]) * (body.com_Inertia_w_[idx] - body.com_Mass_[idx]*skew(r_jc)*skew(r_jc)) * joint.angvel_w_[idx];
  // + force applied from the child : leaf body -> skip
  // + M*[a alpha]
  gen_F += Mc * joint.wrench_[idx];
  // gravity
  gen_F.segment(0,3) -= body.com_Mass_[idx]*gravity;
  gen_F.segment(3,3) -= body.com_COM_w_[idx].cross(body.com_Mass_[idx]*gravity);
  
  F_w_.push_back(gen_F.segment(0,3));
  tau_w_.push_back(gen_F.segment(3,3));
  
  /// slider
  idx = 0;
  // get Spatial Inertia Matrix
  r_jc = body.com_COM_w_[idx] - joint.joint_pos_w_[idx+1];
  Mc = body.SpatialInertiaMat(body.com_Mass_[idx], body.com_Inertia_w_[idx], r_jc);
  // get Generalized Force
  gen_F.segment(0,3) = body.com_Mass_[idx]*skew(joint.angvel_w_[idx])*skew(joint.angvel_w_[idx])*r_jc;
  gen_F.segment(3,3) = skew(joint.angvel_w_[idx]) * (body.com_Inertia_w_[idx] - body.com_Mass_[idx]*skew(r_jc)*skew(r_jc)) * joint.angvel_w_[idx];
  // + force applied from the child
  gen_F.segment(0,3) += F_w_[0];
  gen_F.segment(0,3) += F_w_[1];
  gen_F.segment(3,3) += tau_w_[0] + skew(joint.joint_pos_w_[3] - joint.joint_pos_w_[idx])*F_w_[0]; // index 주의! F_w_는 뒤집혀 있음, skew(joint_pos_child - joint_pos_parent)
  gen_F.segment(3,3) += tau_w_[1] + skew(joint.joint_pos_w_[2] - joint.joint_pos_w_[idx])*F_w_[1];
  // + M*[a alpha]
  gen_F += Mc * joint.wrench_[idx];
  F_w_.push_back(gen_F.segment(0,3));
  tau_w_.push_back(gen_F.segment(3,3));
  
  std::reverse(F_w_.begin(), F_w_.end());
  std::reverse(tau_w_.begin(), tau_w_.end());
  
  Eigen::VectorXd b; b.setZero(3);
  /// Step 2. get Nonlinearities (root to the leaves)
  for (int i=0; i<3; i++) {
    gen_F.setZero(6);
    
    gen_F << F_w_[i], tau_w_[i];
    b.segment(i,1) = joint.s_[i].transpose() * gen_F;
  }
  std::cout << "my nonlinearity answer is : " << b.transpose() << std::endl;
  
  
  /// ABA ///
  

  return Eigen::VectorXd::Ones(3);
}