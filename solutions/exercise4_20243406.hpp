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


class Joint{
public:
  Joint(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv);
  void set_gc(const Eigen::VectorXd& gc); // gc가 update 될 때마다 불러줘야함
  void set_gv(const Eigen::VectorXd& gv);
  void jointKinematics();
  void getWorldAngularVel();
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
  std::vector<Eigen::Vector3d> angvel_w_;   // world frame에서 본 joint angular velocity (fixed joint 제외)
  std::vector<Eigen::MatrixXd> s_;  // motion subspace matrix (fixed joint 제외)
  std::vector<Eigen::MatrixXd> s_dot_;  // time derivative of motion subspace matrix (fixed joint 제외)
  std::vector<Eigen::Vector3d> a_w_, alpha_w_;  // joint acceleration vector expressed in world frame
  std::vector<Eigen::VectorXd> gen_a_; // generalized acceleration vector expressed in world frame
  std::vector<double> joint_sign_;

private:
  static void forwardKinematics(const Eigen::Vector3d& parent_pos_w, const Eigen::Vector3d& parent_ori_w,
                                const Eigen::Vector3d& joint_pos, const Eigen::Matrix3d& joint_ori,
                                Eigen::Vector3d& joint_pos_w, Eigen::Matrix3d& joint_ori_w);
  std::vector<Eigen::VectorXd> joint_ori_, joint_pos_;
  int leg_num_;
  int leg_joint_num_;
};

Joint::Joint(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
  
  /// change urdf !!
  // joint origin, position urdf에 있는 순서로 넣어주기
  // joint sign 넣어주기 (현재 코드는 x축 기준으로 revolute joint만 고려되어 있음)
  joint_sign_ = {1,1,1, 1,-1,-1, -1,1,1, -1,-1,-1}; // lf(haa, hfe, kfe), rf, lh, rh
  leg_num_ = 4;
  leg_joint_num_ = 3;
  joint_leg_idx = {1,4,7, 11,14,17, 21,24,27, 31,34,37}; // lf(haa, hfe, kfe), rf, lh, rh
  
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
  
  set_gc(gc);
  set_gv(gv);
};

void Joint::set_gc(const Eigen::VectorXd& gc) {
  /// change urdf !
  // gc가 update 될 때마다 불러 줘야 함
  gc_ = gc;
  // lf
  joint_ori_[0].segment(3,1) = joint_sign_[0]*gc_.segment(7,1);
  joint_ori_[0].segment(12,1) = joint_sign_[1]*gc_.segment(8,1);
  joint_ori_[0].segment(21,1) = joint_sign_[2]*gc_.segment(9,1);
  // rf
  joint_ori_[1].segment(3,1) = joint_sign_[3]*gc_.segment(10,1);
  joint_ori_[1].segment(12,1) = joint_sign_[4]*gc_.segment(11,1);
  joint_ori_[1].segment(21,1) = joint_sign_[5]*gc_.segment(12,1);
  // lh
  joint_ori_[2].segment(3,1) = joint_sign_[6]*gc_.segment(13,1);
  joint_ori_[2].segment(12,1) = joint_sign_[7]*gc_.segment(14,1);
  joint_ori_[2].segment(21,1) = joint_sign_[8]*gc_.segment(15,1);
  // rh
  joint_ori_[3].segment(3,1) = joint_sign_[9]*gc_.segment(16,1);
  joint_ori_[3].segment(12,1) = joint_sign_[10]*gc_.segment(17,1);
  joint_ori_[3].segment(21,1) = joint_sign_[11]*gc_.segment(18,1);
};

void Joint::set_gv(const Eigen::VectorXd& gv) { gv_=gv; };

void Joint::jointKinematics() {
  // joint 값들 초기화
  joint_pos_w_.clear();
  joint_rot_w_.clear();
  
  // 각 joint의 world frame position, orientation 을 update
  for (int k=0; k < joint_pos_.size(); k++) {
    for (int i=0; i < joint_pos_[k].size(); i+=3) { // i: tree 전체의 joint 개수
      Eigen::Vector3d joint_pos_w; Eigen::Matrix3d joint_rot_w;
      joint_pos_w.setZero(); joint_rot_w.setIdentity();
      // joint pos, rot 구하는 과정
      for (int j = i + 3; j>=3; j-=3) {
        forwardKinematics(joint_pos_[k].segment(j-3,3), joint_ori_[k].segment(j-3,3),
                          joint_pos_w, joint_rot_w,
                          joint_pos_w, joint_rot_w);
      }
      // 최종 world frame 기준 joint pos, rot (urdf에 있는 순서로 들어감)
      joint_pos_w = gc_.head(3) + quat2rot(gc_.segment(3,4))*joint_pos_w;
      joint_rot_w = quat2rot(gc_.segment(3,4)) * joint_rot_w;
      joint_pos_w_.push_back(joint_pos_w);
      joint_rot_w_.push_back(joint_rot_w);
    }
  }
  
  getWorldAngularVel();
  getSmatrix();
  getWorldAcc();
}

void Joint::getWorldAngularVel() {
  angvel_w_.clear();
  angvel_w_.push_back(gv_.segment(3,3));
  
  for (int k=0; k<leg_num_; k++) {
    for (int i=0; i<leg_joint_num_; i++) {
      int idx = k*leg_joint_num_+i;
      Eigen::Vector3d joint_axis = joint_sign_[idx] * joint_rot_w_[joint_leg_idx[idx]].block<3,1>(0,0);

      if (i==0) { angvel_w_.push_back( angvel_w_[0] + joint_axis * gv_[6+idx] ); }
      else { angvel_w_.push_back( angvel_w_.back() + joint_axis * gv_[6+idx] ); }
    }
  }
}

void Joint::getSmatrix() { //TODO : check
  s_.clear();
  s_dot_.clear();
  
  // base
  Eigen::Matrix<double,6,6> s0; s0.setZero();
  s0.block<3,3>(0,0) << Eigen::Matrix3d::Identity();
  s0.block<3,3>(3,3) << Eigen::Matrix3d::Identity();
  s_.push_back(s0);
  Eigen::Matrix<double,6,6> s0_dot; s0.setZero();
  s_dot_.push_back(s0_dot);
  
  // tree
  for (int k=0; k<leg_num_; k++) {
    for (int i=0; i<leg_joint_num_; i++) {
      int idx = k*leg_joint_num_+i;
      Eigen::Vector3d joint_axis = joint_sign_[idx] * joint_rot_w_[joint_leg_idx[idx]].block<3,1>(0,0);
      // get s matrix
      Eigen::VectorXd s; s.setZero(6);
      s.segment(3,3) = joint_axis;
      s_.push_back(s);
      // get time derivative s matrix
      Eigen::VectorXd s_dot; s_dot.setZero(6);
      s_dot.segment(3,3) = skew(angvel_w_[1+idx]) * joint_axis;
      s_dot_.push_back(s_dot);
    }
  }
}

void Joint::getWorldAcc() {
  a_w_.clear();
  alpha_w_.clear();
  gen_a_.clear();
  
  // base
  Eigen::Vector3d a_w, alpha_w;
  a_w.setZero(); alpha_w.setZero();
  Eigen::VectorXd gen_a; gen_a.setZero(6);
  a_w_.push_back(gravity);
  alpha_w_.push_back(alpha_w);
  
  gen_a << a_w_.back(), alpha_w_.back();
  gen_a_.push_back(gen_a);
  
  // tree
  for (int k=0; k<leg_num_; k++) {
    for (int i=0; i<leg_joint_num_; i++) {
      int idx = k*leg_joint_num_+i;
      
      if (i==0) {
        a_w = a_w_[0] + skew(alpha_w_[0]) * (joint_pos_w_[joint_leg_idx[idx]] - gc_.head(3))
              + skew(angvel_w_[0])*skew(angvel_w_[0]) * (joint_pos_w_[joint_leg_idx[idx]] - gc_.head(3));
        alpha_w = alpha_w_[0];
      }
      else {
        a_w = a_w_.back() + skew(alpha_w_.back()) * (joint_pos_w_[joint_leg_idx[idx]] - joint_pos_w_[joint_leg_idx[idx-1]])
              + skew(angvel_w_[idx])*skew(angvel_w_[idx]) * (joint_pos_w_[joint_leg_idx[idx]] - joint_pos_w_[joint_leg_idx[idx-1]]);
        alpha_w = alpha_w_.back();
      }
      
      gen_a << a_w, alpha_w;
      gen_a += s_dot_[1+idx] * gv_(6+idx);
      gen_a_.push_back(gen_a);
      
      a_w << gen_a.segment(0,3);
      alpha_w << gen_a.segment(3,3);
      a_w_.push_back(a_w);
      alpha_w_.push_back(alpha_w);
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
  Body(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv);
  void reset(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv);
  void getBaseCompositeInertia(const Eigen::VectorXd& joint_ori,
                               const Eigen::VectorXd& joint_pos,
                               const Eigen::VectorXd& com_ori,
                               const Eigen::VectorXd& com_pos,
                               const Eigen::VectorXd& mass,
                               const Eigen::VectorXd& inertia) ;
  void getLegComposite(const std::vector<Eigen::Matrix3d>& joint_ori_w,
                              const std::vector<Eigen::Vector3d>& joint_pos_w,
                              const Eigen::VectorXd& com_ori,
                              const Eigen::VectorXd& com_pos,
                              const Eigen::VectorXd& mass,
                              const Eigen::VectorXd& inertia);
  Eigen::Matrix3d getInertia (const Eigen::VectorXd& inertia_vec);
  Eigen::MatrixXd baseMassMatrix (const double& mass, const Eigen::Vector3d& r, const Eigen::Matrix3d& inertia);
  Eigen::VectorXd baseChildMassMatrix (const double& mass, const Eigen::Vector3d& r_com, const Eigen::Matrix3d& I_j,
                                       const Eigen::Vector3d& joint_pos_i, const Eigen::Vector3d& joint_pos_j,
                                       const Eigen::Vector3d& joint_axis_j);
  double childMassMatrix (const double& mass, const Eigen::Vector3d& r_com, const Eigen::Matrix3d& I_j,
                          const Eigen::Vector3d& joint_pos_i, const Eigen::Vector3d& joint_pos_j,
                          const Eigen::Vector3d& joint_axis_i, const Eigen::Vector3d& joint_axis_j);
  void  getComposite (const Eigen::Vector3d& r1_com, const Eigen::Vector3d& r2_com,
                      const double& m1, const double& m2,
                      const Eigen::Matrix3d& I1, const Eigen::Matrix3d& I2,
                      Eigen::Vector3d& r_com, double& mass, Eigen::Matrix3d& I_new);
  void revComposite (const Eigen::Vector3d& r_com, const Eigen::Vector3d& r2_com,
                     const double& m, const double& m2,
                     const Eigen::Matrix3d& I, const Eigen::Matrix3d& I2,
                     Eigen::Vector3d& r1_com, double& m1, Eigen::Matrix3d& I1);

public:
  Eigen::VectorXd gc_, gv_;
  
  // com, mass, inertia of base at world frame (직전 parents 가 base 인 애들 + base 의 composite inertia)
  Eigen::Vector3d base_com;
  double base_mass = 0;
  Eigen::Matrix3d base_inertia;
  
  // com, mass, inertia of legs at world frame (각 다리 들의 com, mass, inertia)
  std::vector<Eigen::Vector3d> leg_COM_w_;
  std::vector<double> leg_Mass_;
  std::vector<Eigen::Matrix3d> leg_Inertia_w_;
  
  // composite com, mass, inertia of legs at world frame (각 다리 들의 composite com, mass, inertia)
  std::vector<Eigen::Vector3d> leg_com_COM_w_;
  std::vector<double> leg_com_Mass_;
  std::vector<Eigen::Matrix3d> leg_com_Inertia_w_;
};

Body::Body(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) { gc_ = gc; gv_ = gv;};

void Body::reset(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
  gc_ = gc;
  gv_ = gv;
  /// 이전에 구해놓은 값들을 초기화, gc, gv가 바뀌었을 경우 다시 구해줘야 하기 때문에..
  base_com.setZero();
  base_mass = 0;
  base_inertia.setZero();
  leg_com_COM_w_.clear();
  leg_com_Mass_.clear();
  leg_com_Inertia_w_.clear();
}

void Body::getBaseCompositeInertia (const Eigen::VectorXd& joint_ori,
                                    const Eigen::VectorXd& joint_pos,
                                    const Eigen::VectorXd& com_ori,
                                    const Eigen::VectorXd& com_pos,
                                    const Eigen::VectorXd& mass,
                                    const Eigen::VectorXd& inertia) {
  /// base 에 붙어있는 fixed joint들을 모투 합쳐 base의 inertia 를 구함
  Eigen::Vector3d gc_pos = gc_.head(3);
  Eigen::Matrix3d gc_rot = quat2rot(gc_.segment(3, 4));
  base_com.setZero();
  base_inertia.setZero();
  
  for (int i = 0; i < mass.size(); i++) {
    Eigen::Vector3d com_pos_w, base_com_new;
    Eigen::Matrix3d com_rot_w, inertia_w;
    // world frame 기준 com의 position, orientation, inertia 구해주기
    com_pos_w = gc_pos + gc_rot * (joint_pos.segment(i * 3, 3) +
                                   rpy2rot(joint_ori.segment(i * 3, 3)) * com_pos.segment(i * 3, 3));
    com_rot_w = gc_rot * rpy2rot(joint_ori.segment(i * 3, 3)) * rpy2rot(com_ori.segment(i * 3, 3));
    inertia_w = com_rot_w * getInertia(inertia.segment(i * 6, 6)) * com_rot_w.transpose();
    
    getComposite(base_com, com_pos_w, base_mass, mass(i), base_inertia, inertia_w,
                 base_com, base_mass,base_inertia);
  };
};

void Body::getLegComposite(const std::vector<Eigen::Matrix3d>& joint_ori_w,
                          const std::vector<Eigen::Vector3d>& joint_pos_w,
                          const Eigen::VectorXd& com_ori,
                          const Eigen::VectorXd& com_pos,
                          const Eigen::VectorXd& mass,
                          const Eigen::VectorXd& inertia) {
  /**
  tree의 composite inertia를 update 해주는 함수 \n
    \n
  { input variables }  \n
  모든 값들은 구하고자 하는 tree 의 전체 값들이 순서대로 들어가 있어야 함. 자료형은 input 자료형 참고 \n
  - joint_ori_w: world frame 기준 joint orientation \n
  - joint_pos_w: world frame 기준 joint position \n
  - com_ori: urdf 상에서 link 의 origin, 즉, parent frame 기준 link 의 origin \n
  - com_pos: urdf 상에서 link 의 position, 즉, parent frame 기준 link 의 position \n
  - mass: urdf 상에서 각 link의 mass \n
  - inertia: urdf 상에서 각 link의 inertia \n
    \n
  { update variables } \n
  - 이 함수를 불러오게 되면 기존에 있던 leg_com_COM_w_, leg_com_Inertia_w_, leg_com_Mass_ std::vector 맨 뒤쪽에 현재 tree의 값들이 추가됨 \n
  **/
  std::vector<Eigen::Vector3d> leg_com_w;
  std::vector<double> leg_mass;
  std::vector<Eigen::Matrix3d> leg_inertia_w;
  
  int num = static_cast<int>(mass.size());
  for (int i = num - 1; i >= 0; i--) {
    Eigen::Vector3d com_pos_w, leg_com_new;
    double leg_mass_new;
    Eigen::Matrix3d com_rot_w, inertia_w, leg_inertia_new;
    com_pos_w = joint_pos_w[i] + joint_ori_w[i] * com_pos.segment(i * 3, 3);
    com_rot_w = joint_ori_w[i] * rpy2rot(com_ori.segment(i * 3, 3));
    inertia_w = com_rot_w * getInertia(inertia.segment(i * 6, 6)) * com_rot_w.transpose();
    
    if (i == num - 1) { // tree의 맨 마지막의 경우, vector에 들어있는게 없으므로 추가만 해줌
      leg_com_w.push_back(com_pos_w);
      leg_mass.push_back(mass(i));
      leg_inertia_w.push_back(inertia_w);
    }
    else {
      getComposite(leg_com_w.back(), com_pos_w, leg_mass.back(), mass(i), leg_inertia_w.back(), inertia_w,
                   leg_com_new, leg_mass_new, leg_inertia_new);
      leg_com_w.push_back(leg_com_new);
      leg_mass.push_back(leg_mass_new);
      leg_inertia_w.push_back(leg_inertia_new);
    }
  }
  
  std::reverse(leg_com_w.begin(), leg_com_w.end());
  std::reverse(leg_mass.begin(), leg_mass.end());
  std::reverse(leg_inertia_w.begin(), leg_inertia_w.end());
  leg_com_COM_w_.insert(leg_com_COM_w_.end(), leg_com_w.begin(), leg_com_w.end());
  leg_com_Mass_.insert(leg_com_Mass_.end(), leg_mass.begin(), leg_mass.end());
  leg_com_Inertia_w_.insert(leg_com_Inertia_w_.end(), leg_inertia_w.begin(), leg_inertia_w.end());
};

Eigen::MatrixXd Body::baseMassMatrix (const double& mass, const Eigen::Vector3d& com, const Eigen::Matrix3d& inertia) {
  /**
   CRBA 를 이용해 floating base의 Mass Matrix를 구해주는 함수 \n
   \n
   { input variables } \n
   - mass: robot의 전체 mass \n
   - com: robot의 전체 com \n
   - inertia: robot의 전체 inertia \n
   **/
  Eigen::Vector3d r = gc_.head(3) - com;
  
  Eigen::Matrix<double, 6, 6> M;
  M.block<3,3>(0,0) = mass*Eigen::Matrix3d::Identity();
  M.block<3,3>(3,0) = -mass*skew(r);
  M.block<3,3>(0,3) = mass*skew(r);
  M.block<3,3>(3,3) = inertia - mass*skew(r)*skew(r);
  
  return M;
}

Eigen::VectorXd Body::baseChildMassMatrix (const double& mass, const Eigen::Vector3d& r_com, const Eigen::Matrix3d& I_j,
                                           const Eigen::Vector3d& joint_pos_i, const Eigen::Vector3d& joint_pos_j,
                                           const Eigen::Vector3d& joint_axis_j) {
  /**
   CRBA 를 이용해 child tree와 floating base 사이의 Mass Matrix를 구해주는 함수 \n
   \n
   { input variables } \n
   **/
  Eigen::Vector3d r_jc = r_com - joint_pos_j;
  Eigen::Vector3d r_ij = joint_pos_j - joint_pos_i;
  Eigen::Matrix<double,6,6> s_i; s_i.setZero();
  s_i.block<3,3>(0,0) = Eigen::Matrix3d::Identity();
  s_i.block<3,3>(3,3) = Eigen::Matrix3d::Identity();
  
  Eigen::VectorXd s_j; s_j.setZero(6);
  s_j.segment(3,3) = joint_axis_j;
  
  Eigen::Matrix<double, 6, 6> Mc;
  Mc.block<3,3>(0,0) = mass*Eigen::Matrix3d::Identity();
  Mc.block<3,3>(3,0) = mass*skew(r_jc);
  Mc.block<3,3>(0,3) = - mass*skew(r_jc);
  Mc.block<3,3>(3,3) = I_j - mass*skew(r_jc)*skew(r_jc);
  
  Eigen::Matrix<double,6,6> A; A.setZero();
  A.block<3,3>(0,0) = Eigen::Matrix3d::Identity();
  A.block<3,3>(3,3) = Eigen::Matrix3d::Identity();
  A.block<3,3>(0,3) = - skew(r_ij);
  
  return s_j.transpose() * Mc * A * s_i;
}

double Body::childMassMatrix (const double& mass, const Eigen::Vector3d& r_com, const Eigen::Matrix3d& I_j,
                              const Eigen::Vector3d& joint_pos_i, const Eigen::Vector3d& joint_pos_j,
                              const Eigen::Vector3d& joint_axis_i, const Eigen::Vector3d& joint_axis_j) {
  /**
   CRBA 를 이용해 child tree의 Mass Matrix를 구해주는 함수 \n
   \n
   { input variables } \n
   **/
  Eigen::Vector3d r_jc = r_com - joint_pos_j;
  Eigen::Vector3d r_ij = joint_pos_j - joint_pos_i;
  Eigen::VectorXd s_i, s_j; s_i.setZero(6); s_j.setZero(6);
  s_i.segment(3,3) = joint_axis_i;
  s_j.segment(3,3) = joint_axis_j;
  
  Eigen::Matrix<double, 6, 6> Mc;
  Mc.block<3,3>(0,0) = mass*Eigen::Matrix3d::Identity();
  Mc.block<3,3>(3,0) = mass*skew(r_jc);
  Mc.block<3,3>(0,3) = - mass*skew(r_jc);
  Mc.block<3,3>(3,3) = I_j - mass*skew(r_jc)*skew(r_jc);
  
  Eigen::Matrix<double,6,6> A; A.setZero();
  A.block<3,3>(0,0) = Eigen::Matrix3d::Identity();
  A.block<3,3>(3,3) = Eigen::Matrix3d::Identity();
  A.block<3,3>(0,3) = - skew(r_ij);
  return s_j.transpose() * Mc * A * s_i;
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
  gravity << 0, 0, 9.81;
  
  // forward kinematics
  Body body(gc, gv);
  body.reset(gc, gv);
  Joint joint(gc, gv);
  joint.gravity << gravity;
  joint.jointKinematics(); // joint kinematics를 구해줌
  
  // Step 1. get center of mass of each body
  /// change urdf!
  // << Base >>
  Eigen::VectorXd base_joint_ori, base_joint_pos, base_com_ori, base_com_pos, base_mass, base_inertia;
  int base_child_num = 11; // base에 붙어있는 fixed joint의 개수
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
  
  body.getBaseCompositeInertia(base_joint_ori, base_joint_pos, base_com_ori, base_com_pos, base_mass, base_inertia);
  
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
  body.getLegComposite(lf_joint_rot_w,lf_joint_pos_w,
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
  
  std::vector<Eigen::Matrix3d> rf_joint_rot_w(joint.joint_rot_w_.begin() + rf_child_num + 2, joint.joint_rot_w_.begin() + rf_child_num*(1+1) + 2);
  std::vector<Eigen::Vector3d> rf_joint_pos_w(joint.joint_pos_w_.begin() + rf_child_num + 2, joint.joint_pos_w_.begin() + rf_child_num*(1+1) + 2);
  body.getLegComposite(rf_joint_rot_w,rf_joint_pos_w,
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
  
  std::vector<Eigen::Matrix3d> lh_joint_rot_w(joint.joint_rot_w_.begin() + lh_child_num*2 + 3, joint.joint_rot_w_.begin() + lh_child_num*(2+1) + 3);
  std::vector<Eigen::Vector3d> lh_joint_pos_w(joint.joint_pos_w_.begin() + lh_child_num*2 + 3, joint.joint_pos_w_.begin() + lh_child_num*(2+1) + 3);
  body.getLegComposite(lh_joint_rot_w,lh_joint_pos_w,
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
  
  std::vector<Eigen::Matrix3d> rh_joint_rot_w(joint.joint_rot_w_.begin() + rh_child_num*3 + 4, joint.joint_rot_w_.begin() + rh_child_num*(3+1) + 4);
  std::vector<Eigen::Vector3d> rh_joint_pos_w(joint.joint_pos_w_.begin() + rh_child_num*3 + 4, joint.joint_pos_w_.begin() + rh_child_num*(3+1) + 4);
  body.getLegComposite(rh_joint_rot_w,rh_joint_pos_w,
                              rh_com_ori, rh_com_pos, rh_mass, rh_inertia);
  
  // CRBA
//  // for debugging : index check!
//  std::cout << body.leg_com_COM_w_.size() << std::endl;
//  for (int i=0; i<body.leg_com_COM_w_.size(); i++){
//    std::cout << "\nidx: " << i << std::endl;
//    std::cout << "pos: " << body.leg_com_COM_w_[i].transpose() << std::endl;
//    std::cout << "mass: " << body.leg_com_Mass_[i] << std::endl;
//  }
  /// change urdf ! check index ! index is revolute joint number
  std::vector<int> body_leg_idx = {0,3,6, 9,12,15, 18,21,24, 27,30,33}; // 최초 revolute joint 이전의 fixed joint 때문에 다름
  std::vector<int> joint_leg_idx = {1,4,7, 11,14,17, 21,24,27, 31,34,37}; // lf(haa, hfe, kfe), rf, lh, rh
  int leg_num = 4; // 다리의 개수
  int leg_joint_num = 3; // 각 다리의 revolute joint 개수
  
  // Base Mass Matrix
  Eigen::Vector3d com = body.base_com;
  double mass = body.base_mass;
  Eigen::Matrix3d inertia = body.base_inertia;
  for (int k=0; k<leg_num; k++) {
    body.getComposite(com, body.leg_com_COM_w_[body_leg_idx[k*leg_joint_num]],
                      mass, body.leg_com_Mass_[body_leg_idx[k*leg_joint_num]], inertia,
                      body.leg_com_Inertia_w_[body_leg_idx[k*leg_joint_num]],
                      com, mass, inertia);
  }
  
  /// change urdf ! 알맞은 mass matrix 크기로 바꿔줘야 함
  Eigen::Matrix<double, 18, 18> Mass_mat; Mass_mat.setZero();
  Mass_mat.block<6,6>(0,0) = body.baseMassMatrix(mass, com, inertia);
  
  // Base Child Mass Matrix
  for (int k=0; k<leg_num; k++) {
    for (int j=0; j<leg_joint_num; j++) {
      Eigen::VectorXd joint_axis_j =
        joint.joint_sign_[k*leg_joint_num+j] * joint.joint_rot_w_[joint_leg_idx[k*leg_joint_num+j]].block<3,1>(0,0);
      Mass_mat.block<1,6>((k*leg_joint_num+j)+6, 0) = body.baseChildMassMatrix
        (body.leg_com_Mass_[body_leg_idx[k*leg_joint_num+j]], body.leg_com_COM_w_[body_leg_idx[k*leg_joint_num+j]],
         body.leg_com_Inertia_w_[body_leg_idx[k*leg_joint_num+j]],
         gc.head(3), joint.joint_pos_w_[joint_leg_idx[k*leg_joint_num+j]], joint_axis_j);
    }
  }
  
  // Child Mass Matrix
  for (int k=0; k<leg_num; k++) {
    for (int i=0; i<leg_joint_num; i++) {
      for (int j=i; j<leg_joint_num; j++) {
        Eigen::Vector3d joint_axis_i = joint.joint_sign_[k*leg_joint_num+i] * joint.joint_rot_w_[joint_leg_idx[k*leg_joint_num+i]].block<3,1>(0,0);
        Eigen::VectorXd joint_axis_j = joint.joint_sign_[k*leg_joint_num+j] * joint.joint_rot_w_[joint_leg_idx[k*leg_joint_num+j]].block<3,1>(0,0);
        Mass_mat((k*leg_joint_num+j)+6, (k*leg_joint_num+i)+6) = body.childMassMatrix(
          body.leg_com_Mass_[body_leg_idx[k*leg_joint_num+j]], body.leg_com_COM_w_[body_leg_idx[k*leg_joint_num+j]],
          body.leg_com_Inertia_w_[body_leg_idx[k*leg_joint_num+j]],
          joint.joint_pos_w_[joint_leg_idx[k*leg_joint_num+i]], joint.joint_pos_w_[joint_leg_idx[k*leg_joint_num+j]],
          joint_axis_i, joint_axis_j);
      }
    }
  }
  
  Mass_mat = Mass_mat.selfadjointView<Eigen::Lower>();
  
  
  // RNE
  std::vector<Eigen::Vector3d> F_w_;
  std::vector<Eigen::Vector3d> tau_w_;
  Eigen::VectorXd base_gen_F; base_gen_F.setZero(6);
  
  for (int k=0; k<leg_num; k++) {
    std::vector<Eigen::Vector3d> F_w; F_w.clear();
    std::vector<Eigen::Vector3d> tau_w; tau_w.clear();
    for (int i=leg_joint_num-1; i>=0; i--) {
      int idx = k*leg_joint_num+i;
      
      Eigen::Vector3d body_com;
      double body_mass;
      Eigen::Matrix3d body_inertia;
      
      // get body com, mass, inertia
      if (i<leg_joint_num-1) {
        body.revComposite(body.leg_com_COM_w_[body_leg_idx[idx]], body.leg_com_COM_w_[body_leg_idx[idx+1]],
                          body.leg_com_Mass_[body_leg_idx[idx]], body.leg_com_Mass_[body_leg_idx[idx+1]],
                          body.leg_com_Inertia_w_[body_leg_idx[idx]], body.leg_com_Inertia_w_[body_leg_idx[idx+1]],
                          body_com, body_mass, body_inertia);
      }
      else {
        body_com = body.leg_com_COM_w_[body_leg_idx[idx]];
        body_mass = body.leg_com_Mass_[body_leg_idx[idx]];
        body_inertia = body.leg_com_Inertia_w_[body_leg_idx[idx]];
      }
      
      // get Spatial Inertia Matrix
      Eigen::VectorXd r_jc = body_com - joint.joint_pos_w_[joint_leg_idx[idx]];
      Eigen::Matrix<double,6,6> Mc; Mc.setZero();
      Mc.block<3,3>(0,0) = body_mass * Eigen::Matrix3d::Identity();
      Mc.block<3,3>(3,0) = body_mass*skew(r_jc);
      Mc.block<3,3>(0,3) = -body_mass*skew(r_jc);
      Mc.block<3,3>(3,3) = body_inertia - body_mass*skew(r_jc)*skew(r_jc);

      // get Generalized Force
      Eigen::Matrix<double,6,1> gen_F;
      Eigen::Vector3d angvel = joint.angvel_w_[1+idx];
      // + b
      gen_F.segment(0,3) = body_mass*skew(angvel)*skew(angvel)*r_jc;
      gen_F.segment(3,3) = skew(angvel) * (body_inertia - body_mass*skew(r_jc)*skew(r_jc)) * angvel;
      // + force applied from the child
      if (i<leg_joint_num-1) { // not leaf body
        gen_F.segment(0,3) += F_w.back();
        gen_F.segment(3,3) += tau_w.back() + skew(joint.joint_pos_w_[joint_leg_idx[idx+1]] - joint.joint_pos_w_[joint_leg_idx[idx]])*F_w.back();
      }
      // + M*[a alpha]
      gen_F += Mc * joint.gen_a_[1+idx];
      F_w.push_back(gen_F.segment(0,3));
      tau_w.push_back(gen_F.segment(3,3));
    }
    // base_gen_F (+ force applied from the child)
    base_gen_F.segment(0,3) += F_w.back();
    base_gen_F.segment(3,3) += tau_w.back() + skew(joint.joint_pos_w_[joint_leg_idx[k*leg_joint_num]] - gc.head(3))*F_w.back();
    // Since we calculated it from child, we turn it over and save it.
    std::reverse(F_w.begin(), F_w.end());
    F_w_.insert(F_w_.end(), F_w.begin(), F_w.end());
    std::reverse(tau_w.begin(), tau_w.end());
    tau_w_.insert(tau_w_.end(), tau_w.begin(), tau_w.end());
  }
  
  // base
  // get Spatial Inertia Matrix
  Eigen::VectorXd r_jc = body.base_com - gc.head(3);
  Eigen::Matrix<double,6,6> Mc; Mc.setZero();
  Mc.block<3,3>(0,0) = body.base_mass * Eigen::Matrix3d::Identity();
  Mc.block<3,3>(3,0) = body.base_mass * skew(r_jc);
  Mc.block<3,3>(0,3) = -body.base_mass * skew(r_jc);
  Mc.block<3,3>(3,3) = body.base_inertia - body.base_mass*skew(r_jc)*skew(r_jc);
  
  // get Generalized Force (base_gen_F) (The force applied by child was already added earlier.)
  Eigen::Vector3d angvel = joint.angvel_w_[0];
  // + b
  base_gen_F.segment(0,3) += body.base_mass*skew(angvel)*skew(angvel)*r_jc;
  base_gen_F.segment(3,3) += skew(angvel) * (body.base_inertia - body.base_mass*skew(r_jc)*skew(r_jc)) * angvel;
  // + M*[a alpha]
  base_gen_F += Mc * joint.gen_a_[0];

  // get Nonlinearities
  Eigen::VectorXd b; b.setZero(18);
  b.head(6) << base_gen_F;
  for (int k=0; k<leg_num; k++) {
    for (int i=0; i<leg_joint_num; i++) {
      int idx = k*leg_joint_num + i;
      
      Eigen::VectorXd gen_F; gen_F.setZero(6);
      
      gen_F << F_w_[idx], tau_w_[idx];
      b.segment(6+idx,1) = joint.s_[1+idx].transpose() * gen_F;
    }
  }
  
  return b;
}