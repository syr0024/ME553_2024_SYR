#pragma once


class Body {
public:
    double mass;

    /// In body frame
    Eigen::Matrix3d inertiaB;
    Eigen::Vector3d comPosB;     // <origin xyz="-0.063 7e-05 0.00046"/> in <inertial>
    Eigen::Matrix3d jointRotB;
    Eigen::Vector3d jointPosB;
    Eigen::Vector3d jointAxisB;
    double rotatingAngle;

    /// In world frame
    Eigen::Matrix3d inertiaW;
    Eigen::Vector3d comPosW;
    Eigen::Matrix3d jointRotW;
    Eigen::Vector3d jointPosW;
    Eigen::Vector3d jointAxisW;

    ///
    double jointVel;
    double jointVelDot;
    Eigen::Vector3d linVelW;
    Eigen::Vector3d angVelW;
    Eigen::Vector3d linAccW;
    Eigen::Vector3d angAccW;
    Eigen::Vector3d forceW;
    Eigen::Vector3d torqueW;

    Eigen::VectorXd genVelW;
    Eigen::VectorXd genAccW;
    Eigen::VectorXd genForceW;

    Eigen::VectorXd s, s_dot;

    ///
    double jointTorque;
    Eigen::MatrixXd XBP;
    Eigen::MatrixXd XBP_dot;
    Eigen::MatrixXd MA;
    Eigen::VectorXd BA;

    Body(double mass_, Eigen::Matrix3d inertiaB_, Eigen::Vector3d comPosB_,
         Eigen::Matrix3d jointRotB_, Eigen::Vector3d jointPosB_, Eigen::Vector3d jointAxisB_, double rotatingAngle_, double jointVel_, double jointTorque_):
            mass(mass_), inertiaB(inertiaB_), comPosB(comPosB_),
            jointRotB(jointRotB_), jointPosB(jointPosB_), jointAxisB(jointAxisB_), rotatingAngle(rotatingAngle_), jointVel(jointVel_), jointTorque(jointTorque_)
    {}
};

Eigen::Matrix3d getInertiaMatrix(double ixx, double ixy, double ixz, double iyy, double iyz, double izz){
    Eigen::Matrix3d inertia;
    inertia << ixx, ixy, ixz,
            ixy, iyy, iyz,
            ixz, iyz, izz;
    return inertia;
}
Eigen::Matrix3d skew(Eigen::Vector3d v){
    Eigen::Matrix3d v_skew;
    v_skew << 0, -v(2), v(1),
            v(2), 0, -v(0),
            -v(1), v(0), 0;
    return v_skew;
}
Eigen::Matrix3d rpyToRotationMatrix(double roll, double pitch, double yaw) {
    Eigen::Matrix3d R_roll, R_pitch, R_yaw, R;
    // Calculate individual rotation matrices
    R_roll << 1, 0, 0,
            0, cos(roll), -sin(roll),
            0, sin(roll), cos(roll);
    R_pitch << cos(pitch), 0, sin(pitch),
            0, 1, 0,
            -sin(pitch), 0, cos(pitch);
    R_yaw << cos(yaw), -sin(yaw), 0,
            sin(yaw), cos(yaw), 0,
            0, 0, 1;
    // Calculate total rotation matrix
    R = R_yaw * R_pitch * R_roll;
    return R;
}
Eigen::Matrix3d getRotatingMatrix(Eigen::Vector3d axis, double angle){
    Eigen::Matrix3d rot;
    rot = rpyToRotationMatrix(0, 0, 0);
    if (axis(0) == 1)
        rot = rpyToRotationMatrix(angle, 0, 0);
    if (axis(0) == -1)
        rot = rpyToRotationMatrix(-angle, 0, 0);
    if (axis(1) == 1)
        rot = rpyToRotationMatrix(0, angle, 0);
    if (axis(2) == 1)
        rot = rpyToRotationMatrix(0, 0, angle);
    return rot;
}
Eigen::Matrix3d quatToRotMat(Eigen::VectorXd quat) {
    Eigen::Matrix3d rotMat;
    rotMat <<-2*quat(2)*quat(2)-2*quat(3)*quat(3)+1, 2*quat(1)*quat(2)-2*quat(0)*quat(3),  2*quat(0)*quat(2)+2*quat(1)*quat(3),
            2*quat(0)*quat(3)+2*quat(1)*quat(2), -2*quat(1)*quat(1)-2*quat(3)*quat(3)+1, 2*quat(2)*quat(3)-2*quat(0)*quat(1),
            2*quat(1)*quat(3)-2*quat(0)*quat(2), 2*quat(0)*quat(1)+2*quat(2)*quat(3), -2*quat(1)*quat(1)-2*quat(2)*quat(2)+1;
    return rotMat;
}

void getForwardKinematicsBase(std::vector<Body>& bodies){
    Eigen::Matrix3d parentJointRotW = bodies[0].jointRotW;
    Eigen::Vector3d parentJointPosW = bodies[0].jointPosW;
    bodies[0].comPosW = bodies[0].jointPosW + bodies[0].jointRotW * bodies[0].comPosB;
    bodies[0].inertiaW = bodies[0].jointRotW * bodies[0].inertiaB * bodies[0].jointRotW.transpose();

    for(int i=1; i<bodies.size(); i++) {
        if(i==6){
            continue;
        }
        bodies[i].jointPosW = parentJointPosW + parentJointRotW * bodies[i].jointPosB;
        bodies[i].jointRotW = parentJointRotW * bodies[i].jointRotB
                              * getRotatingMatrix(bodies[i].jointAxisB, bodies[i].rotatingAngle);
        bodies[i].comPosW = bodies[i].jointPosW + bodies[i].jointRotW * bodies[i].comPosB;
        bodies[i].inertiaW = bodies[i].jointRotW * bodies[i].inertiaB * bodies[i].jointRotW.transpose();
    }
    parentJointPosW = bodies[5].jointPosW; parentJointRotW = bodies[5].jointRotW;
    bodies[6].jointPosW = parentJointPosW + parentJointRotW * bodies[6].jointPosB;
    bodies[6].jointRotW = parentJointRotW * bodies[6].jointRotB
                          * getRotatingMatrix(bodies[6].jointAxisB, bodies[6].rotatingAngle);
    bodies[6].comPosW = bodies[6].jointPosW + bodies[6].jointRotW * bodies[6].comPosB;
    bodies[6].inertiaW = bodies[6].jointRotW * bodies[6].inertiaB * bodies[6].jointRotW.transpose();
}
void getForwardKinematics(std::vector<Body>& bodies, const Eigen::VectorXd& gc){
    Eigen::Matrix3d parentJointRotW = quatToRotMat(gc.segment(3,4));
    Eigen::Vector3d parentJointPosW = gc.head(3);
    for(int i=0; i<bodies.size(); i++) {
        bodies[i].jointPosW = parentJointPosW + parentJointRotW * bodies[i].jointPosB;
        bodies[i].jointRotW = parentJointRotW * bodies[i].jointRotB
                              * getRotatingMatrix(bodies[i].jointAxisB, bodies[i].rotatingAngle);
        parentJointRotW = bodies[i].jointRotW; parentJointPosW = bodies[i].jointPosW;
        bodies[i].comPosW = bodies[i].jointPosW + bodies[i].jointRotW * bodies[i].comPosB;
        bodies[i].inertiaW = bodies[i].jointRotW * bodies[i].inertiaB * bodies[i].jointRotW.transpose();
        bodies[i].jointAxisW = bodies[i].jointRotW * bodies[i].jointAxisB;
    }
}

void getForwardVelocities(std::vector<Body>& bodies, Body& base_Composite){
    Eigen::VectorXd s(6), s_dot(6), genVelW(6);
    Eigen::MatrixXd XBP(6, 6), XBP_dot(6, 6);
    s.setZero(); s_dot.setZero(); genVelW.setZero();
    XBP.setIdentity(); XBP_dot.setZero();

    for(int i=0; i<bodies.size(); i++) {
        if(i==0){
            bodies[i].linVelW = base_Composite.linVelW + skew(base_Composite.angVelW)*(bodies[i].jointPosW - base_Composite.jointPosW);
            bodies[i].angVelW = base_Composite.angVelW + bodies[i].jointAxisW * bodies[i].jointVel;
            genVelW.head(3) = bodies[i].linVelW; genVelW.tail(3) = bodies[i].angVelW;
            bodies[i].genVelW = genVelW;
            s.tail(3) = bodies[i].jointAxisW;
            s_dot.tail(3) = skew(bodies[i].angVelW) * bodies[i].jointAxisW;
            bodies[i].s = s; bodies[i].s_dot = s_dot;
            XBP.bottomLeftCorner(3, 3) = skew(bodies[i].jointPosW - base_Composite.jointPosW);
//            XBP_dot.bottomLeftCorner(3, 3) = skew(base_Composite.angVelW) * XBP.bottomLeftCorner(3, 3);
            XBP_dot.bottomLeftCorner(3, 3) = skew(bodies[i].linVelW - base_Composite.linVelW);
            bodies[i].XBP = XBP;
            bodies[i].XBP_dot = XBP_dot;
        }
        else{
            bodies[i].linVelW = bodies[i-1].linVelW + skew(bodies[i-1].angVelW)*(bodies[i].jointPosW - bodies[i-1].jointPosW);
            bodies[i].angVelW = bodies[i-1].angVelW + bodies[i].jointAxisW * bodies[i].jointVel;
            genVelW.head(3) = bodies[i].linVelW; genVelW.tail(3) = bodies[i].angVelW;
            bodies[i].genVelW = genVelW;
            s.tail(3) = bodies[i].jointAxisW;
            s_dot.tail(3) = skew(bodies[i].angVelW) * bodies[i].jointAxisW;
            bodies[i].s = s; bodies[i].s_dot = s_dot;
            XBP.bottomLeftCorner(3, 3) = skew(bodies[i].jointPosW - bodies[i-1].jointPosW);
//            XBP_dot.bottomLeftCorner(3, 3) = skew(bodies[i-1].angVelW) * XBP.bottomLeftCorner(3, 3);
            XBP_dot.bottomLeftCorner(3, 3) = skew(bodies[i].linVelW - bodies[i-1].linVelW);
            bodies[i].XBP = XBP;
            bodies[i].XBP_dot = XBP_dot;
        }
    }
}
Eigen::MatrixXd getSpatialInertialMatrix(Body& body){
    Eigen::Matrix3d skewRac;
    skewRac = skew(body.comPosW - body.jointPosW);
    Eigen::MatrixXd IS;
    IS.setZero(6, 6);

    IS.topLeftCorner(3,3) = body.mass * Eigen::Matrix3d::Identity();
    IS.topRightCorner(3,3) = -body.mass * skewRac;
    IS.bottomLeftCorner(3,3) = -IS.topRightCorner(3,3);
    IS.bottomRightCorner(3,3) = body.inertiaW - body.mass * skewRac * skewRac;
    return IS;
}
Eigen::MatrixXd getFictitiousForces(Body& body){
    Eigen::Matrix3d skewRac;
    Eigen::Matrix3d skewOmega;
    skewRac = skew(body.comPosW - body.jointPosW);
    skewOmega = skew(body.angVelW);
    Eigen::VectorXd b(6);
    b.setZero();

    b.head(3) = body.mass*skewOmega*skewOmega*(body.comPosW - body.jointPosW);
    b.tail(3) = skewOmega*(body.inertiaW - body.mass*skewRac*skewRac)*body.angVelW;
    return b;
}
Body getCompositeBody(Body& bodyParent, Body& bodyChild){
    double mass = bodyParent.mass + bodyChild.mass;
    Eigen::Matrix3d inertiaW, inertiaB;
    Eigen::Vector3d comPosW, comPosB;
    inertiaB.setZero(); comPosB.setZero();
    comPosW = (bodyParent.mass * bodyParent.comPosW + bodyChild.mass * bodyChild.comPosW)/mass;

    Eigen::Vector3d deltaComParentW = bodyParent.comPosW - comPosW;
    Eigen::Vector3d deltaComChildW = bodyChild.comPosW - comPosW;
    inertiaW = bodyParent.inertiaW + bodyChild.inertiaW
               - bodyParent.mass * skew(deltaComParentW) * skew(deltaComParentW)
               - bodyChild.mass * skew(deltaComChildW) * skew(deltaComChildW);

    Body compositeBody = bodyParent;
    compositeBody.mass = mass;
    compositeBody.inertiaW = inertiaW;
    compositeBody.comPosW = comPosW;
    return compositeBody;
}
std::vector<Body> getConmpositeLeg(std::vector<Body> leg){
    std::vector<Body> leg_bodies;
    Body HIP_Composite_ = getCompositeBody(leg[2], leg[3]);
    Body HIP_Composite = getCompositeBody(leg[1], HIP_Composite_);
    leg_bodies.push_back(HIP_Composite);

    Body THIGH_Composite_ = getCompositeBody(leg[5], leg[6]);
    Body THIGH_Composite = getCompositeBody(leg[4], THIGH_Composite_);
    leg_bodies.push_back(THIGH_Composite);

    Body SHANK_Composite_ = getCompositeBody(leg[8], leg[9]);
    Body SHANK_Composite = getCompositeBody(leg[7], SHANK_Composite_);
    leg_bodies.push_back(SHANK_Composite);

    return leg_bodies;
}
void getArticulatedMandB(std::vector<Body>& leg){
    Eigen::VectorXd BA(6), temp2(6);
    Eigen::MatrixXd MA(6,6), XBP(6,6), XBP_dot(6,6);
    double temp1;
    BA.setZero(); temp2.setZero();
    MA.setIdentity();
    for(int i = leg.size()-1; i>=0; i--) {
        MA = getSpatialInertialMatrix(leg[i]);
        BA = getFictitiousForces(leg[i]);
        if (i < leg.size() - 1) {
            XBP = leg[i+1].XBP; XBP_dot = leg[i+1].XBP_dot;
            MA += XBP * leg[i+1].MA *
               (-leg[i+1].s * (leg[i+1].s.transpose() * leg[i+1].MA * leg[i+1].s).inverse() * leg[i+1].s.transpose()*leg[i+1].MA*XBP.transpose()
               + XBP.transpose());
            temp2 = leg[i+1].s_dot*leg[i+1].jointVel + XBP_dot.transpose()*leg[i].genVelW;
            temp1 = leg[i+1].jointTorque - leg[i+1].s.transpose()*leg[i+1].MA*temp2 - leg[i+1].s.transpose()*leg[i+1].BA;

            BA += XBP * (leg[i+1].MA
                    * (leg[i+1].s * (leg[i+1].s.transpose() * leg[i+1].MA * leg[i+1].s).inverse() * temp1
                    + leg[i+1].s_dot*leg[i+1].jointVel + XBP_dot.transpose()*leg[i].genVelW)
                    + leg[i+1].BA);
//            std::cout << "b: \n" << getFictitiousForces(leg[i]).transpose() << std::endl;
//          std::cout << "\nXBP_dot: \n" << leg[i+1].XBP_dot << std::endl;
//          std::cout << "XBP: \n" << leg[i+1].XBP << std::endl;
//          std::cout << "s: \n" << leg[i+1].s.transpose() << std::endl;
//          std::cout << "s_dot: \n" << leg[i+1].s_dot.transpose() << std::endl;
//          std::cout << "w_AP: \n" << leg[i].genVelW.transpose() << std::endl;
//          std::cout << "temp1: " << temp1 << std::endl;
        }
        leg[i].MA = MA;
        leg[i].BA = BA;
//        std::cout << "body Mc: \n" << getSpatialInertialMatrix(leg[i]) << std::endl;
//        std::cout << "MA: \n" << MA << std::endl;
//        std::cout << "b: \n" << getFictitiousForces(leg[i]).transpose() << std::endl;
        std::cout << "BA: \n" << BA.transpose() << std::endl;
    }
}

void getArticulatedMandBforBase(Body& base, std::vector<Body>& thighs){
    Eigen::VectorXd BA(6), temp2(6);
    Eigen::MatrixXd MA(6,6), XBP(6,6), XBP_dot(6,6);
    double temp1;
    BA.setZero(); temp2.setZero();
    MA.setIdentity();
    MA = getSpatialInertialMatrix(base);
    BA = getFictitiousForces(base);
    for(int i=0; i<thighs.size(); i++) {
        XBP = thighs[i].XBP; XBP_dot = thighs[i].XBP_dot;
        MA += XBP * thighs[i].MA *
              (-thighs[i].s * (thighs[i].s.transpose() * thighs[i].MA * thighs[i].s).inverse() * thighs[i].s.transpose()*thighs[i].MA*XBP.transpose()
               + XBP.transpose());
        temp2 = thighs[i].s_dot*thighs[i].jointVel + XBP_dot.transpose()*base.genVelW;
        temp1 = thighs[i].jointTorque - thighs[i].s.transpose()*thighs[i].MA*temp2 - thighs[i].s.transpose()*thighs[i].BA;

        BA += XBP * (thighs[i].MA
                     * (thighs[i].s * (thighs[i].s.transpose() * thighs[i].MA * thighs[i].s).inverse() * temp1
                         + thighs[i].s_dot*thighs[i].jointVel + XBP_dot.transpose()*base.genVelW)
                     + thighs[i].BA);
    }
    base.MA = MA;
    base.BA = BA;
}
void getForwardAccelerationsforBase(Body& base, Eigen::VectorXd& genAccW){
    Eigen::VectorXd uB_dot(6), temp(6);
    temp.setZero();
    temp = base.genVelW;
    uB_dot.setZero();
    uB_dot = base.MA.inverse()*(base.genForceW - base.BA);
    base.genAccW = uB_dot;
    genAccW.head(6) = uB_dot;
}
void getForwardAccelerations(std::vector<Body>& leg, Body& base, Eigen::VectorXd& genAccW, int legNumber){
    Eigen::VectorXd BA(6);
    Eigen::MatrixXd MA(6,6), XBP(6,6), XBP_dot(6,6);
    BA.setZero(); MA.setIdentity();

    for(int i=0; i<leg.size(); i++){
        if(i==0){
            XBP = leg[i].XBP; XBP_dot = leg[i].XBP_dot;
            leg[i].jointVelDot = 1/(leg[i].s.transpose() * leg[i].MA * leg[i].s)
                                 * (leg[i].jointTorque - leg[i].s.transpose()*leg[i].MA
                                                         *(leg[i].s_dot*leg[i].jointVel + XBP_dot.transpose()*base.genVelW + XBP.transpose()*base.genAccW)
                                    - leg[i].s.transpose() * leg[i].BA);
            leg[i].genAccW = leg[i].s*leg[i].jointVelDot + leg[i].s_dot*leg[i].jointVel
                             + XBP_dot.transpose()*base.genVelW + XBP.transpose()*base.genAccW;
        }
        else{
            XBP = leg[i].XBP; XBP_dot = leg[i].XBP_dot;
            leg[i].jointVelDot = 1/(leg[i].s.transpose() * leg[i].MA * leg[i].s)
                                 * (leg[i].jointTorque - leg[i].s.transpose()*leg[i].MA
                                 *(leg[i].s_dot*leg[i].jointVel + XBP_dot.transpose()*leg[i-1].genVelW + XBP.transpose()*leg[i-1].genAccW)
                                 - leg[i].s.transpose() * leg[i].BA);
            leg[i].genAccW = leg[i].s*leg[i].jointVelDot + leg[i].s_dot*leg[i].jointVel
                             + XBP_dot.transpose()*leg[i-1].genVelW + XBP.transpose()*leg[i-1].genAccW;
        }
        genAccW[6 + legNumber + i] = leg[i].jointVelDot;
    }
}

/// do not change the name of the method
inline Eigen::VectorXd computeGeneralizedAcceleration(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, const Eigen::VectorXd& gf) {
    /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!
    Eigen::VectorXd gc_, gv_, gf_;
    gc_ = gc; gv_ = gv; gf_ = gf;
    double Mag;
    Mag = gc_.segment(3,4).norm();
    gc_.segment(3,4) = gc_.segment(3,4) / Mag;
    Eigen::Matrix3d tempMatrix; tempMatrix.setZero();
    std::vector<Body> base_bodies;

    Body base(6.222, getInertiaMatrix(0.017938806, 0.00387963, 0.001500772, 0.370887745, 6.8963e-05, 0.372497653),
              Eigen::Vector3d{-0.018, -0.002, 0.024}, rpyToRotationMatrix(0, 0, 0),
              Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    base_bodies.push_back(base);
    Body face_front(0.73, getInertiaMatrix(0.005238611, 1.7609e-05, 7.2167e-05, 0.002643098, 1.9548e-05, 0.004325938),
                    Eigen::Vector3d{0.042, -0.001, 0.004}, rpyToRotationMatrix(0, 0, 0),
                    Eigen::Vector3d{0.4145, 0, 0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    base_bodies.push_back(face_front);
    Body face_rear(0.73, getInertiaMatrix(0.005238611, 1.7609e-05, 7.2167e-05, 0.002643098, 1.9548e-05, 0.004325938),
                   Eigen::Vector3d{0.042, -0.001, 0.004}, rpyToRotationMatrix(0, 0, 3.14159265359),
                   Eigen::Vector3d{-0.4145, 0, 0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    base_bodies.push_back(face_rear);
    Body battery(5.53425, getInertiaMatrix(0.00749474794, 0.00016686282, 7.82763e-05, 0.0722338913, 1.42902e-06, 0.07482717535),
                 Eigen::Vector3d{-0.00067, -0.00023, -0.03362}, rpyToRotationMatrix(0, 0, 0),
                 Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    base_bodies.push_back(battery);
    Body docking_hatch_cover(0.065, getInertiaMatrix(0.00063283, 0.0, 3.45e-07, 0.00110971, 0.0, 0.00171883),
                             Eigen::Vector3d{-0.003, 0.0, 0.005}, rpyToRotationMatrix(0, 0, 0),
                             Eigen::Vector3d{0.343, 0.0, -0.07}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    base_bodies.push_back(docking_hatch_cover);
    Body lidar_cage(0.0, getInertiaMatrix(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                    Eigen::Vector3d{0, 0, 0}, rpyToRotationMatrix(0, 0, 0),
                    Eigen::Vector3d{-0.364, 0.0, 0.0735}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    base_bodies.push_back(lidar_cage);
    Body lidar(0.695, getInertiaMatrix(0.000846765, 6.9565e-05, 0.00027111, 0.001367583, 5.8984e-05, 0.001363673),
               Eigen::Vector3d{-0.012, 0.001, -0.008}, rpyToRotationMatrix(0, 0, -1.57079632679),
               Eigen::Vector3d{0.0, 0.0, 0.0687}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    base_bodies.push_back(lidar);
    Body hatch(0.142, getInertiaMatrix(0.001, 0.001, 0.001, 0.001, 0.001, 0.001),
               Eigen::Vector3d{0.116, 0.0, 0.0758}, rpyToRotationMatrix(0, 0, 0),
               Eigen::Vector3d{0.0, 0.0, 0.0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    base_bodies.push_back(hatch);
    Body LF_HAA(2.04, getInertiaMatrix(0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827),
                Eigen::Vector3d{-0.063, 7e-05, 0.00046}, rpyToRotationMatrix(2.61799387799, 0, 0),
                Eigen::Vector3d{0.2999, 0.104, 0.0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    Body RF_HAA(2.04, getInertiaMatrix(0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827),
                Eigen::Vector3d{-0.063, 7e-05, 0.00046}, rpyToRotationMatrix(-2.61799387799, 0, 0),
                Eigen::Vector3d{0.2999, -0.104, 0.0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    Body LH_HAA(2.04, getInertiaMatrix(0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827),
                Eigen::Vector3d{-0.063, 7e-05, 0.00046}, rpyToRotationMatrix(-2.61799387799, 0, -3.14159265359),
                Eigen::Vector3d{-0.2999, 0.104, 0.0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    Body RH_HAA(2.04, getInertiaMatrix(0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827),
                Eigen::Vector3d{-0.063, 7e-05, 0.00046}, rpyToRotationMatrix(2.61799387799, 0, -3.14159265359),
                Eigen::Vector3d{-0.2999, -0.104, 0.0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);

    /// LF leg
    std::vector<Body> LF_leg;
    LF_leg.push_back(LF_HAA);
    Body LF_HIP(0.001, getInertiaMatrix(0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001),
                Eigen::Vector3d{0, 0, 0}, rpyToRotationMatrix(0, 0, 0),
                Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{1, 0, 0}, gc_[7], gv_[6], gf_[6]);
    LF_leg.push_back(LF_HIP);
    Body LF_hip_fixed(0.74, getInertiaMatrix(0.001393106, 8.4012e-05, 2.3378e-05, 0.003798579, 7.1319e-05, 0.003897509),
                      Eigen::Vector3d{0.048, 0.008, -0.003}, rpyToRotationMatrix(-2.61799387799, 0, 0.0),
                      Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    LF_leg.push_back(LF_hip_fixed);
    Body LF_HFE(2.04, getInertiaMatrix(0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827),
                Eigen::Vector3d{-0.063, 7e-05, 0.00046}, rpyToRotationMatrix(0, 0, 1.57079632679),
                Eigen::Vector3d{0.0599, 0.08381, 0.0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    LF_leg.push_back(LF_HFE);

    Body LF_THIGH(0.001, getInertiaMatrix(0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001),
                  Eigen::Vector3d{0, 0, 0}, rpyToRotationMatrix(0, 0, 0),
                  Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{1, 0, 0}, gc_[8], gv_[7], gf_[7]);
    LF_leg.push_back(LF_THIGH);
    Body LF_thigh_fixed(1.03, getInertiaMatrix(0.018644469, 5.2e-08, 1.0157e-05, 0.019312599, 0.002520077, 0.002838361),
                        Eigen::Vector3d{0.0, 0.018, -0.169}, rpyToRotationMatrix(0, 0, -1.57079632679),
                        Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    LF_leg.push_back(LF_thigh_fixed);
    Body LF_KFE(2.04, getInertiaMatrix(0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827),
                Eigen::Vector3d{-0.063, 7e-05, 0.00046}, rpyToRotationMatrix(0, 0, 1.57079632679),
                Eigen::Vector3d{0.0, 0.1003, -0.285}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    LF_leg.push_back(LF_KFE);

    Body LF_SHANK(0.001, getInertiaMatrix(0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001),
                  Eigen::Vector3d{0, 0, 0}, rpyToRotationMatrix(0, 0, 0),
                  Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{1, 0, 0}, gc_[9], gv_[8], gf_[8]);
    LF_leg.push_back(LF_SHANK);
    Body LF_shank_fixed(0.33742, getInertiaMatrix(0.00032748005, 2.142561e-05, 1.33942e-05, 0.00110974122, 7.601e-08, 0.00089388521),
                        Eigen::Vector3d{0.03463, 0.00688, 0.00098}, rpyToRotationMatrix(0, 0, -1.57079632679),
                        Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    LF_leg.push_back(LF_shank_fixed);
    Body LF_FOOT(0.25, getInertiaMatrix(0.00317174097, 2.63048e-06, 6.815581e-05, 0.00317174092, 6.815583e-05, 8.319196e-05),
                 Eigen::Vector3d{0.00948, -0.00948, 0.1468}, rpyToRotationMatrix(0, 0, 0),
                 Eigen::Vector3d{0.08795, 0.01305, -0.33797}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    LF_leg.push_back(LF_FOOT);

    /// RF leg
    std::vector<Body> RF_leg;
    RF_leg.push_back(RF_HAA);
    Body RF_HIP(0.001, getInertiaMatrix(0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001),
                Eigen::Vector3d{0, 0, 0}, rpyToRotationMatrix(0, 0, 0),
                Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{1, 0, 0}, gc_[10], gv_[9], gf_[9]);
    RF_leg.push_back(RF_HIP);
    Body RF_hip_fixed(0.74, getInertiaMatrix(0.001393106, -8.4012e-05, 2.3378e-05, 0.003798579, -7.1319e-05, 0.003897509),
                      Eigen::Vector3d{0.048, -0.008, -0.003}, rpyToRotationMatrix(2.61799387799, 0, 0.0),
                      Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    RF_leg.push_back(RF_hip_fixed);
    Body RF_HFE(2.04, getInertiaMatrix(0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827),
                Eigen::Vector3d{-0.063, 7e-05, 0.00046}, rpyToRotationMatrix(0, 0, -1.57079632679),
                Eigen::Vector3d{0.0599, -0.08381, 0.0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    RF_leg.push_back(RF_HFE);

    Body RF_THIGH(0.001, getInertiaMatrix(0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001),
                  Eigen::Vector3d{0, 0, 0}, rpyToRotationMatrix(0, 0, 0),
                  Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{-1, 0, 0}, gc_[11], gv_[10], gf_[10]);
    RF_leg.push_back(RF_THIGH);
    Body RF_thigh_fixed(1.03, getInertiaMatrix(0.018644469, -5.2e-08, 1.0157e-05, 0.019312599, -0.002520077, 0.002838361),
                        Eigen::Vector3d{0.0, -0.018, -0.169}, rpyToRotationMatrix(0, 0, 1.57079632679),
                        Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    RF_leg.push_back(RF_thigh_fixed);
    Body RF_KFE(2.04, getInertiaMatrix(0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827),
                Eigen::Vector3d{-0.063, 7e-05, 0.00046}, rpyToRotationMatrix(0, 0, -1.57079632679),
                Eigen::Vector3d{0.0, -0.1003, -0.285}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    RF_leg.push_back(RF_KFE);

    Body RF_SHANK(0.001, getInertiaMatrix(0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001),
                  Eigen::Vector3d{0, 0, 0}, rpyToRotationMatrix(0, 0, 0),
                  Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{-1, 0, 0}, gc_[12], gv_[11], gf_[11]);
    RF_leg.push_back(RF_SHANK);
    Body RF_shank_fixed(0.33742, getInertiaMatrix(0.00032748005, -2.142561e-05, 1.33942e-05, 0.00110974122, -7.601e-08, 0.00089388521),
                        Eigen::Vector3d{0.03463, -0.00688, 0.00098}, rpyToRotationMatrix(0, 0, 1.57079632679),
                        Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    RF_leg.push_back(RF_shank_fixed);
    Body RF_FOOT(0.25, getInertiaMatrix(0.00317174097, -2.63048e-06, 6.815581e-05, 0.00317174092, -6.815583e-05, 8.319196e-05),
                 Eigen::Vector3d{0.00948, 0.00948, 0.1468}, rpyToRotationMatrix(0, 0, 0),
                 Eigen::Vector3d{0.08795, -0.01305, -0.33797}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    RF_leg.push_back(RF_FOOT);

    /// LH leg
    std::vector<Body> LH_leg;
    LH_leg.push_back(LH_HAA);
    Body LH_HIP(0.001, getInertiaMatrix(0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001),
                Eigen::Vector3d{0, 0, 0}, rpyToRotationMatrix(0, 0, 0),
                Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{-1, 0, 0}, gc_[13], gv_[12], gf_[12]);
    LH_leg.push_back(LH_HIP);
    Body LH_hip_fixed(0.74, getInertiaMatrix(0.001393106, -8.4012e-05, -2.3378e-05, 0.003798579, 7.1319e-05, 0.003897509),
                      Eigen::Vector3d{-0.048, 0.008, -0.003}, rpyToRotationMatrix(-2.61799387799, 0, -3.14159265359),
                      Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    LH_leg.push_back(LH_hip_fixed);
    Body LH_HFE(2.04, getInertiaMatrix(0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827),
                Eigen::Vector3d{-0.063, 7e-05, 0.00046}, rpyToRotationMatrix(0, 0, 1.57079632679),
                Eigen::Vector3d{-0.0599, 0.08381, 0.0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    LH_leg.push_back(LH_HFE);

    Body LH_THIGH(0.001, getInertiaMatrix(0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001),
                  Eigen::Vector3d{0, 0, 0}, rpyToRotationMatrix(0, 0, 0),
                  Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{1, 0, 0}, gc_[14], gv_[13], gf_[13]);
    LH_leg.push_back(LH_THIGH);
    Body LH_thigh_fixed(1.03, getInertiaMatrix(0.018644469, -5.2e-08, -1.0157e-05, 0.019312599, 0.002520077, 0.002838361),
                        Eigen::Vector3d{-0.0, 0.018, -0.169}, rpyToRotationMatrix(0, 0, -1.57079632679),
                        Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    LH_leg.push_back(LH_thigh_fixed);
    Body LH_KFE(2.04, getInertiaMatrix(0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827),
                Eigen::Vector3d{-0.063, 7e-05, 0.00046}, rpyToRotationMatrix(0, 0, 1.57079632679),
                Eigen::Vector3d{-0.0, 0.1003, -0.285}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    LH_leg.push_back(LH_KFE);

    Body LH_SHANK(0.001, getInertiaMatrix(0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001),
                  Eigen::Vector3d{0, 0, 0}, rpyToRotationMatrix(0, 0, 0),
                  Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{1, 0, 0}, gc_[15], gv_[14], gf_[14]);
    LH_leg.push_back(LH_SHANK);
    Body LH_shank_fixed(0.33742, getInertiaMatrix(0.00032748005, -2.142561e-05, -1.33942e-05, 0.00110974122, 7.601e-08, 0.00089388521),
                        Eigen::Vector3d{-0.03463, 0.00688, 0.00098}, rpyToRotationMatrix(0, 0, -1.57079632679),
                        Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    LH_leg.push_back(LH_shank_fixed);
    Body LH_FOOT(0.25, getInertiaMatrix(0.00317174097, -2.63048e-06, -6.815581e-05, 0.00317174092, 6.815583e-05, 8.319196e-05),
                 Eigen::Vector3d{-0.00948, -0.00948, 0.1468}, rpyToRotationMatrix(0, 0, 0),
                 Eigen::Vector3d{-0.08795, 0.01305, -0.33797}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    LH_leg.push_back(LH_FOOT);

    /// RH leg
    std::vector<Body> RH_leg;
    RH_leg.push_back(RH_HAA);
    Body RH_HIP(0.001, getInertiaMatrix(0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001),
                Eigen::Vector3d{0, 0, 0}, rpyToRotationMatrix(0, 0, 0),
                Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{-1, 0, 0}, gc_[16], gv_[15], gf_[15]);
    RH_leg.push_back(RH_HIP);
    Body RH_hip_fixed(0.74, getInertiaMatrix(0.001393106, 8.4012e-05, -2.3378e-05, 0.003798579, -7.1319e-05, 0.003897509),
                      Eigen::Vector3d{-0.048, -0.008, -0.003}, rpyToRotationMatrix(2.61799387799, 0, -3.14159265359),
                      Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    RH_leg.push_back(RH_hip_fixed);
    Body RH_HFE(2.04, getInertiaMatrix(0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827),
                Eigen::Vector3d{-0.063, 7e-05, 0.00046}, rpyToRotationMatrix(0, 0, -1.57079632679),
                Eigen::Vector3d{-0.0599, -0.08381, 0.0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    RH_leg.push_back(RH_HFE);

    Body RH_THIGH(0.001, getInertiaMatrix(0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001),
                  Eigen::Vector3d{0, 0, 0}, rpyToRotationMatrix(0, 0, 0),
                  Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{-1, 0, 0}, gc_[17], gv_[16], gf_[16]);
    RH_leg.push_back(RH_THIGH);
    Body RH_thigh_fixed(1.03, getInertiaMatrix(0.018644469, 5.2e-08, -1.0157e-05, 0.019312599, -0.002520077, 0.002838361),
                        Eigen::Vector3d{-0.0, -0.018, -0.169}, rpyToRotationMatrix(0, 0, 1.57079632679),
                        Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    RH_leg.push_back(RH_thigh_fixed);
    Body RH_KFE(2.04, getInertiaMatrix(0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827),
                Eigen::Vector3d{-0.063, 7e-05, 0.00046}, rpyToRotationMatrix(0, 0, -1.57079632679),
                Eigen::Vector3d{-0.0, -0.1003, -0.285}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    RH_leg.push_back(RH_KFE);

    Body RH_SHANK(0.001, getInertiaMatrix(0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001),
                  Eigen::Vector3d{0, 0, 0}, rpyToRotationMatrix(0, 0, 0),
                  Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{-1, 0, 0}, gc_[18], gv_[17], gf_[17]);
    RH_leg.push_back(RH_SHANK);
    Body RH_shank_fixed(0.33742, getInertiaMatrix(0.00032748005, 2.142561e-05, -1.33942e-05, 0.00110974122, -7.601e-08, 0.00089388521),
                        Eigen::Vector3d{-0.03463, -0.00688, 0.00098}, rpyToRotationMatrix(0, 0, 1.57079632679),
                        Eigen::Vector3d{0, 0, 0}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    RH_leg.push_back(RH_shank_fixed);
    Body RH_FOOT(0.25, getInertiaMatrix(0.00317174097, 2.63048e-06, -6.815581e-05, 0.00317174092, -6.815583e-05, 8.319196e-05),
                 Eigen::Vector3d{-0.00948, 0.00948, 0.1468}, rpyToRotationMatrix(0, 0, 0),
                 Eigen::Vector3d{-0.08795, -0.01305, -0.33797}, Eigen::Vector3d{0, 0, 0}, 0, 0, 0);
    RH_leg.push_back(RH_FOOT);

    /// Forward Kinematics
    base_bodies[0].jointPosW = gc_.head(3);
    base_bodies[0].jointRotW = quatToRotMat(gc_.segment(3,4));
    getForwardKinematicsBase(base_bodies);
    getForwardKinematics(LF_leg, gc_);
    getForwardKinematics(RF_leg, gc_);
    getForwardKinematics(LH_leg, gc_);
    getForwardKinematics(RH_leg, gc_);

    /// Composite
    Body base_Composite1 = getCompositeBody(base_bodies[0], base_bodies[1]);
    Body base_Composite2 = getCompositeBody(base_Composite1, base_bodies[2]);
    Body base_Composite3 = getCompositeBody(base_Composite2, base_bodies[3]);
    Body base_Composite4 = getCompositeBody(base_Composite3, base_bodies[4]);
    Body lidar_Composite = getCompositeBody(base_bodies[5], base_bodies[6]);
    Body base_Composite5 = getCompositeBody(base_Composite4, lidar_Composite);
    Body base_Composite6 = getCompositeBody(base_Composite5, LF_leg[0]);
    Body base_Composite7 = getCompositeBody(base_Composite6, RF_leg[0]);
    Body base_Composite8 = getCompositeBody(base_Composite7, LH_leg[0]);
    Body base_Composite9 = getCompositeBody(base_Composite8, RH_leg[0]);
    Body base_Composite = getCompositeBody(base_Composite9, base_bodies[7]);

    std::vector<Body> LF_leg_bodies = getConmpositeLeg(LF_leg);
    std::vector<Body> RF_leg_bodies = getConmpositeLeg(RF_leg);
    std::vector<Body> LH_leg_bodies = getConmpositeLeg(LH_leg);
    std::vector<Body> RH_leg_bodies = getConmpositeLeg(RH_leg);

    /// Compute accelerations
    base_Composite.linVelW = gv_.segment(0,3);
    base_Composite.angVelW = gv_.segment(3,3);
    base_Composite.genVelW = gv_.segment(0, 6);

    getForwardVelocities(LF_leg_bodies, base_Composite);
    getForwardVelocities(RF_leg_bodies, base_Composite);
    getForwardVelocities(LH_leg_bodies, base_Composite);
    getForwardVelocities(RH_leg_bodies, base_Composite);

    /// Compute Articulated Mass and b
    getArticulatedMandB(LF_leg_bodies);
    getArticulatedMandB(RF_leg_bodies);
    getArticulatedMandB(LH_leg_bodies);
    getArticulatedMandB(RH_leg_bodies);

    std::vector<Body> thigh_bodies;
    thigh_bodies.push_back(LF_leg_bodies[0]); thigh_bodies.push_back(RF_leg_bodies[0]);
    thigh_bodies.push_back(LH_leg_bodies[0]); thigh_bodies.push_back(RH_leg_bodies[0]);
    getArticulatedMandBforBase(base_Composite, thigh_bodies);

    /// Compute u_dot and a
    Eigen::VectorXd genAccW = Eigen::VectorXd::Ones(18);

    Eigen::VectorXd tempGenAccW(6); tempGenAccW.setZero();
    base_Composite.genAccW = tempGenAccW;
    base_Composite.genForceW = gf_.segment(0, 6);

    getForwardAccelerationsforBase(base_Composite, genAccW);
    genAccW[2] += -9.81;   // gravity

    getForwardAccelerations(LF_leg_bodies, base_Composite, genAccW, 0);
    getForwardAccelerations(RF_leg_bodies, base_Composite, genAccW, 3);
    getForwardAccelerations(LH_leg_bodies, base_Composite, genAccW, 6);
    getForwardAccelerations(RH_leg_bodies, base_Composite, genAccW, 9);

    return genAccW;
}


