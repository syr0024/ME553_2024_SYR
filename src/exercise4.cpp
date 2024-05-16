//
// Created by Jemin Hwangbo on 2024/05/12.
//

#include "raisim/RaisimServer.hpp"
#include "exercise4_20243406.hpp"

#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world

  // anymal
  auto anymal = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/anymal_c/urdf/anymal.urdf");

  // anymal configuration
  Eigen::VectorXd gc(anymal->getGeneralizedCoordinateDim()), gv(anymal->getDOF());
  gc << 0, 0, 0.54, 1.0, 0.0, 0.0, 0.0, 0.03, 0.4, -0.8, -0.03, 0.4, -0.8, 0.03, -0.4, 0.8, -0.03, -0.4, 0.8; /// Jemin: I'll randomize the gc, gv when grading
  gv << 0.1, 0.2, 0.3, 0.1, 0.4, 0.3, 0.2,0.2,0.2, 0.2,0.2,0.2, 0.2,0.2,0.2, 0.2,0.2,0.2;
  anymal->setState(gc, gv);

  /// if you are using an old version of Raisim, you need this line
  world.integrate1();
  anymal->getMassMatrix();

  std::cout<<"nonlinearities should be \n"<< anymal->getNonlinearities({0,0,-9.81}).e().transpose()<<std::endl;

  if((getNonlinearities(gc, gv) - anymal->getNonlinearities({0,0,-9.81}).e()).norm() < 1e-8)
    std::cout<<"passed "<<std::endl;
  else
    std::cout<<"failed "<<std::endl;
  
  
  raisim::Vec<3> pos;
  anymal->getFramePosition("LF_HAA", pos);
  std::cout << "pos " << pos.e().transpose() << std::endl;
  anymal->getFramePosition("LF_HFE", pos);
  std::cout << "pos " << pos.e().transpose() << std::endl;
  
  raisim::Mat<3,3> ori;
  anymal->getFrameOrientation("LF_HAA", ori);
  std::cout << "ori " << ori.e() << std::endl;
  
  raisim::Vec<3> angVel;
  anymal->getFrameAngularVelocity("LF_HAA", angVel);
  std::cout << "angVel " << angVel.e().transpose() << std::endl;
  anymal->getFrameAngularVelocity("LF_HFE", angVel);
  std::cout << "angVel " << angVel.e().transpose() << std::endl;
  anymal->getFrameAngularVelocity("RH_shank_fixed_RH_FOOT", angVel);
  std::cout << "angVel " << angVel.e().transpose() << std::endl;
  
  raisim::Vec<3> linAcc;
  anymal->getFrameAcceleration("LF_HAA", linAcc);
  std::cout << "linAcc " << linAcc.e().transpose() << std::endl;
  anymal->getFrameAcceleration("LF_HFE", linAcc);
  std::cout << "linAcc " << linAcc.e().transpose() << std::endl;
  anymal->getFrameAcceleration("LF_KFE", linAcc);
  std::cout << "linAcc " << linAcc.e().transpose() << std::endl;
  anymal->getFrameAcceleration("RF_HAA", linAcc);
  std::cout << "linAcc " << linAcc.e().transpose() << std::endl;
  anymal->getFrameAcceleration("RF_KFE", linAcc);
  std::cout << "linAcc " << linAcc.e().transpose() << std::endl;
  
  return 0;
}
