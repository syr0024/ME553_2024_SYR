//
// Created by Jemin Hwangbo on 2022/04/08.
//

#include "raisim/RaisimServer.hpp"
#include "exercise3_20243406.hpp"

#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world);

  // kinova
  auto anymal = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/anymal_c/urdf/anymal.urdf");

  // kinova configuration
  Eigen::VectorXd gc(anymal->getGeneralizedCoordinateDim()), gv(anymal->getDOF());
  gc << 0, 0, 0.54, 1.0, 0.0, 0.0, 0.0, 0.03, 0.4, -0.8, -0.03, 0.4, -0.8, 0.03, -0.4, 0.8, -0.03, -0.4, 0.8; /// Jemin: I'll randomize the gc, gv when grading
//  gc.head(3) << Eigen::Vector3d::Random(); gc(2) + 0.4;
//  Eigen::Vector4d quat = Eigen::Vector4d::Random();
//  gc.segment(3,4) << quat/quat.norm();
//  gc.segment(7, 12) << Eigen::VectorXd::Random(12);
  gv << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8;
  anymal->setState(gc, gv);

  /// if you are using an old version of Raisim, you need this line
  world.integrate1();
  
  std::cout<<"mass matrix should be \n"<< anymal->getMassMatrix().e()<<std::endl;
  
//  raisim::Vec<3> pos;
//  raisim::Mat<3,3> rot;
//  anymal->getFramePosition("RH_shank_fixed_RH_FOOT", pos);
//  anymal->getFrameOrientation("RH_shank_fixed_RH_FOOT", rot);
//  std::cout << "RH_shank_fixed_RH_FOOT pos " << pos.e().transpose() << std::endl;
//  std::cout << "RH_shank_fixed_RH_FOOT rot \n" << rot.e() << std::endl;
  
  if((getMassMatrix(gc) - anymal->getMassMatrix().e()).norm() < 1e-8)
    std::cout<<"passed "<<std::endl;
  else
    std::cout<<"failed "<<std::endl;

  return 0;
}
