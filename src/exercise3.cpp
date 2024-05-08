//
// Created by Jemin Hwangbo on 2022/04/08.
//

#include "raisim/RaisimServer.hpp"
#include "exercise3_STUDENTID.hpp"

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
  gv << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8;
  anymal->setState(gc, gv);

  /// if you are using an old version of Raisim, you need this line
  world.integrate1();
  
  std::vector<raisim::Vec<3>> &comVec = anymal->getBodyCOM_W();
  std::vector<raisim::Mat<3,3>> &inertiaMat = anymal->getInertia();
  std::vector<double> &mass = anymal->getMass();

  for (size_t i = 0; i < comVec.size(); ++i) {
    std::cout << "Body " << i << "의 질량 중심 위치와 질량 : ";
    std::cout << "(" << comVec[i][0] << ", " << comVec[i][1] << ", " << comVec[i][2] << ")" << "    " << mass[i] << std::endl;
    std::cout << "Inertia matrix" << std::endl;
    std::cout << inertiaMat[i].e() << std::endl;
  }
//int i=0;
//  std::cout << "Body " << i << "의 질량 중심 위치와 질량 : ";
//  std::cout << "(" << comVec[i][0] << ", " << comVec[i][1] << ", " << comVec[i][2] << ")" << "    " << mass[i] << std::endl;
//  std::cout << "Inertia matrix" << std::endl;
//  std::cout << inertiaMat[i].e() << std::endl;
  
  raisim::Vec<3> pos;
  raisim::Mat<3,3> rot;
  anymal->getFramePosition("RH_shank_fixed_RH_FOOT", pos);
  anymal->getFrameOrientation("RH_shank_fixed_RH_FOOT", rot);
  std::cout << "RH_shank_fixed_RH_FOOT pos " << pos.e().transpose() << std::endl;
  std::cout << "RH_shank_fixed_RH_FOOT rot \n" << rot.e() << std::endl;

//  std::cout<<"mass matrix should be \n"<< anymal->getMassMatrix().e()<<std::endl;

  if((getMassMatrix(gc) - anymal->getMassMatrix().e()).norm() < 1e-8)
    std::cout<<"passed "<<std::endl;
  else
    std::cout<<"failed "<<std::endl;  

  return 0;
}
