//
// Created by Jemin Hwangbo on 2022/03/17.
//

#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

#include "exercise1_20243406.hpp"
#include "raisim/RaisimServer.hpp"


int main(int argc, char* argv[]) {
  
  std::cout << "Eigen version: " << EIGEN_WORLD_VERSION << "." << EIGEN_MAJOR_VERSION << "." << EIGEN_MINOR_VERSION << std::endl;
  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world); // visualization server
  world.addGround();

  // anymal
  auto anymal = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/anymal_c/urdf/anymal.urdf");
  anymal->setName("anymal");
  server.focusOn(anymal);

  // anymal configuration
  Eigen::VectorXd jointNominalConfig(anymal->getGeneralizedCoordinateDim());
  jointNominalConfig << 0, 0, 0.54, 1.0, 0.0, 0.0, 0.0, 0.03, 0.4, -0.8, -0.03, 0.4, -0.8, 0.03, -0.4, 0.8, -0.03, -0.4, 0.8;
  anymal->setGeneralizedCoordinate(jointNominalConfig);
//
//  raisim::Mat<3,3> joint1_ori, joint2_ori, joint3_ori;
//  raisim::Vec<3> base_pos, joint1_pos, joint2_pos, joint3_pos, foot_pos;
//
//  anymal->getFrameOrientation(anymal->getFrameIdxByName("LH_HAA"), joint1_ori);
//  anymal->getFrameOrientation(anymal->getFrameIdxByName("LH_HFE"), joint2_ori);
//  anymal->getFrameOrientation(anymal->getFrameIdxByName("LH_KFE"), joint3_ori);
//
//  anymal->getFramePosition(anymal->getFrameIdxByName("base"), base_pos);
//  anymal->getFramePosition(anymal->getFrameIdxByName("LH_HAA"), joint1_pos);
//  anymal->getFramePosition(anymal->getFrameIdxByName("LH_HFE"), joint2_pos);
//  anymal->getFramePosition(anymal->getFrameIdxByName("LH_KFE"), joint3_pos);
//  anymal->getFramePosition(anymal->getFrameIdxByName("LH_shank_fixed_LH_FOOT"), foot_pos);
//  std::cout << foot_pos << std::endl;

//  std::cout << "relative joint orientation" << std::endl;
//  std::cout << joint1_ori.e() <<std::endl<<std::endl;
//  std::cout << joint1_ori.e().transpose() * joint2_ori.e() <<std::endl<<std::endl;
//  std::cout << joint2_ori.e().transpose() * joint3_ori.e() <<std::endl<<std::endl;
//
//  std::cout << "relative joint position" << std::endl;
//  std::cout << joint1_pos.e() - base_pos.e() <<std::endl<<std::endl;
//  std::cout << joint2_pos.e() - joint1_pos.e() <<std::endl<<std::endl;
//  std::cout << joint3_pos.e() - joint2_pos.e() <<std::endl<<std::endl;

  anymal->updateKinematics();

  // debug sphere
  auto debugSphere = server.addVisualSphere("debug_sphere", 0.04);
  debugSphere->setColor(1,0,0,1);
  debugSphere->setPosition(getEndEffectorPosition(jointNominalConfig));

  // solution sphere
  auto answerSphere = server.addVisualSphere("answer_sphere", 0.02);
  answerSphere->setColor(0,1,0,1);
  raisim::Vec<3> pos;
  anymal->getFramePosition("LH_shank_fixed_LH_FOOT", pos);
  answerSphere->setPosition(pos.e());

  // visualization
  server.launchServer(1234);
  for (int i=0; i<2000000; i++)
    std::this_thread::sleep_for(std::chrono::microseconds(1000));

  server.killServer();
}
