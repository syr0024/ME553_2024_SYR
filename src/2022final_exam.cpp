//
// Created by jemin on 6/15/22.
//


#include "raisim/World.hpp"
#include "2022final_exam_STUDENTID.hpp"
#include "raisim/RaisimServer.hpp"

#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

int main(int argc, char* argv[]) {
  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world);
  world.setTimeStep(0.002);
  auto cartpole_double = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/cartPole/cartpole_double.urdf");

  // a1_simplified configuration
  Eigen::VectorXd gc(cartpole_double->getGeneralizedCoordinateDim()), gv(cartpole_double->getDOF());
  gc << 0, 1, 2;
  gv << 0.1, 0.2, 0.3;
  cartpole_double->setState(gc, gv);

  /// this function updates internal variables in raisim (such as the mass matrix and nonlinearities)
  world.integrate1();
  cartpole_double->getMassMatrix(); /// due to the current bug in RaiSim, we have to compute the mass matrix first to get the bias term
  
  
//  raisim::Vec<3> pos;
//  cartpole_double->getFramePosition("world_to_sliderBar", pos);
//  std::cout << "pos " << pos.e().transpose() << std::endl;
//  cartpole_double->getFramePosition("slider", pos);
//  std::cout << "pos " << pos.e().transpose() << std::endl;
//  cartpole_double->getFramePosition("rod_revolute", pos);
//  std::cout << "pos " << pos.e().transpose() << std::endl;
//  cartpole_double->getFramePosition("rod_revolute2", pos);
//  std::cout << "pos " << pos.e().transpose() << std::endl;

//  raisim::Mat<3,3> ori;
//  cartpole_double->getFrameOrientation("world_to_sliderBar", ori);
//  std::cout << "ori\n " << ori.e() << std::endl;
//  cartpole_double->getFrameOrientation("slider", ori);
//  std::cout << "ori\n " << ori.e() << std::endl;
//  cartpole_double->getFrameOrientation("rod_revolute", ori);
//  std::cout << "ori\n " << ori.e() << std::endl;
//  cartpole_double->getFrameOrientation("rod_revolute2", ori);
//  std::cout << "ori\n " << ori.e() << std::endl;

//  raisim::Vec<3> angVel;
//  cartpole_double->getFrameAngularVelocity("slider", angVel);
//  std::cout << "angVel " << angVel.e().transpose() << std::endl;
//  cartpole_double->getFrameAngularVelocity("rod_revolute", angVel);
//  std::cout << "angVel " << angVel.e().transpose() << std::endl;
//  cartpole_double->getFrameAngularVelocity("rod_revolute2", angVel);
//  std::cout << "angVel " << angVel.e().transpose() << std::endl;

//  raisim::Vec<3> linVel;
//  cartpole_double->getFrameVelocity("slider", linVel);
//  std::cout << "angVel " << linVel.e().transpose() << std::endl;
//  cartpole_double->getFrameVelocity("rod_revolute", linVel);
//  std::cout << "angVel " << linVel.e().transpose() << std::endl;
//  cartpole_double->getFrameVelocity("rod_revolute2", linVel);
//  std::cout << "angVel " << linVel.e().transpose() << std::endl;

//  raisim::Vec<3> linAcc;
//  cartpole_double->getFrameAcceleration("slider", linAcc);
//  std::cout << "\nlinAcc " << linAcc.e().transpose() << std::endl;
//  cartpole_double->getFrameAcceleration("rod_revolute", linAcc);
//  std::cout << "linAcc " << linAcc.e().transpose() << std::endl;
//  cartpole_double->getFrameAcceleration("rod_revolute2", linAcc);
//  std::cout << "linAcc " << linAcc.e().transpose() << std::endl;

//  std::vector<raisim::Vec<3>> comVec = cartpole_double->getCompositeCOM(); //getBodyCOM_W(); // cartpole_double->getCompositeCOM()
//  std::vector<raisim::Mat<3,3>> inertia = cartpole_double->getCompositeInertia();
//  std::vector<double> mass = cartpole_double->getCompositeMass(); // cartpole_double->getCompositeMass()
//
//   for (size_t i = 1; i < comVec.size(); ++i) {
//     std::cout << "Body " << i << " COM and mass : ";
//     std::cout << "(" << comVec[i][0] << ", " << comVec[i][1] << ", " << comVec[i][2] << ")" << "    " << mass[i] << std::endl;
//     std::cout << "Body " << i << " inertia : \n";
//     std::cout << inertia[i].e() << std::endl;
//   }



//  std::cout<<"mass matrix should be \n"<< cartpole_double->getMassMatrix().e()<<std::endl;
  std::cout<<"nonlinearities should be \n"<< cartpole_double->getNonlinearities({0,0,-9.81}).e().transpose() <<std::endl;
//  std::cout<<"The acceleration should be \n"<< (massMatrix.inverse() * (gf-nonlinearity)).transpose() <<std::endl;

  
  if((getNonlinearities(gc, gv) - cartpole_double->getNonlinearities({0,0,-9.81}).e()).norm() < 1e-8)
    std::cout<<"passed "<<std::endl;
  else
    std::cout<<"failed "<<std::endl;

  // uncomment these lines if you want to visualize the system
//  server.launchServer();
//
//  while(true) {
//    server.integrateWorldThreadSafe();
//    std::this_thread::sleep_for(std::chrono::microseconds(int(world.getTimeStep() * 1e6)));
//  }

  return 0;
}
