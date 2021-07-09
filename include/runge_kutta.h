//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//Output:
//  q - set q to the updated generalized coordinate using Runge-Kutta time integration
//  qdot - set qdot to the updated generalized velocity using Runge-Kutta time integration
#include <iostream>
template<typename FORCE> 
inline void runge_kutta(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force) {
  Eigen::MatrixXd A_inv(2, 2);
  A_inv(0, 0) = 1/mass;
  A_inv(0, 1) = 0;
  A_inv(1, 0) = 0;
  A_inv(1, 1) = 1;
  
  Eigen::MatrixXd y(2, 1);
  y(0, 0) = qdot(0);
  y(1, 0) = q(0);
  
  Eigen::VectorXd f;
  force(f, q, qdot);
  
  
  
  
  
  Eigen::MatrixXd f_y(2, 1);
  f_y(0, 0) = f(0);
  f_y(1, 0) = qdot(0);
  
  
  
  Eigen::MatrixXd y_next(2, 1);
  
  
  
  static bool once = true;
  if (once) {
    std::cout << "---" << std::endl;
    std::cout << qdot << std::endl;
    std::cout << "---" << std::endl;
    std::cout << q << std::endl;
    std::cout << "---" << std::endl;
    std::cout << y << std::endl;
    std::cout << "---" << std::endl;
    once = false;
  }
}
