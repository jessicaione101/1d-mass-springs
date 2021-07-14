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
Eigen::VectorXd f(Eigen::VectorXd& y, FORCE &force) {
  Eigen::VectorXd y0(1);
  y0 << y(0);
  
  Eigen::VectorXd y1(1);
  y1 << y(1);
  
  Eigen::VectorXd _force;
  force(_force, y1, y0);
  
  Eigen::VectorXd _f(2);
  _f << _force(0), y0(0);
  
  return std::move(_f);
}

template<typename FORCE> 
inline void runge_kutta(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass, FORCE &force) {
  Eigen::MatrixXd A_inv(2, 2);
  A_inv(0, 0) = 1/mass;
  A_inv(0, 1) = 0;
  A_inv(1, 0) = 0;
  A_inv(1, 1) = 1;
  
  Eigen::VectorXd y(2);
  y << qdot(0), q(0);
  
  Eigen::VectorXd k1 = A_inv * f(y, force);
  
  Eigen::VectorXd _y(2);
  
  _y = y + dt/2 * k1;
  Eigen::VectorXd k2 = A_inv * f(_y, force);
  
  _y = y + dt/2 * k2;
  Eigen::VectorXd k3 = A_inv * f(_y, force);
  
  _y = y + dt * k3;
  Eigen::VectorXd k4 = A_inv * f(_y, force);
  
  _y = y + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
  
  q(0) = _y(1);
  qdot(0) = _y(0);
}
