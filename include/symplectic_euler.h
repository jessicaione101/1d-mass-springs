//force - function which computes the gradient of the potential energy function as a function of the state (q, qdot).
template<typename FORCE> 
inline void symplectic_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force) {

}