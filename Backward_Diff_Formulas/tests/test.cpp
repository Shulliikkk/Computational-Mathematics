#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <../src/BDF.hpp>

Eigen::Vector<double, 2> f(double t, Eigen::Vector<double, 2> x) {
  return Eigen::Vector<double, 2>(x[1], -x[0]);
}

int main() {
  Eigen::Vector<double, 2> init_state(0.0, 1.0);
  std::vector<Point<2>> sol = BDF<2>(f, init_state, 5., 1e-4, 1e-15, 100);
  std::cout << sol[sol.size() - 1].state[0] << std::endl;
}
