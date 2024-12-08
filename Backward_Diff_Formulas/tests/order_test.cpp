#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include <../src/BDF.hpp>

Eigen::Vector<double, 2> F(const double t, const Eigen::Vector<double, 2> x) {
  return Eigen::Vector<double, 2>(x[1], -x[0]);
}

Eigen::Matrix<double, 2, 2> J(const double t, const Eigen::Vector<double, 2> x) {
  Eigen::Matrix<double, 2, 2> ans;
  ans << 0., 1., -1., 0;
  return ans;
}

int main() {
  std::ofstream out("BDF_order.csv");
  Eigen::Vector<double, 2> init_state(0.0, 1.0);
  double h = 0.1, max_err;
  out << "err" << ", " << "h" << std::endl;
  for (std::size_t i = 1; i < 20; i += 1) {
    //std::vector<Point<2>> sol = BDF<2, BT_RK5, CF_BDF5>(F, J, init_state, 5., h, 1e-15, 100);
    std::vector<Point<2>> sol = BDF<2, BT_RK5, CF_BDF5>(F, init_state, 5., h, 1e-15, 100);
    max_err = 0.0;
    for (std::size_t j = 1; j < sol.size(); j += 1) {
      max_err = std::max(max_err, std::abs(sol[i].state[0] - std::sin(i * h)));
    }
    out.precision(15);
    out << max_err << ", " << h << std::endl;
    h /= 2;
  }
}
