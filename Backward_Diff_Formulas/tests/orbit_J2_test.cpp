#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <cmath>
#include <../src/BDF.hpp>

Eigen::Vector<double, 6> f(const double t, const Eigen::Vector<double, 6> x) {
  const double mu = 398600.4418, R_e = 6378.137, J_2 = 1.08262668e-3;
  const double eps = 3. * mu * R_e * R_e * J_2 / 2.;
  double r = std::sqrt(x[0] * x[0] + x[2] * x[2] + x[4] * x[4]);
  return Eigen::Vector<double, 6>(x[1],
                                  x[0] * (5 * eps * x[4] * x[4] / std::pow(r, 7) - eps / std::pow(r, 5) - mu / std::pow(r, 3)),
                                  x[3],
                                  x[2] * (5 * eps * x[4] * x[4] / std::pow(r, 7) - eps / std::pow(r, 5) - mu / std::pow(r, 3)),
                                  x[5],
                                  5 * eps * x[4] * x[4] / std::pow(r, 7) - 3 * eps * x[4] / std::pow(r, 5) - mu * x[4] / std::pow(r, 3));
}

int main() {
  std::ofstream out("orbit_J2.csv");
  const double mu = 398600.4418, R_e = 6378.137, H = 400;;
  Eigen::Vector<double, 6> init_state{R_e + H, 0.0, 0.0, std::sqrt(mu / (R_e + H)), 0.0, 0.0};
  std::vector<Point<6>> sol = BDF<6, BT_RK5, CF_BDF5>(f, init_state, 5. * 24. * 60. * 60., 1e-1);
  out << 'x' << ',' << 'y' << ',' << 'z' << std::endl;
  for (std::size_t i = 0; i < sol.size(); i++) {
    out << sol[i].state[0] << ',' << sol[i].state[2] << ',' << sol[i].state[4]  << std::endl;
  }
}
