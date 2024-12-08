#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <../src/odeint.hpp>

std::array<double, 4> f(double t, std::array<double, 4> x) {
  const double m = 0.012277471, M = 1 - m;
  double R_1 = std::pow(((x[0] + m) * (x[0] + m) + x[2] * x[2]), 3. / 2.),
         R_2 = std::pow(((x[0] - M) * (x[0] - M) + x[2] * x[2]), 3. / 2.);
  std::array<double, 4> rhs{x[1],
                            x[0] + 2 * x[3] - M * (x[0] + m) / R_1 - m * (x[0] - M) / R_2,
                            x[3],
                            x[2] - 2 * x[1] - M * x[2] / R_1 - m * x[2] / R_2};
  return rhs;
}

int main() {
  std::ofstream out("orbit.csv");
  std::array<double, 4> init_condition{0.994, 0.0, 0.0, -2.031732629557337};
  std::vector<Point<4>> sol = DP45<4>(f, init_condition, 11.124340337, 0.0001, 1e-15);
  out << 'x' << ',' << 'y' << std::endl;
  for (std::size_t i = 0; i < sol.size(); i++) {
    out << sol[i].state[0] << ',' << sol[i].state[2] << std::endl;
  }
}
