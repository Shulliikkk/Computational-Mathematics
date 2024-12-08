#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <../src/odeint.hpp>

std::array<double, 2> osc(const double t, const std::array<double, 2> x) {
  return std::array<double, 2>{x[1], -x[0]};
}

int main() {
  std::array<double, 2> init_state{0.0, 1.0};
  std::vector<Point<2>> sol = DP45<2>(osc, init_state, 5., 1e-4, 1e-15);
  std::cout.precision(15);
  std::cout << sol.size() << std::endl;
  std::cout << sol[sol.size() - 1].time << ' ' << sol[sol.size() - 1].state[0] << std::endl;
}
