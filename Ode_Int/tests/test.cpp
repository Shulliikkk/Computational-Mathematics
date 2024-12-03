#include <array>
#include <vector>
#include <iostream>
#include <../src/odeint.hpp>

std::array<double, 2> f(double t, std::array<double, 2> x) {
  std::array<double, 2> rhs{x[1], - x[0]};
  return rhs;
}

int main() {
  std::array<double, 2> init_condition{0.0, 1.0};
  std::vector<Point<std::array<double, 2>>> sol = odeint<2>(f, init_condition, 5., 0.0001, 1e-12);
  std::cout << sol[sol.size() - 1].state[0] << std::endl;
  std::cout << sol.size();
}
