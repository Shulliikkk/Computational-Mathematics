#include "Hermit.hpp"
#include <array>
#include <iostream>
#include <cmath>

template<typename T>
T func(T x) {
    return x * sin(x);
}

template<typename T>
T dfunc(T x) {
    return sin(x) + x * cos(x);
}

int main() {
  const std::size_t N = 4,  M = 6, N_err = 1000;
  double H, h, err;

  std::array<double, N> x{};
  std::array<double, N> f{};
  std::array<double, N> df{};
  std::array<double, M> Err;

  double start = 0.;
  double end = 2.;

  for (std::size_t i = 0; i < M; i++) {
    h = (end - start) / N;
    for (std::size_t j = 0; j < N; j++) {
      x[j] = start + j * h;
      f[j] = func(x[j]);
      df[j] = dfunc(x[j]);
    }

    Hermit_Interpolator<double, double, N> interpolator(x, f, df);

    err = 0;
    H = (end - start) / N_err;
    for (std::size_t k = 0; k < N_err; k++) {
      err = std::max(err, std::abs(func(start + k * H) - interpolator.interpolate(start + k * H)));
    }
    Err[i] = err;
    end /= 2;
  }

  std::cout << "Errors ";
  for (std::size_t i = 0; i < M; i++) {
    std::cout << Err[i] << ' ';
  }
  return 0;
}
