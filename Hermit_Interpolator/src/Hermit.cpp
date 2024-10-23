#include "Hermit.hpp"
#include <array>

template<typename T_x, typename T_y, std::size_t N>
std::pair<std::array<T_x, 2 * N>,std::array<T_y, 2 * N>> Hermit_Interpolator<T_x, T_y, N>::stretch(const std::array<T_x, N>& x, const std::array<T_y, N>& f) noexcept {
  std::array<T_x, 2 * N> z{};
  std::array<T_y, 2 * N> dvd_dif{};
  for (std::size_t i = 0; i < N; i++) {
    z[2 * i] = x[i];
    z[2 * i + 1] = x[i];
    dvd_dif[2 * i] = f[i];
    dvd_dif[2 * i + 1] = f[i];
  }
  return {z, dvd_dif};
}

template<typename T_x, typename T_y, std::size_t N>
Hermit_Interpolator<T_x, T_y, N>::Hermit_Interpolator(const std::array<T_x, N> &x, const std::array<T_y, N>& f, const std::array<T_y, N>& df) noexcept {
  std::pair<std::array<T_x, 2 * N>,std::array<T_y, 2 * N>> strch = stretch(x, f);
  z = strch.first;
  dvd_dif = strch.second;
  //auto [z, dvd_dif] = stretch(x, f);

  T_y t, dvd_dif_j_prev = dvd_dif[0];
  for (std::size_t j = 1; j < 2 * N; j++) {
    if (z[j] == z[j - 1]) {
      dvd_dif_j_prev = dvd_dif[j];
      dvd_dif[j] = df[(j - j % 2) / 2];
    }
    else {
      t = (dvd_dif[j] - dvd_dif_j_prev) / (z[j] - z[j - 1]);
      dvd_dif_j_prev = dvd_dif[j];
      dvd_dif[j] = t;
    }
  }
  for (std::size_t i = 2; i < 2 * N; i++) {
    dvd_dif_j_prev = dvd_dif[i - 1];
    for (std::size_t j = i; j < 2 * N; j++) {
      t = (dvd_dif[j] - dvd_dif_j_prev) / (z[j] - z[j - i]);
      dvd_dif_j_prev = dvd_dif[j];
      dvd_dif[j] = t;
    }
  }
}

template<typename T_x, typename T_y, std::size_t N>
T_y Hermit_Interpolator<T_x, T_y, N>::interpolate(const T_x& x_0) const noexcept {
  T_y f_x_0 = dvd_dif[2 * N - 2] + dvd_dif[2 * N - 1] * (x_0 - z[2 * N - 2]);
  for (std::size_t i = 2 * N - 2; i > 0; i--) {
    f_x_0 *= x_0 - z[i - 1];
    f_x_0 += dvd_dif[i - 1];
  }
  return f_x_0;
}
