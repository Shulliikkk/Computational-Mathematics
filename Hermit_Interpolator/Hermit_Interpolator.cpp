#include <array>
#include <iostream>

template<typename T_x, typename T_y, std::size_t N>
class Hermit_Interpolator {
  private:
    std::array<T_x, 2 * N> z{};
    std::array<T_y, 2 * N> dvd_dif{};

    std::pair<std::array<T_x, 2 * N>,std::array<T_y, 2 * N>> stretch(const std::array<T_x, N>& x, const std::array<T_y, N>& f) noexcept {
      std::array<T_x, 2 * N> z{};
      std::array<T_y, 2 * N> dvd_dif{};
      for (std::size_t i = 0; i < N; i++) {
        z[2 * i] = x[i];
        z[2 * i + 1] = x[i];
        dvd_dif[2 * i] = f[i];
        dvd_dif[2 * i + 1] = f[i];
      }
      return std::pair<std::array<T_x, 2 * N>, std::array<T_y, 2 * N>>(z, dvd_dif);
    }
  public:
    Hermit_Interpolator(const std::array<T_x, N> &x, const std::array<T_y, N>& f, const std::array<T_y, N>& df) noexcept {
      std::pair<std::array<T_x, 2 * N>,std::array<T_y, 2 * N>> strch = stretch(x, f);
      z = strch.first;
      dvd_dif = strch.second;
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

    T_y interpolate(const T_x& x_0) const noexcept {
      T_y f_x_0 = dvd_dif[2 * N - 2] + dvd_dif[2 * N - 1] * (x_0 - z[2 * N - 2]);
      for (std::size_t i = 2 * N - 2; i > 0; i--) {
        f_x_0 *= x_0 - z[i - 1];
        f_x_0 += dvd_dif[i - 1];
      }
      return f_x_0;
    }
};

int main() {
  const std::size_t N = 3;
  std::array<double, N> x{0, 1, 2};
  std::array<double, N> f{1, 2, 33};
  std::array<double, N> df{0, 5, 80};
  const double x_0 = 2;
  const Hermit_Interpolator<double, double, N> intrpolator(x, f, df);
  std::cout << intrpolator.interpolate(x_0);
}
