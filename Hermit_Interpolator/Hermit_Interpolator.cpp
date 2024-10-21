#include <array>
#include <iostream>
#include <cmath>

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

template<typename T>
T func(T x) {
    return exp(x);
}

template<typename T>
T dfunc(T x) {
    return exp(x);
}

int main() {
  const std::size_t N = 5,  M = 6, N_err = 1000;
  double H, err;

  std::array<double, N> x{};
  std::array<double, N> f{};
  std::array<double, N> df{};
  std::array<double, M> Err;

  double start = 0.;
  double end = 2.;
  double h = (end - start) / 2.;
  for (std::size_t i = 0; i < M; i++) {
    for (std::size_t j = 0; j < N; j++) {
      x[j] = start + j * h;
      f[j] = func(x[j]);
      df[j] = dfunc(x[j]);
    }

    const Hermit_Interpolator<double, double, N> interpolator(x, f, df);

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
