#include <array>
#include <iostream>

template<std::size_t N>
constexpr std::size_t dob() {
  return 2 * N;
}

template<typename T_x, typename T_y, std::size_t N>
class Ermit_Interpolator {
  private:
    const std::array<T_x, N> x;
    const std::array<T_y, N> f;
    const std::array<T_y, N> df;

    std::pair<std::array<T_x, dob<N>()>,std::array<T_y, dob<N>()>> stretch(const std::array<T_x, N>& x, const std::array<T_y, N>& f) const noexcept {
      std::array<T_x, dob<N>()> z{};
      std::array<T_y, dob<N>()> f_z{};
      for (std::size_t i = 0; i < N; i++) {
        z[2 * i] = x[i];
        z[2 * i + 1] = x[i];
        f_z[2 * i] = f[i];
        f_z[2 * i + 1] = f[i];
      }
      return std::pair<std::array<T_x, dob<N>()>,std::array<T_y, dob<N>()>>(z, f_z);
    }

    std::array<T_y, dob<N>()> divided_differences(const std::array<T_x, dob<N>()>& z, const std::array<T_y, dob<N>()>& f_z, const std::array<T_y, N>& d_f) const noexcept {
      std::array<T_y, dob<N>()> dv_df_z{};
      dv_df_z[0] = f_z[0];
      for (std::size_t j = 1; j < dob<N>(); j++) {
        if (z[j] == z[j - 1]) {
          dv_df_z[j] = d_f[(j - j % 2) / 2];
        }
        else {
          dv_df_z[j] = (f_z[j] - f_z[j - 1]) / (z[j] - z[j - 1]);
        }
      }
      for (std::size_t i = 2; i < dob<N>(); i++) {
        T_y dv_df_z_j_prev = dv_df_z[i - 1];
        for (std::size_t j = i; j < dob<N>(); j++) {
          T_y t = (dv_df_z[j] - dv_df_z_j_prev) / (z[j] - z[j - i]);
          dv_df_z_j_prev = dv_df_z[j];
          dv_df_z[j] = t;
        }
      }
      return dv_df_z;
    }
  public:
    Ermit_Interpolator(const std::array<T_x, N> &x, const std::array<T_y, N>& f, const std::array<T_y, N>& df) noexcept : x(x), f(f), df(df) {};

    T_y interpolate(const T_x& x_0) const noexcept {
      std::pair<std::array<T_x, dob<N>()>,std::array<T_y, dob<N>()>> strch = stretch(x, f);
      const std::array<T_x, dob<N>()> z = strch.first;
      const std::array<T_y, dob<N>()> f_z = strch.second;
      const std::array<T_y, dob<N>()> dv_df_z = divided_differences(z, f_z, df);
      T_y f_x_0 = dv_df_z[dob<N>() - 2] + dv_df_z[dob<N>() - 1] * (x_0 - z[dob<N>() - 2]);
      for (std::size_t i = dob<N>() - 2; i > 0; i--) {
        f_x_0 *= x_0 - z[i - 1];
        f_x_0 += dv_df_z[i - 1];
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
  Ermit_Interpolator<double, double, N> intrpolator(x, f, df);
  std::cout << intrpolator.interpolate(x_0);
}
