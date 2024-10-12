include <vector>
#include <type_traits>

template<typename T_num, typename T_denom>
using T_div = decltype(std::declval<T_num>() / std::declval<T_denom>());

/** класс для работы со СЛАУ с трехдиагональной матрицей **/
template<typename T_mat, typename T_col, std::size_t N>
class Tridiagonal_System {
  private:
    std::array<T_mat, N> a{}, c{};
    std::array<T_mat, N + 1> b{};
    std::array<T_col, N + 1> column{};

  public:
    Tridiagonal_System(const std::array<T_mat, N>& a, const std::array<T_mat, N + 1>& b, const std::array<T_mat, N>& c, const std::array<T_col, N + 1>& column) : a{a}, b{b}, c{c}, column{coloumn};
    /** метод прогонки **/
    std::array<T_div<T_col, T_mat>, N> solve_tridig_mat_alg() {
      std::array<T_mat, N + 1> p{};
      std::array<T_div<T_col, T_mat>, N + 1> q{}, x{};
      p[1] = -c[0] / b[0];
      q[1] = coloumn[0] / b[0];
      for (std::size_t i = 1; i < N; i++) {
        p[i + 1] = -c[i] / (a[i] * p[i] + b[i]);
        q[i + 1] = (column[i] - a[i] * q[i]) / (a[i] * p[i] + b[i]);
      }
      x[N] = (column[N] - a[N] * q[N]) / (p[N] * a[N] + b[N]);
      for (std::size_t i = N - 1; i < -1; i--) {
        x[i] = x[i + 1] * p[i + 1] + q[i + 1];
      }
      return x;
    }
};

/** класс естественного кубического сплайна**/
template<typename T_x, typename T_y, std::size_t N>
class CubicSpline {
  private:
    std::array<T_div<T_y, T_col>, N> A{}, B{}, C{};
    std::array<T_x, N> h{};
    std::array<T_y, N + 1> u{};

    Tridiagonal_System<T_x, T_y, N - 2> create_system(std::array<T_x, N>& h, std::array<T_y, N + 1>& u) {
      std::array<T_x, N - 2> a{}, c{};
      std::array<T_x, N - 1> b{};
      std::array<T_y, N - 1> column{};
      for (std::size_t i = 2; i < N + 1; i++) {
        dvd_dif_1_0 = (u[i - 1] - u[i - 2]) / h[i - 2];
        dvd_dif_1_1 = (u[i] - u[i - 1]) / h[i - 1];
        column[i - 2] = 6 * (dvd_dif_1_1 - dvd_dif_1_0) / (h[i - 1] - h[i - 2]);
        b[i - 2] = 2;
        if (i < N) {
          c[i - 2] = h[i - 1] / (h[i - 2] + h[i - 1]);
        }
        if (i > 2) {
          a[i - 3] = h[i - 2] / (h[i - 2] - h[i - 1]);
        }
      }
      Tridiagonal_System<T_x, T_y, N - 2> system(a, b, c, column);
      return system;
    }

    public:
    CubicSpline(const std::array<T_x, N + 1> &x, const std::array<T_y, N + 1>& y) : u{y} {
      for (std::size_t i = 1; i < N + 1; i++) {
        h[i - 1] = x[i] - x[i - 1];
        A[i - 1] = u[i];
      }
      Tridiagonal_System<T_x, T_y, N - 2> system = create_system(h, u);
      std::array<T_div<T_col, T_mat>, N - 2> C = system.solve_tridig_mat_alg();

    }

    yType interpolate(const xType& x) const noexcept;
};
