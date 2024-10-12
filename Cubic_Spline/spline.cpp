#include <array>
#include <iostream>

/** класс для работы со СЛАУ с трехдиагональной матрицей **/
template<typename T, std::size_t N>
class Tridiagonal_System {
  private:
    std::array<T, N - 1> a{}, c{};
    std::array<T, N> b{}, column{};

  public:
    Tridiagonal_System(const std::array<T, N - 1>& a, const std::array<T, N>& b, const std::array<T, N - 1>& c, const std::array<T, N>& column) : a{a}, b{b}, c{c}, column{column} {};
    /** метод прогонки **/
    std::array<T, N> solve_tridig_mat_alg() {
      std::array<T, N - 1> p{}, q{};
      std::array<T, N> x{};
      p[0] = -c[0] / b[0];
      q[0] = column[0] / b[0];
      for (std::size_t i = 1; i < N - 1; i++) {
        p[i] = -c[i] / (a[i] * p[i - 1] + b[i]);
        q[i] = (column[i] - a[i] * q[i - 1]) / (a[i] * p[i - 1] + b[i]);
      }
      x[N - 1] = (column[N - 1] - a[N - 1] * q[N - 2]) / (p[N - 2] * a[N - 1] + b[N - 1]);
      for (std::size_t i = N - 2; i < -1; i--) {
        x[i] = x[i + 1] * p[i] + q[i];
      }
      return x;
    }
};

/** класс естественного кубического сплайна **/
template<typename T, std::size_t N>
class Cubic_Spline {
  private:
    std::array<T, N> A{}, B{}, C{}, D{};
    std::array<T, N> h{};
    std::array<T, N + 1> u{}, x{};

    Tridiagonal_System<T, N - 1> create_system(std::array<T, N>& h, std::array<T, N + 1>& u) {
      std::array<T, N - 2> a{}, c{};
      std::array<T, N - 1> b{}, column{};;
      T dvd_dif_1_0, dvd_dif_1_1;
      for (std::size_t i = 2; i < N + 1; i++) {
        dvd_dif_1_0 = (u[i - 1] - u[i - 2]) / h[i - 2];
        dvd_dif_1_1 = (u[i] - u[i - 1]) / h[i - 1];
        column[i - 2] = 6 * (dvd_dif_1_1 - dvd_dif_1_0) / (h[i - 1] + h[i - 2]);
        b[i - 2] = 2;
        if (i < N) {
          c[i - 2] = h[i - 1] / (h[i - 2] + h[i - 1]);
        }
        if (i > 2) {
          a[i - 2] = h[i - 2] / (h[i - 2] + h[i - 1]);
        }
      }
      Tridiagonal_System<T, N - 1> system(a, b, c, column);
      return system;
    }

    public:
    Cubic_Spline(const std::array<T, N + 1> &x, const std::array<T, N + 1>& y) :x{x}, u{y} {
      for (std::size_t i = 1; i < N + 1; i++) {
        h[i - 1] = x[i] - x[i - 1];
        A[i - 1] = u[i];
      }
      Tridiagonal_System<T, N - 1> system = create_system(h, u);
      std::array<T, N - 1> C_0 = system.solve_tridig_mat_alg();
      for (std::size_t i = 0; i < N - 1; i++) {
        C[i] = C_0[i];
      }
      C[N - 1] = 0;
      D[0] = C[0] / h[0];
      for (std::size_t i = 1; i < N; i++) {
        D[i] = (C[i] - C[i - 1]) / h[i];
      }
      B[0] = A[0] / h[0] + C[0] / 2 * h[0] - D[0] / 6 * h[0] * h[0] - u[0] / h[0];
      for (std::size_t i = 1; i < N; i++) {
        B[i] = B[i - 1] + C[i] * h[i] - D[i] / 2 * h[i] * h[i];
      }
    }

    T interpolate(const T& x_0) const noexcept {
      T y_0;
      for (int i = 0; i < N; i++) {
        if (x_0 >= x[i] && x_0 <= x[i + 1]) {
          y_0 = A[i] + B[i] * (x_0 - x[i + 1]) + C[i] / 2 * (x_0 - x[i + 1]) * (x_0 - x[i + 1]) + D[i] / 6 * (x_0 - x[i + 1]) * (x_0 - x[i + 1]) * (x_0 - x[i + 1]);
          break;
        }
      }
      return y_0;
    }
};

int main() {
  std::array<double, 5> x{1, 2, 3, 4, 5};
  std::array<double, 5> y{1, 4, 9, 16, 25};
  Cubic_Spline<double, 4> spline(x, y);
  double x_0 = 3.5;
  std::cout << "Интерполяция:" << spline.interpolate(x_0) << '\n';
  std::cout << "Точное значение:" << x_0 * x_0 << '\n';
}
