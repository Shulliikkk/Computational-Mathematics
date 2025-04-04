#include <iostream>
#include <utility>
#include <array>
#include <stdexcept>

template<typename T, std::size_t N>
struct System {
  std::array<std::array<T, N + 1>, N + 1> matrix{};
  std::array<T, N + 1> col{};
};

template<typename T, unsigned int N>
struct Derivative_coef {
  T central_coef;
  std::array<T, N> other_coefs{};
};

template<std::size_t... Ints>
constexpr std::size_t mult(std::index_sequence<Ints...>) noexcept {
  return ((Ints + 1) * ...);
}

template<std::size_t I>
constexpr std::size_t fact() noexcept {
  if constexpr (I == 0) {
    return 1;
  }
  else {
    return mult(std::make_index_sequence<I>());
  }
}

template<typename T, std::size_t N>
std::array<std::array<T, N - 1>, N - 1> minor(const std::array<std::array<T, N>, N>& matrix, const std::size_t i, const std::size_t j) noexcept {
  std::array<std::array<T, N - 1>, N - 1> minor_res;
  std::size_t row, col;
  for (std::size_t k = 0; k < N; k++) {
    for (std::size_t m = 0; m < N; m++) {
      if (k == i || m == j) {
        continue;
      }
      row = (k > i) ? k - 1 : k;
      col = (m > j) ? m - 1 : m;
      minor_res[row][col] = matrix[k][m];
    }
  }
  return minor_res;
}

template<typename T, std::size_t N>
T determinant(const std::array<std::array<T, N>, N>& matrix) noexcept {
  T det_res = 0;
  const std::size_t i = 0;
  if constexpr (N == 1) {
    return matrix[0][0];
  }
  else {
    for (std::size_t j = 0; j < N; j++) {
      det_res += ((i + j) % 2 == 0 ? 1 : -1) * matrix[i][j] * determinant<T, N - 1>(minor<T, N>(matrix, i, j));
    }
  }
  return det_res;
}

template<typename T, size_t N>
T det_change_col(std::array<std::array<T, N>, N> matrix, const std::array<T, N>& col, const std::size_t j) noexcept {
  for (std::size_t i = 0; i < N; i++) {
    matrix[i][j] = col[i];
  }
  return determinant<T, N>(matrix);
}

template<typename T, std::size_t N, std::size_t L>
System<T, N> create_system(std::array<T, N>& points) noexcept {
  static_assert(N > L, "N must be greater than L");
  System<T, N> system;
  for (std::size_t j = 0; j < N + 1; j++) {
    system.matrix[j][0] = (j == 0) ? 1 : 0;
  }
  for (std::size_t i = 1; i < N + 1; i++) {
    system.matrix[0][i] = 1;
  }
  for (std::size_t i = 1; i < N + 1; i++) {
    system.matrix[1][i] = points[i - 1];
  }
  for (std::size_t i = 1; i < N + 1; i++) {
  	for (std::size_t j = 2;  j < N + 1; j++) {
      system.matrix[j][i] = points[i - 1] * system.matrix[j - 1][i];
    }
  }
  system.col[L] = fact<L>();
  return system;
}

template<typename T, std::size_t N>
std::array<T, N> solve(const System<T, N>& system) {
  std::array<T, N> res;
  const T det = determinant<T, N + 1>(system.matrix);
  if (det == 0) {
    throw std::runtime_error("Determinant is zero");
  }
  for (std::size_t i = 0; i < N + 1; i++) {
    res[i] = det_change_col<T, N + 1>(system.matrix, system.col, i) / det;
  }
  return res;
}

template<typename T, std::size_t N, std::size_t L>
Derivative_coef<T, N> calc_derivative_coef(std::array<T, N> points) {
  Derivative_coef<T, N> derivative_coef;
  System<T, N> system = create_system<T, N, L>(points);
  std::array<T, N> solution = solve<T, N>(system);
  derivative_coef.central_coef = solution[0];
  for (unsigned int i = 0; i < N; i++) {
    derivative_coef.other_coefs[i] = solution[i + 1];
  }
  return derivative_coef;
}

int main(){
  try {
    const std::size_t N = 3, L = 2;
    std::array<double, N> points = std::array<double, N>{1, 2, 3};
    Derivative_coef<double, N> derivative_coef = calc_derivative_coef<double, N, L>(points);
    std::cout << "central_coef = " << derivative_coef.central_coef << std::endl;
    std::cout << "other_coefs = {";
    for (std::size_t i = 0; i < N - 1; i++) {
        std::cout << derivative_coef.other_coefs[i] << ", ";
    }
    std::cout << derivative_coef.other_coefs[N - 1] << '}';
  }
  catch (const std::runtime_error& e) {
    std::cerr << "Eror: " << e.what() << std::endl;
  }
}
