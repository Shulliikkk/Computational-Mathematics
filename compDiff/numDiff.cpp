#include <iostream>
#include <array>
#include <cmath>

template<typename RealType>
RealType fact(RealType L) {
  RealType fact = 1;
  for (unsigned int i = 2; i < L + 1; i++) {
    fact *= i;
  }
  return fact;
}

template<typename RealType, unsigned int N>
struct DerivativeCoef {
    RealType centralCoef;
    std::array<RealType, N> otherCoefs;
};

template<typename RealType, unsigned int N, unsigned int L>
DerivativeCoef<RealType, N> calcDerivativeCoef(const std::array<RealType, N>& points) noexcept {
  // создание и заполнения матрицы СЛАУ для метода неопределенных коэффициентов
  unsigned int i_0;

  if (points[0] > 0) {
    i_0 = 0;
  }
  else if (points[N - 1] < 0) {
    i_0 = N;
  }
  else {
    for (unsigned int i = 0; i < N - 1; i++) {
      if (points[i] < 0 && points[i + 1] > 0) { // points[i] + 1 == 0
        i_0 = i + 1;
        break;
      }
    }
  }

  std::array<std::array<RealType, N + 2>, N + 1> A = {};
  for (unsigned int i = 0; i < N + 1; i++) {
		for (unsigned int j = 0;  j < N + 1; j++) {
      if (i < i_0) {
        A[j][i] = std::pow(points[i], j);
      }
      else if (i == i_0) {
        A[j][i] = std::pow(0, j);
      }
      else {
        A[j][i] = std::pow(points[i - 1], j);
      }
    }
  }
  for (unsigned int j = 0; j < N + 2; j++) {
    A[L][N + 1] = fact(L);
  }

  // решение СЛАУ методом Гаусса
  for (unsigned int i = 0; i < N + 1; i++) { // выбор ведущей строки
    float A_ii = A[i][i];
    for (unsigned int j = 0; j < N + 2; j++) { // нормировка ведущей строки
      A[i][j] /= A_ii;
    }
    for (unsigned int n = 0; n < N + 1; n++) { // проход по строкам
      if (n != i) {
        int A_ni = A[n][i];
        for (unsigned int k = 0; k < N + 2; k++) { // проход по выбранной строке
          A[n][k] -= A_ni * A[i][k];
        }
      }
    }
  }

  // формирование ответа
  DerivativeCoef<RealType, N> derivativeCoef;
  for (unsigned int i = 0; i < N + 1; i++) {
    if (i < i_0) {
      derivativeCoef.otherCoefs[i] = A[i][N + 1];

    }
    else if (i == i_0) {
      derivativeCoef.centralCoef = A[i][N + 1];
    }
    else {
      derivativeCoef.otherCoefs[i - 1] = A[i][N + 1];
    }
  }
  return derivativeCoef;
}

int main() {
  const unsigned int N = 3, L = 2;
  std::array<float, N> points = {1, 2, 3};
  DerivativeCoef<float, N> derivativeCoef = calcDerivativeCoef<float, N, L>(points);
  std::cout << "centralCoef = " << derivativeCoef.centralCoef << std::endl;
  std::cout << "otherCoefs = {";
  for (unsigned int i = 0; i < N - 1; i++) {
      std::cout << derivativeCoef.otherCoefs[i] << ", ";
  }
  std::cout << derivativeCoef.otherCoefs[N - 1] << '}';
}
