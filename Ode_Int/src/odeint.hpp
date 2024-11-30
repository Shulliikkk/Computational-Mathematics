#ifndef ODEINT_HPP
#define ODEINT_HPP

#include <array>
#include <vector>
#include <functional>
#include <cmath>

struct DP45 {
  static constexpr std::size_t stages = 7;
  static constexpr std::array<std::array<double, stages>, stages> tabel{{{0, 0, 0, 0, 0, 0, 0},
                                                                           {1. / 5., 0, 0, 0, 0, 0, 0},
                                                                           {3. / 40., 9. / 40., 0, 0, 0, 0, 0},
                                                                           {44. / 45., -56. / 15., 32. / 9., 0, 0, 0, 0},
                                                                           {19372. / 6561., -25360. / 2187., 64448. / 6561., -212. / 729., 0, 0, 0},
                                                                           {9017. / 3168., -355. / 33., 46732. / 5247., 49. / 176., -5103. / 18656., 0, 0},
                                                                           {35. / 384., 0, 500. / 1113., 125. / 192., -2187. / 6784., 11. / 84., 0}}};
  static constexpr std::array<double, stages> c_column{0, 1. / 5., 3. / 10., 4. / 5., 8. / 9., 1., 1.};
  static constexpr std::array<double, stages> b_string_1{35. / 384., 0, 500. / 1113., 125. / 192., -2187. / 6784., 11. / 84., 0};
  static constexpr std::array<double, stages> b_string_2{5179. / 57600., 0, 7571. / 16695., 393. / 640., -92097. / 339200., 187. / 2100., 1. / 40.};
};

struct RKF45 {
  static constexpr std::size_t stages = 6;
  static constexpr std::array<std::array<double, stages>, stages> tabel{{{0, 0, 0, 0, 0, 0},
                                                                           {1. / 4., 0, 0, 0, 0, 0},
                                                                           {3. / 32., 9. / 32., 0, 0, 0, 0},
                                                                           {1932. / 2197., -7200. / 2197., 7296. / 2197., 0, 0, 0},
                                                                           {439. / 216., -8., 3680. / 513., -845. / 4104., 0, 0},
                                                                           {-8. / 27., -2., 3544. / 2565., 1859. / 4104., -11. / 40., 0}}};
  static constexpr std::array<double, stages> c_column{0, 1. / 4., 3. / 8., 12. / 13., 1., 1. / 2.};
  static constexpr std::array<double, stages> b_string_1{25. / 216., 0, 1408. / 2565., 2197. / 4104., -1. / 5., 0.};
  static constexpr std::array<double, stages> b_string_2{16. / 135., 0, 6656. / 12825., 28561. / 56430., -9. / 50., 2. / 55.};
};

struct BS23 {
  static constexpr std::size_t stages = 4;
  static constexpr std::array<std::array<double, stages>, stages> tabel{{{0, 0, 0, 0},
                                                                           {1. / 2., 0, 0},
                                                                           {0., 3. / 4., 0, 0},
                                                                           {2. / 9., 1. / 3., 4. / 9., 0}}};
  static constexpr std::array<double, stages> c_column{0, 1. / 2, 3. / 4., 1.};
  static constexpr std::array<double, stages> b_string_1{2. / 9., 1. / 3., 4. / 9., 0.};
  static constexpr std::array<double, stages> b_string_2{7. / 24., 1. / 4., 1. / 3., 1. / 8.};
};

template<std::size_t N>
double norm_error(std::array<double, N>& next_state_1, std::array<double, N>& next_state_2) {
  double diff, sum_sqr = 0;
  for (std::size_t i = 0; i < N; i++) {
    diff = std::abs(next_state_1[i] - next_state_2[i]);
    sum_sqr += diff * diff;
  }
  return std::sqrt(sum_sqr);
}

template<class Butcher_Table, std::size_t N>
std::vector<std::array<double, N>> odeint(const std::function<std::array<double, N>(double, std::array<double, N>)> rhs,
  std::array<double, N>& init_state,
  const double end_time,
  const double init_h,
  double tolerance
) {
  static constexpr Butcher_Table BT;
  static constexpr std::size_t stages = BT.stages;
  double curr_time = 0, curr_h = init_h;
  std::array<double, N> curr_state = init_state, next_state_1, next_state_2, next_arg_k_j;;
  std::array<std::array<double, N>, stages> k, sum_l2j;
  std::array<double, N> sum_j2s_b_1,  sum_j2s_b_2;
  std::vector<std::array<double, N>> solution;
  double error;
  while (curr_time <= end_time) {
    solution.push_back(curr_state);
    while(true) {
      // calc k_1, ..., k_j, ..., k_s-1
      for (std::size_t j = 0; j < stages; j++) {
        for (std::size_t m = 0; m < N; m++) {
          sum_l2j[j][m] = 0;
          for (std::size_t l = 0; l < j; l++) {
            sum_l2j[j][m] += BT.tabel[j][l] * k[l][m];
          }
          next_arg_k_j[m] = curr_state[m] + curr_h  * sum_l2j[j][m];
        }
        for (std::size_t m = 0; m < N; m++) {
          k[j][m] = rhs(curr_time + BT.c_column[j] * curr_h,  next_arg_k_j)[m];
        }
      }
      // calc u_i+1^(1), u_i+1^(2)
      for (std::size_t m = 0; m < N; m++) {
        sum_j2s_b_1[m] = 0;
        sum_j2s_b_2[m] = 0;
        for (std::size_t j = 0; j < stages; j++) {
          sum_j2s_b_1[m] += BT.b_string_1[j] * k[j][m];
          sum_j2s_b_2[m] += BT.b_string_2[j] * k[j][m];
        }
        next_state_1[m] = curr_state[m] + curr_h * sum_j2s_b_1[m];
        next_state_2[m] = curr_state[m] + curr_h * sum_j2s_b_2[m];
      }
      error = norm_error(next_state_1, next_state_2);
      curr_h *=  0.9 * std::pow(tolerance / error, 1. / 5.);
      if (error <= tolerance) {
        break;
      }
    }

    if (curr_time + curr_h >= end_time) {
      curr_h = end_time - curr_time;
      if (curr_h - 0.0 < 1e-16) {
        curr_state = next_state_2;
        solution.push_back(curr_state);
        break;
      }
    }
    curr_time += curr_h;
    curr_state = next_state_2;
  }
  return solution;
}

#endif //ODEINT_HPP
