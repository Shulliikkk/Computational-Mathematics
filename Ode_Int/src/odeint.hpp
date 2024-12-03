#ifndef ODEINT_HPP
#define ODEINT_HPP

#include <array>
#include <vector>
#include <functional>
#include <cmath>

struct  BT_DP45 {
  static constexpr std::size_t stages = 7;
  static constexpr std::size_t order = 5;
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

template<typename T>
struct Point {
  double time;
  T state;
};

template<std::size_t N>
double norm_error(const std::array<double, N>& next_state_1, const std::array<double, N>& next_state_2) {
  double diff, sum_sqr = 0;
  for (std::size_t i = 0; i < N; i++) {
    diff = next_state_1[i] - next_state_2[i];
    sum_sqr += diff * diff;
  }
  return std::sqrt(sum_sqr);
}

template<std::size_t N>
  std::vector<Point<std::array<double, N>>> odeint(const std::function<std::array<double, N>(double, std::array<double, N>)> rhs,
  const std::array<double, N>& init_state,
  const double end_time,
  const double init_h,
  const double tolerance
) {
  static constexpr std::size_t stages = BT_DP45::stages, oreder_p_1 = BT_DP45::order;;
  double curr_time = 0, curr_h = init_h, error;
  std::array<double, N> curr_state, next_state_1, next_state_2, next_arg_k_j;;
  std::array<std::array<double, N>, stages> k, sum_l2j;
  std::array<double, N> sum_j2s_b_1, sum_j2s_b_2;
  std::vector<Point<std::array<double, N>>> solution;
  for (std::size_t i = 0; i < N; i++) {
    curr_state[i] = init_state[i];
  }
  solution.push_back(Point<std::array<double, N>>{curr_time, curr_state});
  while (curr_time < end_time) {
    //while(true) {
      // calc k_1, ..., k_j, ..., k_s-1
      for (std::size_t j = 0; j < stages; j++) {
        for (std::size_t m = 0; m < N; m++) {
          sum_l2j[j][m] = 0;
          for (std::size_t l = 0; l < j; l++) {
            sum_l2j[j][m] += BT_DP45::tabel[j][l] * k[l][m];
          }
          next_arg_k_j[m] = curr_state[m] + curr_h  * sum_l2j[j][m];
        }
        for (std::size_t m = 0; m < N; m++) {
          k[j][m] = rhs(curr_time + BT_DP45::c_column[j] * curr_h,  next_arg_k_j)[m];
        }
      }
      // calc u_i+1^(1), u_i+1^(2)
      for (std::size_t m = 0; m < N; m++) {
        sum_j2s_b_1[m] = 0;
        sum_j2s_b_2[m] = 0;
        for (std::size_t j = 0; j < stages; j++) {
          sum_j2s_b_1[m] += BT_DP45::b_string_1[j] * k[j][m];
          sum_j2s_b_2[m] += BT_DP45::b_string_2[j] * k[j][m];
        }
        next_state_1[m] = curr_state[m] + curr_h * sum_j2s_b_1[m];
        next_state_2[m] = curr_state[m] + curr_h * sum_j2s_b_2[m];
      }
      error = norm_error(next_state_1, next_state_2);
      if (error > tolerance) {
          curr_h *=  0.9 * std::pow(tolerance / error, 1. / oreder_p_1);
      }
      else {
          curr_time += curr_h;
          curr_state = next_state_2;
          solution.push_back(Point<std::array<double, N>>{curr_time, curr_state});
          curr_h *= 0.9 * std::pow(tolerance / error, 1. / oreder_p_1);
          if (curr_time + curr_h > end_time) {
            curr_h = end_time - curr_time;
          }
      }
  }
  return solution;
}

#endif //ODEINT_HPP
