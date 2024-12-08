
#ifndef BDF_HPP
#define BDF_HPP

#include <Eigen/Dense>
#include <vector>

template<std::size_t N>
struct Point {
  double time;
  Eigen::Vector<double, N> state;
};

struct  BT_RK4{
  static constexpr std::size_t stages = 4;
  static constexpr std::array<std::array<double, stages>, stages> tabel = {{{0, 0, 0, 0}, {1. / 2., 0, 0, 0}, {0, 1. / 2., 0, 0}, {0, 0, 1, 0}}};
  static constexpr std::array<double, stages> c_column = {0, 1. / 2., 1. / 2., 1};
  static constexpr std::array<double, stages> b_string = {1. / 6., 1. / 3., 1. / 3., 1. / 6.};
};

struct CF_BDF4 {
  static constexpr std::size_t order = 4;
  static constexpr std::array<double, order + 1> a{10. / 3., 6., -2., 1. / 3., 4.};
  static constexpr std::array<double, order + 1> b{48. / 25., -36. / 25., 16. / 25., -3. / 25., 12. / 25.};
};

struct  BT_RK5 {
  static constexpr std::size_t stages = 6;
  static constexpr std::array<std::array<double, stages>, stages> tabel = {{{0, 0, 0, 0, 0, 0},
                                                                            {1. / 3., 0, 0, 0, 0, 0},
                                                                            {4. / 25., 6. / 25., 0, 0, 0, 0},
                                                                            {1. / 4., -3., 15. / 4., 0, 0, 0},
                                                                            {2. / 27., 10. / 9., -50. / 81., 8. / 81., 0, 0},
                                                                            {2. / 25., 12. / 25., 2. / 15., 8. / 75., 0., 0}}};
  static constexpr std::array<double, stages> c_column = {0, 1. / 3., 2. / 5., 1., 2. / 3., 4. / 5.};
  static constexpr std::array<double, stages> b_string = {23. / 192., 0, 125. / 192., 0., -27. / 64., 125. / 192.};
};

struct  BT_RK6 {
  static constexpr std::size_t stages = 7;
  static constexpr std::array<std::array<double, stages>, stages> tabel = {{{0, 0, 0, 0, 0, 0, 0},
                                                                            {2. / 5., 0, 0, 0, 0, 0, 0},
                                                                            {0., 4. / 5., 0, 0, 0, 0, 0},
                                                                            {169. / 1458., 110. / 729., -65. / 1458., 0, 0, 0, 0},
                                                                            {-44. / 675., -88. / 135., 76. / 351., 336 / 325., 0, 0, 0},
                                                                            {21. / 106., 0., -105. / 689., -324. / 689., 45. / 106., 0, 0},
                                                                            {-2517. / 4864., -55. / 38., 10615. / 31616., 567. / 7904., 7245. / 4864., 2597. / 2432., 0}}};
  static constexpr std::array<double, stages> c_column = {0, 2. / 5., 4. / 5., 2. / 9., 8. / 15., 0., 1.};
  static constexpr std::array<double, stages> b_string = {0, 0, 1375. / 4992., 6561. / 20384., 3375. / 12544., 53. / 768., 19. / 294.};
};

struct CF_BDF5 {
  static constexpr std::size_t order = 5;
  static constexpr std::array<double, order + 1> a{-65. / 12., 10., -5., 5. / 3., -1. / 4., 5.};
  static constexpr std::array<double, order + 1> b{300. / 137., -300. / 137., 200. / 137., -75. / 137., 12. / 137., 60. / 137.};
};


struct CF_BDF6 {
  static constexpr std::size_t order = 6;
  static constexpr std::array<double, order + 1> a{-77. / 10., 15., -10., 5., -3. / 2., 1. / 5., 6.};
  static constexpr std::array<double, order + 1> b{60. / 147., -450. / 147., 400. / 147., -225. / 147., 72. / 147., -10. / 147., 60. / 147.};
};

/*
template<std::size_t N>
Eigen::Vector<double, N> Newton(const std::function<Eigen::Vector<double, N>(const double, const Eigen::Vector<double, N>)> F,
                                const std::function<Eigen::Matrix<double, N, N>(const double, const Eigen::Vector<double, N>)> J,
                                const double curr_time,
                                const Eigen::Vector<double, N> init_approx,
                                const double tolerance,
                                const std::size_t max_iter,
                                const Eigen::Vector<double, N> add_coef,
                                const double mult_coef
) {
  Eigen::Vector<double, N> curr_state = init_approx;
  for (std::size_t i = 0; i < max_iter; i++) {
    Eigen::Vector<double, N> F_curr = add_coef + mult_coef * F(curr_time, curr_state) - curr_state;
    Eigen::Matrix<double, N, N> J_curr = mult_coef * J(curr_time, curr_state) - Eigen::Matrix<double, N, N>::Identity();
    Eigen::Vector<double, N> delta_curr_state = -J_curr.lu().solve(F_curr);
    curr_state += delta_curr_state;
    if (delta_curr_state.norm() < tolerance) {
        return curr_state;
    }
  }
  return curr_state;
}
*/

template<std::size_t N>
Eigen::Vector<double, N> FPI(const std::function<Eigen::Vector<double, N>(double, Eigen::Vector<double, N>)> rhs,
                             const double curr_time,
                             const Eigen::Vector<double, N> init_approx,
                             const double tolerance,
                             const std::size_t max_iter,
                             const Eigen::Vector<double, N> add_coef,
                             const double mult_coef
) {
  Eigen::Vector<double, N> curr_state = init_approx, next_state;
  for (std::size_t i = 0; i < max_iter; i++) {
    next_state = add_coef + mult_coef * rhs(curr_time, curr_state);
    if ((curr_state - next_state).norm() < tolerance) {
      return curr_state;
    }
    curr_state = next_state;
  }
  return curr_state;
}

template<std::size_t N, class BT>
std::vector<Eigen::Vector<double, N>> RK(const std::function<Eigen::Vector<double, N>(double, Eigen::Vector<double, N>)> rhs,
                                         const Eigen::Vector<double, N> init_state,
                                         const double end_time,
                                         const double h
) {
  static constexpr std::size_t stages = BT::stages;
  Eigen::Vector<double, N> curr_state = init_state, next_state;
  Eigen::Vector<Eigen::Vector<double, N>, stages> k;
  Eigen::Vector<double, N> sum_j2s;
  Eigen::Vector<Eigen::Vector<double, N>, stages> sum_l2j;
  std::vector<Eigen::Vector<double, N>> solution;
  for (double curr_time = 0; curr_time <= end_time; curr_time += h) {
    solution.push_back(curr_state);
      // calc k_1, ..., k_j, ..., k_s-1
    for (std::size_t j = 0; j < stages; j++) {
      sum_l2j[j].setZero();
      for (std::size_t l = 0; l < j; l++) {
        sum_l2j[j] += BT::tabel[j][l] * k[l];
      }
      k[j] = rhs(curr_time + BT::c_column[j] * h, curr_state + h  * sum_l2j[j]);
    }
    // calc u_i+1^(1), u_i+1^(2)
    sum_j2s.setZero();
    for (std::size_t j = 0; j < stages; j++) {
      sum_j2s += BT::b_string[j] * k[j];
    }
    next_state = curr_state + h * sum_j2s;
    curr_state = next_state;
  }
  return solution;
}

template<std::size_t N, class BT, class CF>
std::vector<Point<N>> BDF(const std::function<Eigen::Vector<double, N>(double, Eigen::Vector<double, N>)> F,
                          //const std::function<Eigen::Matrix<double, N, N>(const double, const Eigen::Vector<double, N>)> J,
                          const Eigen::Vector<double, N> init_state,
                          const double end_time,
                          const double h,
                          const double tolerance = 1e-15,
                          const std::size_t max_iter = 100
) {
  static constexpr std::size_t order = CF::order;
  static constexpr std::array<double, order + 1> a = CF::a;
  static constexpr std::array<double, order + 1> b = CF::b;
  Eigen::Vector<double, N> next_state_pridict, next_state_correct, add_coef;
  std::vector<Eigen::Vector<double, N>> prev_states = RK<N, BT>(F, init_state, (order - 1) * h, h);
  std::vector<Point<N>> solution;
  const double mult_coef =  b[order] * h;
  for (std::size_t i = 0; i < order; i++) {
    solution.push_back(Point<N>{i * h, prev_states[i]});
  }
  for (double curr_time = (order - 1) * h; curr_time < end_time; curr_time += h) {
    add_coef.setZero();
    for (std::size_t j = 0; j < order; j++) {
      add_coef += b[j] * prev_states[order - 1 - j];
    }
    next_state_pridict.setZero();
    for (std::size_t j = 0; j < order; j++) {
      next_state_pridict += a[j] * prev_states[order - 1 - j];
    }
    next_state_pridict +=  a[order] * h * F(curr_time, prev_states[order - 1]);
    //next_state_correct = Newton<N>(F, J, curr_time + h, next_state_pridict,
                                //tolerance, max_iter, add_coef, mult_coef);
    next_state_correct = FPI<N>(F, curr_time + h, next_state_pridict,
                                tolerance, max_iter, add_coef, mult_coef);
    for (std::size_t j = 0; j < order - 1; j++) {
      prev_states[j] = prev_states[j + 1];
    }
    prev_states[order - 1] = next_state_correct;
    solution.push_back(Point<N>{curr_time + h, next_state_correct});
  }
  return solution;
}

#endif //BDF_HPP
