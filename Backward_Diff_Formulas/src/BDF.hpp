#ifndef BDF_HPP
#define BDF_HPP

#include <Eigen/Dense>
#include <vector>

template<std::size_t N>
struct Point {
  double time;
  Eigen::Vector<double, N> state;
};

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

template<std::size_t N, std::size_t M>
Eigen::Vector<Eigen::Vector<double, N>, M> RK4(const std::function<Eigen::Vector<double, N>(double, Eigen::Vector<double, N>)> rhs,
  const Eigen::Vector<double, N> init_state,
  const double h
) {
  Eigen::Vector<double, 3> a(1. / 2., 1. / 6., 1. / 3.);
  Eigen::Vector<double, N> k_1, k_2, k_3, k_4;
  Eigen::Vector<double, N> curr_state = init_state;
  Eigen::Vector<Eigen::Vector<double, N>, M> solution;
  double curr_time = 0;
  for (std::size_t j = 0; j < M; j++) {
    solution[j] = curr_state;
    k_1 = rhs(curr_time, curr_state);
    k_2 = rhs(curr_time + a[0] * h, curr_state + a[0] * h * k_1);
    k_3 = rhs(curr_time + a[0] * h, curr_state + a[0] * h * k_2);
    k_4 = rhs(curr_time + h, curr_state +  h * k_3);
    curr_state += h * (a[1] * k_1 + a[2] * k_2 + a[2] * k_3 + a[1] * k_4);
    curr_time += h;
  }
  return solution;
}

template<std::size_t N>
std::vector<Point<N>> BDF(const std::function<Eigen::Vector<double, N>(double, Eigen::Vector<double, N>)> rhs,
                               const Eigen::Vector<double, N> init_state,
                               const double end_time,
                               const double h,
                               const double tolerance,
                               const std::size_t max_iter
) {
  Eigen::Vector<double, 7> a(-77. / 10., 15., -10., 5., -3. / 2., 1. / 5., 6.);
  Eigen::Vector<double, 7> b(360. / 147., -450. / 147., 400. / 147., -225. / 147., 72. / 147., -10. / 147., 60. / 147.);
  Eigen::Vector<double, N> next_state_pridict, next_state_correct, add_coef;
  Eigen::Vector<Eigen::Vector<double, N>, 6> prev_states = RK4<N, 6>(rhs, init_state, h);
  std::vector<Point<N>> solution;
  const double mult_coef =  b[6] * h;
  for (std::size_t i = 0; i < 6; i++) {
    solution.push_back(Point<N>{i * h, prev_states[i]});
  }
  for (double curr_time = 5 * h; curr_time <= end_time; curr_time += h) {
    add_coef = b[0] * prev_states[5] + b[1] * prev_states[4] + b[2] * prev_states[3] +
               b[3] * prev_states[2] + b[4] * prev_states[1] + b[5] * prev_states[0];
    next_state_pridict = a[0] * prev_states[5] + a[1] * prev_states[4] + a[2] * prev_states[3] +
                         a[3] * prev_states[2] + a[4] * prev_states[1] + a[5] * prev_states[0] +
                         a[6] * h * rhs(curr_time, prev_states[5.]);
    next_state_correct = FPI<N>(rhs, curr_time, next_state_pridict,
                                tolerance, max_iter, add_coef, mult_coef);
    prev_states[0] = prev_states[1];
    prev_states[1] = prev_states[2];
    prev_states[2] = prev_states[3];
    prev_states[3] = prev_states[4];
    prev_states[4] = prev_states[5];
    prev_states[5] = next_state_correct;
    solution.push_back(Point<N>{curr_time, prev_states[5]});
  }
  return solution;
}

#endif //BDF_HPP
