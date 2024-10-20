#include <iostream>
#include <functional>
#include <cmath>

template<typename T>
T f(T x) {
    return sin(x);
}

template<typename T, std::size_t N>
T elem_gauss_quad(const std::function<T(T)> func, const T x_prev, const T x_next) {
    constexpr std::array<T, N> nodes = []{
        if constexpr (N == 5) {
            return std::array<T, 5>{-0.906113631, -0.53846931, 0, 0.53846931, 0.906113631};
        }
        else if constexpr (N == 8) {
            return std::array<T, 8>{-0.96028986, -0.79666648, -0.52553242, -0.18343464, 0.18343464, 0.52553242, 0.79666648, 0.96028986};
        }
        else {
            static_assert(N == 5 || N == 8, "N must be either 5 or 8");
        }
     }();

    constexpr std::array<T, N> weights = []{
        if constexpr (N == 5) {
            return std::array<T, 5>{0.23692688, 0.47862868, 0.56888889, 0.47862868, 0.23692688};
        } else if constexpr (N == 8) {
            return std::array<T, 8>{0.10122854, 0.22238104, 0.31370664, 0.36268379, 0.36268379, 0.31370664, 0.22238104, 0.10122854};
        } else {
            static_assert(N == 5 || N == 8, "N must be either 5 or 8");
        }
    }();

    const T half_h = (x_next - x_prev) / 2.;
    const T centre_h = (x_next + x_prev) / 2.;
    T I = 0;
    for (std::size_t i = 0; i < N; i++) {
        I += weights[i] * func(centre_h + half_h * nodes[i]);
    }
    return half_h * I;
}

template<typename T, std::size_t N>
T gauss_quad(const std::function<T(T)> func, const T start, const T end, const std::size_t M) {
    const T h = (end - start) / M;
    T x_prev = start;
    T x_next = start + h;
    T I = 0;
    for (std::size_t i = 0; i < M; i++) {
        I += elem_gauss_quad<T, N>(func, x_prev, x_next);
        x_prev += h;
        x_next += h;
    }
    return I;
}

template<std::size_t N>
constexpr std::size_t deg_2() {
    std::size_t res = 1;
    for (std::size_t i = 0; i < N; i++) {
        res *= 2;
    }
    return res;
}

template<typename T, std::size_t N>
T runge_gauss_quad(const std::function<T(T)> func, const T start, const T end, const T err) {
    std::size_t M = 10;
    T I_h, I_h2, Delta_I, I;
    I_h = gauss_quad<T, N>(func, start, end, M);
    while (true) {
        I_h2 = gauss_quad<T, N>(func, start, end, 2 * M);
        Delta_I = (I_h - I_h2) / (deg_2<N + 1>() - 1);
        if (Delta_I < err || M > 100000) {
            I = I_h2 - Delta_I;
            break;
        }
        M *= 2;
        I_h = I_h2;
    }
    return I;
}

int main() {
    double err, eps = 1.0;
    double I, I_accur = 1.839071529076452452258864;
    for (std::size_t i = 0; i < 15; i++) {
        I = runge_gauss_quad<double, 5>(f<double>, 0, 10, eps);
        err = std::abs(I - I_accur);
        std::cout << err << ' ';
        eps *= 0.1;
    }
}
