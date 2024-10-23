#ifndef HERMIT_HPP
#define HERMIT_HPP

#include <array>

template<typename T_x, typename T_y, std::size_t N>
class Hermit_Interpolator {
  private:
    std::array<T_x, 2 * N> z;
    std::array<T_y, 2 * N> dvd_dif;
    std::pair<std::array<T_x, 2 * N>,std::array<T_y, 2 * N>> stretch(const std::array<T_x, N>& x, const std::array<T_y, N>& f) noexcept;

  public:
    Hermit_Interpolator(const std::array<T_x, N> &x, const std::array<T_y, N>& f, const std::array<T_y, N>& df) noexcept;
    T_y interpolate(const T_x& x_0) const noexcept;
};

#include "Hermit.cpp"

#endif //HERMIT_HPP
