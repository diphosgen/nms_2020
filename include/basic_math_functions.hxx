#ifndef BASIC_MATH_FUNCTIONS_HXX_INCLUDED
#define BASIC_MATH_FUNCTIONS_HXX_INCLUDED

#include "basic_libs.hxx"

namespace Math_structures
{

template<
    typename T_lhs,
    typename T_rhs,
    typename T_abs = decltype(std::abs(std::declval<T_lhs>() -
                                       std::declval<T_rhs>())),
    typename = std::enable_if_t<
                    std::is_floating_point_v<std::decay_t<T_lhs>> ||
                    std::is_floating_point_v<std::decay_t<T_rhs>>
                    >
    >
constexpr bool compare_floating_numbers(const T_lhs& lhs, const T_rhs& rhs,
                                        const T_abs& eps_multplier = T_abs(100.0)) noexcept
{
    using common_type_res = std::common_type_t<T_lhs, T_rhs>;

    constexpr auto eps_base =
            std::max(static_cast<common_type_res>(std::numeric_limits<T_lhs>::epsilon()),
                     static_cast<common_type_res>(std::numeric_limits<T_rhs>::epsilon()));

    const T_abs eps = eps_multplier * static_cast<T_abs>(eps_base);

    const T_abs abs_diff = std::abs(lhs - rhs);

    if (abs_diff < eps) {

        return true;

    } else {

        const T_abs abs_sum = std::abs(lhs) + std::abs(rhs);
        const T_abs rel_diff = 2 * abs_diff / abs_sum;

        if (rel_diff < eps) {
            return true;
        } else {
            return false;
        }

    }
}

template<
    typename T,
    typename = std::enable_if_t<std::is_arithmetic_v<std::decay_t<T>>>
    >
constexpr auto max_abs(const T& value) noexcept
{
    return std::abs(value);
}

template<
    typename T,
    typename = std::enable_if_t<std::is_floating_point_v<std::decay_t<T>>>
    >
constexpr auto max_abs(const std::complex<T>& value) noexcept
{
    return std::abs(value);
}

}

#endif // BASIC_MATH_FUNCTIONS_HXX_INCLUDED
