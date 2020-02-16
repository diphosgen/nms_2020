#ifndef VECTOR_FIXED_HXX_INCLUDED
#define VECTOR_FIXED_HXX_INCLUDED

#include "basic_libs.hxx"
#include "basic_math_functions.hxx"

namespace Math_structures
{

template<
    typename T,
    int capacity
    >
class Vector_fixed
{
public:

    static_assert(capacity > 0);

    using value_type = T;

    constexpr Vector_fixed() noexcept = default;

    constexpr Vector_fixed(const Vector_fixed<T, capacity>& constr_vector) noexcept
        :   vector_data{constr_vector.vector_data}
    { }

    constexpr Vector_fixed(Vector_fixed<T, capacity>&& constr_vector) noexcept
        :   vector_data{std::move(constr_vector.vector_data)}
    { }

    constexpr Vector_fixed(const int data_size) noexcept
    {
        assert(data_size == capacity);
    }

    constexpr Vector_fixed(const int data_size, const T& initial_value) noexcept
        :   vector_data{std::array<T, capacity>{}.fill(initial_value)}
    {
        assert(data_size == capacity);
    }

    constexpr Vector_fixed(const T (&data)[capacity]) noexcept
    {
        for (int i = 0; i < capacity; ++i) {
            this->vector_data[i] = data[i];
        }
    }

    constexpr Vector_fixed(std::initializer_list<T> data) noexcept
    {
        assert(data.size() == capacity);

        int n = 0;
        for (auto const& elem: data) {
            this->vector_data[n++] = elem;
        }
    }

    ~Vector_fixed() noexcept = default;

    constexpr Vector_fixed<T, capacity>& operator=(const Vector_fixed<T, capacity>& assign_vector) noexcept
    {
        if (&assign_vector != this) {
            this->vector_data = assign_vector.vector_data;
        } else {}

        return *this;
    }

    constexpr Vector_fixed<T, capacity>& operator=(Vector_fixed<T, capacity>&& assign_vector) noexcept
    {
        if (&assign_vector != this) {
            this->vector_data = std::move(assign_vector.vector_data);
        } else {}

        return *this;
    }

    constexpr T& operator[](const int n) noexcept
    {
        assert(n >= 0 && n < capacity);
        return this->vector_data[n];
    }

    constexpr const T& operator[](const int n) const noexcept
    {
        assert(n >= 0 && n < capacity);
        return this->vector_data[n];
    }

    constexpr int size() const noexcept
    {
        return capacity;
    }

    constexpr void set_null_value() noexcept
    {
        for (auto& elem: this->vector_data) {

            constexpr bool is_simple_standard_math_type =
                                        std::is_arithmetic_v<std::decay_t<T>> ||
                                        std::is_same_v<std::decay_t<T>, std::complex<float>> ||
                                        std::is_same_v<std::decay_t<T>, std::complex<double>> ||
                                        std::is_same_v<std::decay_t<T>, std::complex<long double>>;

            if constexpr (is_simple_standard_math_type) {
                elem = T(0);
            } else {
                elem.set_null_value();
            }
        }
    }

    constexpr Vector_fixed<T, capacity>& operator+=(const Vector_fixed<T, capacity>& add_vector) noexcept
    {
        for (int i = 0; i < capacity; ++i) {
            this->vector_data[i] += add_vector.vector_data[i];
        }

        return *this;
    }

    constexpr Vector_fixed<T, capacity>& operator-=(const Vector_fixed<T, capacity>& sub_vector) noexcept
    {
        for (int i = 0; i < capacity; ++i) {
            this->vector_data[i] -= sub_vector.vector_data[i];
        }

        return *this;
    }

    constexpr Vector_fixed<T, capacity>& operator*=(const T& multiplier) noexcept
    {
        for (int i = 0; i < capacity; ++i) {
            this->vector_data[i] *= multiplier;
        }

        return *this;
    }

    constexpr Vector_fixed<T, capacity>& operator/=(const T& denom) noexcept
    {
        for (int i = 0; i < capacity; ++i) {
            this->vector_data[i] /= denom;
        }

        return *this;
    }

private:

    std::array<T, capacity> vector_data{};

};

template<
    typename T_lhs,
    typename T_rhs,
    int capacity
    >
constexpr bool operator==(const Vector_fixed<T_lhs, capacity>& lhs,
                          const Vector_fixed<T_rhs, capacity>& rhs) noexcept
{
    auto comp = [](const T_lhs& lhs, const T_rhs& rhs) constexpr -> bool
    {
        constexpr bool is_floating_type =
                        std::is_floating_point_v<std::decay_t<T_lhs>> ||
                        std::is_floating_point_v<std::decay_t<T_rhs>>;

        if constexpr (is_floating_type) {
            return compare_floating_numbers(lhs, rhs);
        } else {
            return lhs == rhs;
        }
    };

    for (int i = 0; i < capacity; ++i) {
        if (!comp(lhs[i], rhs[i])) {
            return false;
        }
    }

    return true;
}

template<
    typename T_lhs,
    typename T_rhs,
    int capacity
    >
constexpr Vector_fixed<decltype(std::declval<T_lhs>() + std::declval<T_rhs>()), capacity>
    operator+(const Vector_fixed<T_lhs, capacity>& lhs,
              const Vector_fixed<T_rhs, capacity>& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() +
                              std::declval<T_rhs>());

    Vector_fixed<res_type, capacity> res{};

    for (int i = 0; i < capacity; ++i) {
        res[i] = lhs[i] + rhs[i];
    }

    return res;
}

template<
    typename T_lhs,
    typename T_rhs,
    int capacity
    >
constexpr Vector_fixed<decltype(std::declval<T_lhs>() - std::declval<T_rhs>()), capacity>
    operator-(const Vector_fixed<T_lhs, capacity>& lhs,
              const Vector_fixed<T_rhs, capacity>& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() -
                              std::declval<T_rhs>());

    Vector_fixed<res_type, capacity> res{};

    for (int i = 0; i < capacity; ++i) {
        res[i] = lhs[i] - rhs[i];
    }

    return res;
}

template<
    typename T_lhs,
    typename T_rhs,
    int capacity
    >
constexpr Vector_fixed<decltype(std::declval<T_lhs>() * std::declval<T_rhs>()), capacity>
    operator*(const Vector_fixed<T_lhs, capacity>& lhs, const T_rhs& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() *
                              std::declval<T_rhs>());

    Vector_fixed<res_type, capacity> res{};

    for (int i = 0; i < capacity; ++i) {
        res[i] = lhs[i] * rhs;
    }

    return res;
}

template<
    typename T_lhs,
    typename T_rhs,
    int capacity
    >
constexpr Vector_fixed<decltype(std::declval<T_lhs>() * std::declval<T_rhs>()), capacity>
    operator*(const T_lhs& lhs, const Vector_fixed<T_rhs, capacity>& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() *
                              std::declval<T_rhs>());

    Vector_fixed<res_type, capacity> res{};

    for (int i = 0; i < capacity; ++i) {
        res[i] = lhs * rhs[i];
    }

    return res;
}

template<
    typename T_lhs,
    typename T_rhs,
    int capacity
    >
constexpr decltype(std::declval<decltype(std::declval<T_lhs>() * std::declval<T_lhs>())>() +
                   std::declval<decltype(std::declval<T_rhs>() * std::declval<T_rhs>())>())
    operator*(const Vector_fixed<T_lhs, capacity>& lhs,
              const Vector_fixed<T_rhs, capacity>& rhs) noexcept
{
    using res_type = decltype(std::declval<
                                decltype(std::declval<T_lhs>() *
                                         std::declval<T_lhs>())
                                >() +
                              std::declval<
                                decltype(std::declval<T_rhs>() *
                                         std::declval<T_rhs>())
                                >());

    res_type res = res_type();

    for (int i = 0; i < capacity; ++i) {
        res += lhs[i] * rhs[i];
    }

    return res;
}

template<
    typename T_lhs,
    typename T_rhs,
    int capacity
    >
constexpr Vector_fixed<decltype(std::declval<T_lhs>() / std::declval<T_rhs>()), capacity>
    operator/(const Vector_fixed<T_lhs, capacity>& lhs, const T_rhs& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() /
                              std::declval<T_rhs>());

    Vector_fixed<res_type, capacity> res{};

    for (int i = 0; i < capacity; ++i) {
        res[i] = lhs[i] / rhs;
    }

    return res;
}

template<
    typename T,
    int capacity
    >
constexpr auto max_abs(const Vector_fixed<T, capacity>& v) noexcept
{
    constexpr bool is_simple_standard_math_type =
                                std::is_arithmetic_v<std::decay_t<T>> ||
                                std::is_same_v<std::decay_t<T>, std::complex<float>> ||
                                std::is_same_v<std::decay_t<T>, std::complex<double>> ||
                                std::is_same_v<std::decay_t<T>, std::complex<long double>>;

    if constexpr (is_simple_standard_math_type) {
        const int vector_size = v.size();
        using max_value_t = decltype(std::abs(std::declval<decltype(v[0])>()));
        max_value_t max_value = max_value_t(0);
        for (int i = 0; i < vector_size; ++i) {
            max_value = std::max(max_value, std::abs(v[i]));
        }
        return max_value;
    } else {
        const int vector_size = v.size();
        using max_value_t = decltype(max_abs(std::declval<decltype(v[0])>()));
        max_value_t max_value = max_value_t();
        for (int i = 0; i < vector_size; ++i) {
            max_value = std::max(max_value, max_abs(v[i]));
        }
        return max_value;
    }
}

}

#endif // VECTOR_FIXED_HXX_INCLUDED
