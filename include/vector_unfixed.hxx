#ifndef VECTOR_UNFIXED_HXX_INCLUDED
#define VECTOR_UNFIXED_HXX_INCLUDED

#include "basic_libs.hxx"
#include "basic_math_functions.hxx"

namespace Math_structures
{

template<typename T = double>
class Vector_unfixed
{
public:

    using value_type = T;

    constexpr Vector_unfixed() noexcept = default;

    constexpr Vector_unfixed(const Vector_unfixed<T>& constr_vector) noexcept
        :   vector_size{constr_vector.vector_size},
            vector_data{constr_vector.vector_data}
    { }

    constexpr Vector_unfixed(Vector_unfixed<T>&& constr_vector) noexcept
        :   vector_size(std::move(constr_vector.vector_size)),
            vector_data{std::move(constr_vector.vector_data)}
    { }

    constexpr Vector_unfixed(const int data_size) noexcept
        :   vector_size{data_size},
            vector_data(data_size)
    { }

    constexpr Vector_unfixed(const int data_size, const T& initial_value) noexcept
        :   vector_size{data_size},
            vector_data(data_size, initial_value)
    { }

    template<int data_size>
    constexpr Vector_unfixed(const T (&data)[data_size]) noexcept
        :   vector_size(data_size),
            vector_data(data_size)
    {
        for (int i = 0; i < data_size; ++i) {
            this->vector_data[i] = data[i];
        }
    }

    constexpr Vector_unfixed(std::initializer_list<T> data) noexcept
        :   vector_size(data.size()),
            vector_data(data)
    { }

    ~Vector_unfixed() noexcept = default;

    constexpr Vector_unfixed<T>& operator=(const Vector_unfixed<T>& assign_vector) noexcept
    {
        if (&assign_vector != this) {
            if (this->vector_size == 0) {

                this->vector_size = assign_vector.vector_size;
                this->vector_data = assign_vector.vector_data;

            } else {

                assert(this->vector_size == assign_vector.vector_size);
                this->vector_data = assign_vector.vector_data;

            }
        } else {}

        return *this;
    }

    constexpr Vector_unfixed<T>& operator=(Vector_unfixed<T>&& assign_vector) noexcept
    {
        if (&assign_vector != this) {
            if (this->vector_size == 0) {

                this->vector_size = std::move(assign_vector.vector_size);
                this->vector_data = std::move(assign_vector.vector_data);

            } else {

                assert(this->vector_size == assign_vector.vector_size);
                this->vector_data = std::move(assign_vector.vector_data);

            }
        } else {}

        return *this;
    }

    constexpr T& operator[](const int n) noexcept
    {
        assert(n >= 0 && n < this->vector_size);
        return this->vector_data[n];
    }

    constexpr const T& operator[](const int n) const noexcept
    {
        assert(n >= 0 && n < this->vector_size);
        return this->vector_data[n];
    }

    constexpr int size() const noexcept
    {
        return this->vector_size;
    }

    constexpr void set_null_value() noexcept
    {
        if (this->vector_size != 0) {
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
        } else {}
    }

    constexpr Vector_unfixed<T>& operator+=(const Vector_unfixed<T>& add_vector) noexcept
    {
        if (this->vector_size == 0) {

            this->vector_size = add_vector.vector_size;
            this->vector_data = add_vector.vector_data;

        } else {
            assert(this->vector_size == add_vector.vector_size);

            for (int i = 0; i < this->vector_size; ++i) {
                this->vector_data[i] += add_vector.vector_data[i];
            }
        }

        return *this;
    }

    constexpr Vector_unfixed<T>& operator-=(const Vector_unfixed<T>& sub_vector) noexcept
    {
        if (this->vector_size == 0) {
            this->vector_size = sub_vector.vector_size;
            this->set_null_value();
        } else {
            assert(this->vector_size == sub_vector.vector_size);
        }

        for (int i = 0; i < this->vector_size; ++i) {
            this->vector_data[i] -= sub_vector.vector_data[i];
        }

        return *this;
    }

    constexpr Vector_unfixed<T>& operator*=(const T& multiplier) noexcept
    {
        for (int i = 0; i < this->vector_size; ++i) {
            this->vector_data[i] *= multiplier;
        }

        return *this;
    }

    constexpr Vector_unfixed<T>& operator/=(const T& denom) noexcept
    {
        for (int i = 0; i < this->vector_size; ++i) {
            this->vector_data[i] /= denom;
        }

        return *this;
    }

private:

    int vector_size = 0;
    std::vector<T> vector_data{};

};

template<typename T_lhs, typename T_rhs>
constexpr bool operator==(const Vector_unfixed<T_lhs>& lhs,
                          const Vector_unfixed<T_rhs>& rhs) noexcept
{
    assert(lhs.size() == rhs.size());

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

    const int vector_size = lhs.size();

    for (int i = 0; i < vector_size; ++i) {
        if (!comp(lhs[i], rhs[i])) {
            return false;
        }
    }

    return true;
}

template<typename T_lhs, typename T_rhs>
constexpr Vector_unfixed<decltype(std::declval<T_lhs>() +
                                  std::declval<T_rhs>())>
    operator+(const Vector_unfixed<T_lhs>& lhs,
              const Vector_unfixed<T_rhs>& rhs) noexcept
{
    if (rhs.size() == 0) {

        assert(lhs.size() != 0);
        return lhs;

    } else {

        assert(lhs.size() == rhs.size());

        const int vector_size = lhs.size();

        using res_type = decltype(std::declval<T_lhs>() +
                                  std::declval<T_rhs>());

        Vector_unfixed<res_type> res(vector_size);

        for (int i = 0; i < vector_size; ++i) {
            res[i] = lhs[i] + rhs[i];
        }

        return res;

    }
}

template<typename T_lhs, typename T_rhs>
constexpr Vector_unfixed<decltype(std::declval<T_lhs>() -
                                  std::declval<T_rhs>())>
    operator-(const Vector_unfixed<T_lhs>& lhs,
              const Vector_unfixed<T_rhs>& rhs) noexcept
{
    if (rhs.size() == 0) {

        assert(lhs.size() != 0);
        return lhs;

    } else {

        assert(lhs.size() == rhs.size());

        const int vector_size = lhs.size();

        using res_type = decltype(std::declval<T_lhs>() -
                                  std::declval<T_rhs>());

        Vector_unfixed<res_type> res(vector_size);

        for (int i = 0; i < vector_size; ++i) {
            res[i] = lhs[i] - rhs[i];
        }

        return res;

    }
}

template<typename T_lhs, typename T_rhs>
constexpr Vector_unfixed<decltype(std::declval<T_lhs>() *
                                  std::declval<T_rhs>())>
    operator*(const Vector_unfixed<T_lhs>& lhs, const T_rhs& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() *
                              std::declval<T_rhs>());

    const int vector_size = lhs.size();

    Vector_unfixed<res_type> res(vector_size);
    for (int i = 0; i < vector_size; ++i) {
        res[i] = lhs[i] * rhs;
    }

    return res;
}

template<typename T_lhs, typename T_rhs>
constexpr Vector_unfixed<decltype(std::declval<T_lhs>() *
                                  std::declval<T_rhs>())>
    operator*(const T_lhs& lhs, const Vector_unfixed<T_rhs>& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() *
                              std::declval<T_rhs>());

    const int vector_size = rhs.size();

    Vector_unfixed<res_type> res(vector_size);
    for (int i = 0; i < vector_size; ++i) {
        res[i] = lhs * rhs[i];
    }

    return res;
}

template<typename T_lhs, typename T_rhs>
constexpr decltype(std::declval<decltype(std::declval<T_lhs>() * std::declval<T_lhs>())>() +
                   std::declval<decltype(std::declval<T_rhs>() * std::declval<T_rhs>())>())
    operator*(const Vector_unfixed<T_lhs>& lhs,
              const Vector_unfixed<T_rhs>& rhs) noexcept
{
    assert(lhs.size() == rhs.size());

    const int vector_size = lhs.size();

    using res_type = decltype(std::declval<
                                decltype(std::declval<T_lhs>() *
                                         std::declval<T_lhs>())
                                >() +
                              std::declval<
                                decltype(std::declval<T_rhs>() *
                                         std::declval<T_rhs>())
                                >());

    res_type res = res_type();

    for (int i = 0; i < vector_size; ++i) {
        res += lhs[i] * rhs[i];
    }

    return res;
}

template<typename T_lhs, typename T_rhs>
constexpr Vector_unfixed<decltype(std::declval<T_lhs>() /
                                  std::declval<T_rhs>())>
    operator/(const Vector_unfixed<T_lhs>& lhs, const T_rhs& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() /
                              std::declval<T_rhs>());

    const int vector_size = lhs.size();

    Vector_unfixed<res_type> res(vector_size);
    for (int i = 0; i < vector_size; ++i) {
        res[i] = lhs[i] / rhs;
    }

    return res;
}

template<typename T>
constexpr auto max_abs(const Vector_unfixed<T>& v) noexcept
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

#endif // VECTOR_UNFIXED_HXX_INCLUDED
