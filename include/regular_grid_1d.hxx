#ifndef REGULAR_GRID_1D_HXX_INCLUDED
#define REGULAR_GRID_1D_HXX_INCLUDED

#include "basic_libs.hxx"
#include "basic_math_functions.hxx"

namespace Math_structures
{

template<typename T = double>
class Regular_grid_1d
{
public:

    using value_type = T;

    constexpr Regular_grid_1d() noexcept = default;

    constexpr Regular_grid_1d(const Regular_grid_1d<T>& constr_grid) noexcept
        :   grid_size{constr_grid.grid_size},
            grid_data{constr_grid.grid_data}
    { }

    constexpr Regular_grid_1d(Regular_grid_1d<T>&& constr_grid) noexcept
        :   grid_size(std::move(constr_grid.grid_size)),
            grid_data{std::move(constr_grid.grid_data)}
    { }

    constexpr Regular_grid_1d(const int data_size) noexcept
        :   grid_size{data_size},
            grid_data(data_size)
    { }

    constexpr Regular_grid_1d(const int data_size, const T& initial_value) noexcept
        :   grid_size{data_size},
            grid_data(data_size, initial_value)
    { }

    template<int data_size>
    constexpr Regular_grid_1d(const T (&data)[data_size]) noexcept
        :   grid_size(data_size),
            grid_data(data_size)
    {
        for (int i = 0; i < data_size; ++i) {
            this->grid_data[i] = data[i];
        }
    }

    constexpr Regular_grid_1d(std::initializer_list<T> data) noexcept
        :   grid_size(data.size()),
            grid_data(data)
    { }

    ~Regular_grid_1d() noexcept = default;

    constexpr Regular_grid_1d<T>& operator=(const Regular_grid_1d<T>& assign_grid) noexcept
    {
        if (&assign_grid != this) {
            if (this->grid_size == 0) {

                this->grid_size = assign_grid.grid_size;
                this->grid_data = assign_grid.grid_data;

            } else {

                assert(this->grid_size == assign_grid.grid_size);
                this->grid_data = assign_grid.grid_data;

            }
        } else {}

        return *this;
    }

    constexpr Regular_grid_1d<T>& operator=(Regular_grid_1d<T>&& assign_grid) noexcept
    {
        if (&assign_grid != this) {
            if (this->grid_size == 0) {

                this->grid_size = std::move(assign_grid.grid_size);
                this->grid_data = std::move(assign_grid.grid_data);

            } else {

                assert(this->grid_size == assign_grid.grid_size);
                this->grid_data = std::move(assign_grid.grid_data);

            }
        } else {}

        return *this;
    }

    constexpr T& operator[](const int n) noexcept
    {
        assert(n >= 0 && n < this->grid_size);
        return this->grid_data[n];
    }

    constexpr const T& operator[](const int n) const noexcept
    {
        assert(n >= 0 && n < this->grid_size);
        return this->grid_data[n];
    }

    constexpr int size() const noexcept
    {
        return this->grid_size;
    }

    constexpr void set_null_value() noexcept
    {
        if (this->grid_size != 0) {
            for (auto& elem: this->grid_data) {

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

    constexpr Regular_grid_1d<T>& operator+=(const Regular_grid_1d<T>& add_grid) noexcept
    {
        if (this->grid_size == 0) {

            this->grid_size = add_grid.grid_size;
            this->grid_data = add_grid.grid_data;

        } else {
            assert(this->grid_size == add_grid.grid_size);

            for (int i = 0; i < this->grid_size; ++i) {
                this->grid_data[i] += add_grid.grid_data[i];
            }
        }

        return *this;
    }

    constexpr Regular_grid_1d<T>& operator-=(const Regular_grid_1d<T>& sub_grid) noexcept
    {
        if (this->grid_size == 0) {
            this->grid_size = sub_grid.grid_size;
            this->set_null_value();
        } else {
            assert(this->grid_size == sub_grid.grid_size);
        }

        for (int i = 0; i < this->grid_size; ++i) {
            this->grid_data[i] -= sub_grid.grid_data[i];
        }

        return *this;
    }

    constexpr Regular_grid_1d<T>& operator*=(const T& multiplier) noexcept
    {
        for (int i = 0; i < this->grid_size; ++i) {
            this->grid_data[i] *= multiplier;
        }

        return *this;
    }

    constexpr Regular_grid_1d<T>& operator/=(const T& denom) noexcept
    {
        for (int i = 0; i < this->grid_size; ++i) {
            this->grid_data[i] /= denom;
        }

        return *this;
    }

private:

    int grid_size = 0;
    std::vector<T> grid_data{};

};

template<typename T_lhs, typename T_rhs>
constexpr bool operator==(const Regular_grid_1d<T_lhs>& lhs,
                          const Regular_grid_1d<T_rhs>& rhs) noexcept
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

    const int grid_size = lhs.size();

    for (int i = 0; i < grid_size; ++i) {
        if (!comp(lhs[i], rhs[i])) {
            return false;
        }
    }

    return true;
}

template<typename T_lhs, typename T_rhs>
constexpr Regular_grid_1d<decltype(std::declval<T_lhs>() +
                                   std::declval<T_rhs>())>
    operator+(const Regular_grid_1d<T_lhs>& lhs,
              const Regular_grid_1d<T_rhs>& rhs) noexcept
{
    if (rhs.size() == 0) {

        assert(lhs.size() != 0);
        return lhs;

    } else {

        assert(lhs.size() == rhs.size());

        const int grid_size = lhs.size();

        using res_type = decltype(std::declval<T_lhs>() +
                                  std::declval<T_rhs>());

        Regular_grid_1d<res_type> res(grid_size);

        for (int i = 0; i < grid_size; ++i) {
            res[i] = lhs[i] + rhs[i];
        }

        return res;

    }
}

template<typename T_lhs, typename T_rhs>
constexpr Regular_grid_1d<decltype(std::declval<T_lhs>() -
                                   std::declval<T_rhs>())>
    operator-(const Regular_grid_1d<T_lhs>& lhs,
              const Regular_grid_1d<T_rhs>& rhs) noexcept
{
    if (rhs.size() == 0) {

        assert(lhs.size() != 0);
        return lhs;

    } else {

        assert(lhs.size() == rhs.size());

        const int grid_size = lhs.size();

        using res_type = decltype(std::declval<T_lhs>() -
                                  std::declval<T_rhs>());

        Regular_grid_1d<res_type> res(grid_size);

        for (int i = 0; i < grid_size; ++i) {
            res[i] = lhs[i] - rhs[i];
        }

        return res;

    }
}

template<typename T_lhs, typename T_rhs>
constexpr Regular_grid_1d<decltype(std::declval<T_lhs>() *
                                   std::declval<T_rhs>())>
    operator*(const Regular_grid_1d<T_lhs>& lhs, const T_rhs& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() *
                              std::declval<T_rhs>());

    const int grid_size = lhs.size();

    Regular_grid_1d<res_type> res(grid_size);
    for (int i = 0; i < grid_size; ++i) {
        res[i] = lhs[i] * rhs;
    }

    return res;
}

template<typename T_lhs, typename T_rhs>
constexpr Regular_grid_1d<decltype(std::declval<T_lhs>() *
                                   std::declval<T_rhs>())>
    operator*(const T_lhs& lhs, const Regular_grid_1d<T_rhs>& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() *
                              std::declval<T_rhs>());

    const int grid_size = rhs.size();

    Regular_grid_1d<res_type> res(grid_size);
    for (int i = 0; i < grid_size; ++i) {
        res[i] = lhs * rhs[i];
    }

    return res;
}

template<typename T_lhs, typename T_rhs>
constexpr decltype(std::declval<decltype(std::declval<T_lhs>() * std::declval<T_lhs>())>() +
                   std::declval<decltype(std::declval<T_rhs>() * std::declval<T_rhs>())>())
    operator*(const Regular_grid_1d<T_lhs>& lhs,
              const Regular_grid_1d<T_rhs>& rhs) noexcept
{
    assert(lhs.size() == rhs.size());

    const int grid_size = lhs.size();

    using res_type = decltype(std::declval<
                                decltype(std::declval<T_lhs>() *
                                         std::declval<T_lhs>())
                                >() +
                              std::declval<
                                decltype(std::declval<T_rhs>() *
                                         std::declval<T_rhs>())
                                >());

    res_type res = res_type();

    for (int i = 0; i < grid_size; ++i) {
        res += lhs[i] * rhs[i];
    }

    return res;
}

template<typename T_lhs, typename T_rhs>
constexpr Regular_grid_1d<decltype(std::declval<T_lhs>() /
                                   std::declval<T_rhs>())>
    operator/(const Regular_grid_1d<T_lhs>& lhs, const T_rhs& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() /
                              std::declval<T_rhs>());

    const int grid_size = lhs.size();

    Regular_grid_1d<res_type> res(grid_size);
    for (int i = 0; i < grid_size; ++i) {
        res[i] = lhs[i] / rhs;
    }

    return res;
}

template<typename T>
constexpr auto max_abs(const Regular_grid_1d<T>& v) noexcept
{
    constexpr bool is_simple_standard_math_type =
                            std::is_arithmetic_v<std::decay_t<T>> ||
                            std::is_same_v<std::decay_t<T>, std::complex<float>> ||
                            std::is_same_v<std::decay_t<T>, std::complex<double>> ||
                            std::is_same_v<std::decay_t<T>, std::complex<long double>>;

    if constexpr (is_simple_standard_math_type) {
        const int grid_size = v.size();
        using max_value_t = decltype(std::abs(std::declval<decltype(v[0])>()));
        max_value_t max_value = max_value_t(0);
        for (int i = 0; i < grid_size; ++i) {
            max_value = std::max(max_value, std::abs(v[i]));
        }
        return max_value;
    } else {
        const int grid_size = v.size();
        using max_value_t = decltype(max_abs(std::declval<decltype(v[0])>()));
        max_value_t max_value = max_value_t();
        for (int i = 0; i < grid_size; ++i) {
            max_value = std::max(max_value, max_abs(v[i]));
        }
        return max_value;
    }
}

}

#endif // REGULAR_GRID_1D_HXX_INCLUDED
