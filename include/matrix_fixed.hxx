#ifndef MATRIX_FIXED_HXX_INCLUDED
#define MATRIX_FIXED_HXX_INCLUDED

#include "basic_libs.hxx"
#include "vector_fixed.hxx"
#include "basic_math_functions.hxx"

namespace Math_structures
{

template<
    typename T,
    int rows_amount,
    int cols_amount = rows_amount
    >
class Matrix_fixed
{
private:

    class Index_proxy_variable;
    class Index_proxy_constant;

    friend class Index_proxy_variable;
    friend class Index_proxy_constant;

public:

    static_assert(rows_amount > 0);
    static_assert(cols_amount > 0);

    using value_type = T;

    constexpr Matrix_fixed() noexcept = default;

    constexpr Matrix_fixed(const Matrix_fixed<T, rows_amount, cols_amount>& constr_matrix) noexcept
        :   matrix_data{constr_matrix.matrix_data}
    { }

    constexpr Matrix_fixed(Matrix_fixed<T, rows_amount, cols_amount>&& constr_matrix) noexcept
        :   matrix_data{std::move(constr_matrix.matrix_data)}
    { }

    constexpr Matrix_fixed(int constr_cols_amount,
                           int constr_rows_amount) noexcept
    {
        assert(constr_rows_amount == rows_amount);
        assert(constr_cols_amount == cols_amount);
    }

    constexpr Matrix_fixed(int constr_cols_amount,
                           int constr_rows_amount,
                           const T& initial_value) noexcept
        :   matrix_data{std::array<T, container_size>{}.fill(initial_value)}
    {
        assert(constr_rows_amount == rows_amount);
        assert(constr_cols_amount == cols_amount);
    }

    constexpr Matrix_fixed(const T (&data)[rows_amount][cols_amount]) noexcept
        :   matrix_data(container_size)
    {
        for (int i = 0; i < cols_amount; ++i) {
            for (int j = 0; j < rows_amount; ++j) {
                const int index = this->get_converted_index(i, j);
                this->matrix_data[index] = data[i][j];
            }
        }
    }

    constexpr Matrix_fixed(std::initializer_list<std::initializer_list<T>> data) noexcept
    {
        assert(data.size() == rows_amount);
        int n = 0;
        for (auto const& row: data) {
            assert(row.size() == cols_amount);
            for (auto const& elem: row) {
                this->matrix_data[n++] = elem;
            }
        }
    }

    ~Matrix_fixed() noexcept = default;

    constexpr Matrix_fixed<T, rows_amount, cols_amount>&
            operator=(const Matrix_fixed<T, rows_amount, cols_amount>& assign_matrix) noexcept
    {
        if (&assign_matrix != this) {
            this->matrix_data = assign_matrix.matrix_data;
        } else {}

        return *this;
    }

    constexpr Matrix_fixed<T, rows_amount, cols_amount>&
            operator=(Matrix_fixed<T, rows_amount, cols_amount>&& assign_matrix) noexcept
    {
        if (&assign_matrix != this) {
            this->matrix_data = std::move(assign_matrix.matrix_data);
        } else {}

        return *this;
    }

    constexpr Index_proxy_variable operator[](const int row_number) noexcept
    {
        assert(row_number >= 0 && row_number < rows_amount);
        return Index_proxy_variable{*this, row_number};
    }

    constexpr Index_proxy_constant operator[](const int row_number) const noexcept
    {
        assert(row_number >= 0 && row_number < rows_amount);
        return Index_proxy_constant{*this, row_number};
    }

    constexpr int get_rows_amount() const noexcept
    {
        return rows_amount;
    }

    constexpr int get_cols_amount() const noexcept
    {
        return cols_amount;
    }

    constexpr void set_null_value() noexcept
    {
        for (auto& elem: this->matrix_data) {

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

    constexpr Matrix_fixed<T, rows_amount, cols_amount>&
            operator+=(const Matrix_fixed<T, rows_amount, cols_amount>& add_matrix) noexcept
    {
        for (int i = 0; i < container_size; ++i) {
            this->matrix_data[i] += add_matrix.matrix_data[i];
        }

        return *this;
    }

    constexpr Matrix_fixed<T, rows_amount, cols_amount>&
            operator-=(const Matrix_fixed<T, rows_amount, cols_amount>& sub_matrix) noexcept
    {
        for (int i = 0; i < container_size; ++i) {
            this->matrix_data[i] -= sub_matrix.matrix_data[i];
        }

        return *this;
    }

    constexpr Matrix_fixed<T, rows_amount, cols_amount>& operator*=(const T& multiplier) noexcept
    {
        for (int i = 0; i < container_size; ++i) {
            this->matrix_data[i] *= multiplier;
        }

        return *this;
    }

    constexpr Matrix_fixed<T, rows_amount, cols_amount>& operator/=(const T& denom) noexcept
    {
        for (int i = 0; i < container_size; ++i) {
            this->matrix_data[i] /= denom;
        }

        return *this;
    }

private:

    constexpr int get_converted_index(const int row_number,
                                      const int col_number) const noexcept
    {
        assert(row_number >= 0 && row_number < rows_amount);
        assert(col_number >= 0 && col_number < cols_amount);

        const int index = row_number * cols_amount + col_number;

        assert(index >= 0 && index < container_size);

        return index;
    }

    static constexpr const int container_size = rows_amount * cols_amount;

    std::array<T, container_size> matrix_data{};

};

template<
    typename T,
    int rows_amount,
    int cols_amount
    >
class Matrix_fixed<T, rows_amount, cols_amount>::Index_proxy_variable
{
public:

    constexpr Index_proxy_variable() noexcept = default;

    constexpr Index_proxy_variable(const Index_proxy_variable&) noexcept = default;
    constexpr Index_proxy_variable(Index_proxy_variable&&) noexcept = default;

    constexpr Index_proxy_variable(Matrix_fixed<T, rows_amount, cols_amount>& matrix,
                                   const int row_number) noexcept
        :   matrix{matrix},
            row_number{row_number}
    {
        assert(this->row_number >= 0 && this->row_number < rows_amount);
    }

    ~Index_proxy_variable() noexcept = default;

    constexpr T& operator[](const int col_number) noexcept
    {
        assert(col_number >= 0 && col_number < cols_amount);
        const int index = this->matrix.get_converted_index(row_number, col_number);
        return this->matrix.matrix_data[index];
    }

    constexpr operator std::array<T, cols_amount>() noexcept
    {
        std::array<T, cols_amount> res{};
        for (int n = 0; n < cols_amount; ++n) {
            const int index = this->matrix.get_converted_index(this->row_number, n);
            res[n] = this->matrix.matrix_data[index];
        }
        return res;
    }

private:

    Matrix_fixed<T, rows_amount, cols_amount>& matrix;
    const int row_number = 0;

};

template<
    typename T,
    int rows_amount,
    int cols_amount
    >
class Matrix_fixed<T, rows_amount, cols_amount>::Index_proxy_constant
{
public:

    constexpr Index_proxy_constant() noexcept = default;

    constexpr Index_proxy_constant(const Index_proxy_constant&) noexcept = default;
    constexpr Index_proxy_constant(Index_proxy_constant&&) noexcept = default;

    constexpr Index_proxy_constant(const Matrix_fixed<T, rows_amount, cols_amount>& matrix,
                                   const int row_number) noexcept
        :   matrix{matrix},
            row_number{row_number}
    {
        assert(this->row_number >= 0 && this->row_number < rows_amount);
    }

    ~Index_proxy_constant() noexcept = default;

    constexpr const T& operator[](const int col_number) const noexcept
    {
        assert(col_number >= 0 && col_number < cols_amount);
        const int index = this->matrix.get_converted_index(this->row_number, col_number);
        return this->matrix.matrix_data[index];
    }

    constexpr operator std::array<T, cols_amount>() const noexcept
    {
        std::array<T, cols_amount> res{};
        for (int n = 0; n < cols_amount; ++n) {
            const int index = this->matrix.get_converted_index(this->row_number, n);
            res[n] = this->matrix.matrix_data[index];
        }
        return res;
    }

private:

    const Matrix_fixed<T, rows_amount, cols_amount>& matrix;
    const int row_number = 0;

};

template<
    typename T_lhs,
    typename T_rhs,
    int rows_amount,
    int cols_amount
    >
constexpr bool operator==(const Matrix_fixed<T_lhs, rows_amount, cols_amount>& lhs,
                          const Matrix_fixed<T_rhs, rows_amount, cols_amount>& rhs) noexcept
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

    for (int i = 0; i < rows_amount; ++i) {
        for (int j = 0; j < cols_amount; ++j) {
            if (!comp(lhs[i][j], rhs[i][j])) {
                return false;
            }
        }
    }

    return true;
}

template<
    typename T_lhs,
    typename T_rhs,
    int rows_amount,
    int cols_amount
    >
constexpr Matrix_fixed<
        decltype(std::declval<T_lhs>() +
                 std::declval<T_rhs>()),
        rows_amount,
        cols_amount
        >
    operator+(const Matrix_fixed<T_lhs, rows_amount, cols_amount>& lhs,
              const Matrix_fixed<T_rhs, rows_amount, cols_amount>& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() +
                              std::declval<T_rhs>());

    Matrix_fixed<res_type, rows_amount, cols_amount> res{};

    for (int i = 0; i < rows_amount; ++i) {
        for (int j = 0; j < cols_amount; ++j) {
            res[i][j] = lhs[i][j] + rhs[i][j];
        }
    }

    return res;
}

template<
    typename T_lhs,
    typename T_rhs,
    int rows_amount,
    int cols_amount
    >
constexpr Matrix_fixed<
        decltype(std::declval<T_lhs>() -
                 std::declval<T_rhs>()),
        rows_amount,
        cols_amount
        >
    operator-(const Matrix_fixed<T_lhs, rows_amount, cols_amount>& lhs,
              const Matrix_fixed<T_rhs, rows_amount, cols_amount>& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() -
                              std::declval<T_rhs>());

    Matrix_fixed<res_type, rows_amount, cols_amount> res{};

    for (int i = 0; i < rows_amount; ++i) {
        for (int j = 0; j < cols_amount; ++j) {
            res[i][j] = lhs[i][j] - rhs[i][j];
        }
    }

    return res;
}

template<
    typename T_lhs,
    typename T_rhs,
    int rows_amount,
    int cols_amount
    >
constexpr Matrix_fixed<
        decltype(std::declval<T_lhs>() *
                 std::declval<T_rhs>()),
        rows_amount,
        cols_amount
        >
    operator*(const Matrix_fixed<T_lhs, rows_amount, cols_amount>& lhs, const T_rhs& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() *
                              std::declval<T_rhs>());

    Matrix_fixed<res_type, rows_amount, cols_amount> res{};

    for (int i = 0; i < rows_amount; ++i) {
        for (int j = 0; j < cols_amount; ++j) {
            res[i][j] = lhs[i][j] * rhs;
        }
    }

    return res;
}

template<
    typename T_lhs,
    typename T_rhs,
    int rows_amount,
    int cols_amount
    >
constexpr Matrix_fixed<
        decltype(std::declval<T_lhs>() *
                 std::declval<T_rhs>()),
        rows_amount,
        cols_amount
        >
    operator*(const T_lhs& lhs, const Matrix_fixed<T_rhs, rows_amount, cols_amount>& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() *
                              std::declval<T_rhs>());

    Matrix_fixed<res_type, rows_amount, cols_amount> res{};

    for (int i = 0; i < rows_amount; ++i) {
        for (int j = 0; j < cols_amount; ++j) {
            res[i][j] = lhs * rhs[i][j];
        }
    }

    return res;
}

template<
    typename T_lhs,
    typename T_rhs,
    int rows_amount_lhs,
    int cols_amount_lhs,
    int rows_amount_rhs,
    int cols_amount_rhs
    >
constexpr Matrix_fixed<
        decltype(std::declval<T_lhs>() *
                 std::declval<T_rhs>()),
        rows_amount_lhs,
        cols_amount_rhs
        >
    operator*(const Matrix_fixed<T_lhs, rows_amount_lhs, cols_amount_lhs>& lhs,
              const Matrix_fixed<T_rhs, rows_amount_rhs, cols_amount_rhs>& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() *
                              std::declval<T_rhs>());

    static_assert(cols_amount_lhs == rows_amount_rhs);

    Matrix_fixed<res_type, rows_amount_lhs, cols_amount_rhs> res{};

    constexpr bool is_simple_standard_math_type =
                            std::is_arithmetic_v<std::decay_t<res_type>> ||
                            std::is_same_v<std::decay_t<res_type>, std::complex<float>> ||
                            std::is_same_v<std::decay_t<res_type>, std::complex<double>> ||
                            std::is_same_v<std::decay_t<res_type>, std::complex<long double>>;

    if constexpr (is_simple_standard_math_type) {
        for (int i = 0; i < rows_amount_lhs; ++i) {
            for (int j = 0; j < cols_amount_rhs; ++j) {
                res[i][j] = res_type(0);
            }
        }
    } else {
        res.set_null_value();
    }

    for (int i = 0; i < rows_amount_lhs; ++i) {
        for (int j = 0; j < cols_amount_rhs; ++j) {
            for (int k = 0; k < cols_amount_lhs; ++k) {
                res[i][j] += lhs[i][k] * rhs[k][j];
            }
        }
    }

    return res;
}

template<
    typename T_lhs,
    typename T_rhs,
    int rows_amount_matrix,
    int cols_amount_matrix
    >
constexpr Vector_fixed<
        decltype(std::declval<T_lhs>() *
                 std::declval<T_rhs>()),
        rows_amount_matrix
        >
    operator*(const Matrix_fixed<T_lhs, rows_amount_matrix, cols_amount_matrix>& lhs,
              const Vector_fixed<T_rhs, cols_amount_matrix>& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() *
                              std::declval<T_rhs>());

    Vector_fixed<res_type, rows_amount_matrix> res{};

    constexpr bool is_simple_standard_math_type =
                            std::is_arithmetic_v<std::decay_t<res_type>> ||
                            std::is_same_v<std::decay_t<res_type>, std::complex<float>> ||
                            std::is_same_v<std::decay_t<res_type>, std::complex<double>> ||
                            std::is_same_v<std::decay_t<res_type>, std::complex<long double>>;

    if constexpr (is_simple_standard_math_type) {
        for (int i = 0; i < rows_amount_matrix; ++i) {
            res[i] = res_type(0);
        }
    } else {
        res.set_null_value();
    }

    for (int i = 0; i < rows_amount_matrix; ++i) {
        for (int j = 0; j < cols_amount_matrix; ++j) {
            res[i] += lhs[i][j] * rhs[j];
        }
    }

    return res;
}

template<
    typename T_lhs,
    typename T_rhs,
    int rows_amount,
    int cols_amount
    >
constexpr Matrix_fixed<
        decltype(std::declval<T_lhs>() /
                 std::declval<T_rhs>()),
        rows_amount,
        cols_amount
        >
    operator/(const Matrix_fixed<T_lhs, rows_amount, cols_amount>& lhs, const T_rhs& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() /
                              std::declval<T_rhs>());

    Matrix_fixed<res_type, rows_amount, cols_amount> res{};

    for (int i = 0; i < rows_amount; ++i) {
        for (int j = 0; j < cols_amount; ++j) {
            res[i][j] = lhs[i][j] / rhs;
        }
    }

    return res;
}

template<
    typename T,
    int rows_amount,
    int cols_amount
    >
constexpr auto max_abs(const Matrix_fixed<T, rows_amount, cols_amount>& m) noexcept
{
    constexpr bool is_simple_standard_math_type =
                            std::is_arithmetic_v<std::decay_t<T>> ||
                            std::is_same_v<std::decay_t<T>, std::complex<float>> ||
                            std::is_same_v<std::decay_t<T>, std::complex<double>> ||
                            std::is_same_v<std::decay_t<T>, std::complex<long double>>;

    if constexpr (is_simple_standard_math_type) {

        using max_value_t = decltype(std::abs(std::declval<decltype(m[0][0])>()));

        max_value_t max_value = max_value_t(0);

        for (int i = 0; i < rows_amount; ++i) {
            for (int j = 0; j < cols_amount; ++j) {
                max_value = std::max(max_value, std::abs(m[i][j]));
            }
        }

        return max_value;

    } else {

        using max_value_t = decltype(max_abs(std::declval<decltype(m[0][0])>()));

        max_value_t max_value = max_value_t();

        for (int i = 0; i < rows_amount; ++i) {
            for (int j = 0; j < cols_amount; ++j) {
                max_value = std::max(max_value, max_abs(m[i][j]));
            }
        }

        return max_value;

    }
}

}

#endif // MATRIX_FIXED_HXX_INCLUDED
