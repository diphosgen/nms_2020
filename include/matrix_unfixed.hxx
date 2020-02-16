#ifndef MATRIX_UNFIXED_HXX_INCLUDED
#define MATRIX_UNFIXED_HXX_INCLUDED

#include "basic_libs.hxx"
#include "vector_unfixed.hxx"
#include "basic_math_functions.hxx"

namespace Math_structures
{

template<typename T = double>
class Matrix_unfixed
{
private:

    class Index_proxy_variable;
    class Index_proxy_constant;

    friend class Index_proxy_variable;
    friend class Index_proxy_constant;

public:

    using value_type = T;

    constexpr Matrix_unfixed() noexcept = default;

    constexpr Matrix_unfixed(const Matrix_unfixed<T>& constr_matrix) noexcept
        :   rows_amount{constr_matrix.rows_amount},
            cols_amount{constr_matrix.cols_amount},
            container_size{constr_matrix.container_size},
            matrix_data{constr_matrix.matrix_data}
    { }

    constexpr Matrix_unfixed(Matrix_unfixed<T>&& constr_matrix) noexcept
        :   rows_amount{std::move(constr_matrix.rows_amount)},
            cols_amount{std::move(constr_matrix.cols_amount)},
            container_size{std::move(constr_matrix.container_size)},
            matrix_data{std::move(constr_matrix.matrix_data)}
    { }

    constexpr Matrix_unfixed(int cols_amount, int rows_amount) noexcept
        :   rows_amount{rows_amount},
            cols_amount{cols_amount},
            container_size{rows_amount * cols_amount},
            matrix_data(container_size)
    { }

    constexpr Matrix_unfixed(int cols_amount, int rows_amount, const T& initial_value) noexcept
        :   rows_amount{rows_amount},
            cols_amount{cols_amount},
            container_size{rows_amount * cols_amount},
            matrix_data(container_size, initial_value)
    { }

    template<
        int rows_amount,
        int cols_amount
        >
    constexpr Matrix_unfixed(const T (&data)[rows_amount][cols_amount]) noexcept
        :   rows_amount{rows_amount},
            cols_amount{cols_amount},
            container_size{rows_amount * cols_amount},
            matrix_data(container_size)
    {
        for (int i = 0; i < cols_amount; ++i) {
            for (int j = 0; j < rows_amount; ++j) {
                const int index = this->get_converted_index(i, j);
                this->matrix_data[index] = data[i][j];
            }
        }
    }

    constexpr Matrix_unfixed(std::initializer_list<std::initializer_list<T>> data) noexcept
        :   rows_amount(data.size()),
            cols_amount(data.begin()->size()),
            container_size(rows_amount * cols_amount),
            matrix_data(container_size)
    {
        int n = 0;
        for (auto const& row: data) {
            assert(int(row.size()) == this->cols_amount);
            for (auto const& elem: row) {
                this->matrix_data[n++] = elem;
            }
        }
    }

    ~Matrix_unfixed() noexcept = default;

    constexpr Matrix_unfixed<T>& operator=(const Matrix_unfixed<T>& assign_matrix) noexcept
    {
        if (&assign_matrix != this) {
            if (this->rows_amount == 0 &&
                this->cols_amount == 0) {

                this->rows_amount = assign_matrix.rows_amount;
                this->cols_amount = assign_matrix.cols_amount;
                this->container_size = assign_matrix.container_size;
                this->matrix_data = assign_matrix.matrix_data;

            } else {

                assert(this->rows_amount == assign_matrix.rows_amount);
                assert(this->cols_amount == assign_matrix.cols_amount);
                assert(this->container_size == assign_matrix.container_size);

                this->matrix_data = assign_matrix.matrix_data;

            }
        } else {}

        return *this;
    }

    constexpr Matrix_unfixed<T>& operator=(Matrix_unfixed<T>&& assign_matrix) noexcept
    {
        if (&assign_matrix != this) {
            if (this->rows_amount == 0 &&
                this->cols_amount == 0) {

                this->rows_amount = std::move(assign_matrix.rows_amount);
                this->cols_amount = std::move(assign_matrix.cols_amount);
                this->container_size = std::move(assign_matrix.container_size);
                this->matrix_data = std::move(assign_matrix.matrix_data);

            } else {

                assert(this->rows_amount == assign_matrix.rows_amount);
                assert(this->cols_amount == assign_matrix.cols_amount);
                assert(this->container_size == assign_matrix.container_size);

                this->matrix_data = std::move(assign_matrix.matrix_data);

            }
        } else {}

        return *this;
    }

    constexpr Index_proxy_variable operator[](const int row_number) noexcept
    {
        assert(row_number >= 0 && row_number < this->rows_amount);
        return Index_proxy_variable{*this, row_number};
    }

    constexpr Index_proxy_constant operator[](const int row_number) const noexcept
    {
        assert(row_number >= 0 && row_number < this->rows_amount);
        return Index_proxy_constant{*this, row_number};
    }

    constexpr int get_rows_amount() const noexcept
    {
        return this->rows_amount;
    }

    constexpr int get_cols_amount() const noexcept
    {
        return this->cols_amount;
    }

    constexpr void set_null_value() noexcept
    {
        if (this->container_size != 0) {
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
        } else {}
    }

    constexpr Matrix_unfixed<T>& operator+=(const Matrix_unfixed<T>& add_matrix) noexcept
    {
        if (this->rows_amount == 0 &&
            this->cols_amount == 0) {

            this->rows_amount = add_matrix.rows_amount;
            this->cols_amount = add_matrix.cols_amount;
            this->container_size = add_matrix.container_size;
            this->matrix_data = add_matrix.matrix_data;

        } else {

            assert(this->rows_amount == add_matrix.rows_amount);
            assert(this->cols_amount == add_matrix.cols_amount);
            assert(this->container_size == add_matrix.container_size);

            for (int i = 0; i < this->container_size; ++i) {
                this->matrix_data[i] += add_matrix.matrix_data[i];
            }
        }

        return *this;
    }

    constexpr Matrix_unfixed<T>& operator-=(const Matrix_unfixed<T>& sub_matrix) noexcept
    {
        if (this->rows_amount == 0 &&
            this->cols_amount == 0) {

            this->rows_amount = sub_matrix.rows_amount;
            this->cols_amount = sub_matrix.cols_amount;
            this->container_size = sub_matrix.container_size;

            this->set_null_value();

        } else {

            assert(this->rows_amount == sub_matrix.rows_amount);
            assert(this->cols_amount == sub_matrix.cols_amount);
            assert(this->container_size == sub_matrix.container_size);

        }

        for (int i = 0; i < this->container_size; ++i) {
            this->matrix_data[i] -= sub_matrix.matrix_data[i];
        }

        return *this;
    }

    constexpr Matrix_unfixed<T>& operator*=(const T& multiplier) noexcept
    {
        for (int i = 0; i < this->container_size; ++i) {
            this->matrix_data[i] *= multiplier;
        }

        return *this;
    }

    constexpr Matrix_unfixed<T>& operator/=(const T& denom) noexcept
    {
        for (int i = 0; i < this->container_size; ++i) {
            this->matrix_data[i] /= denom;
        }

        return *this;
    }

private:

    constexpr int get_converted_index(const int row_number,
                                      const int col_number) const noexcept
    {
        assert(row_number >= 0 && row_number < this->rows_amount);
        assert(col_number >= 0 && col_number < this->cols_amount);

        const int index = row_number * this->cols_amount + col_number;

        assert(index >= 0 && index < this->container_size);

        return index;
    }

    int rows_amount = 0;
    int cols_amount = 0;

    int container_size = 0;

    std::vector<T> matrix_data{};

};

template<typename T>
class Matrix_unfixed<T>::Index_proxy_variable
{
public:

    constexpr Index_proxy_variable() noexcept = default;

    constexpr Index_proxy_variable(const Index_proxy_variable&) noexcept = default;
    constexpr Index_proxy_variable(Index_proxy_variable&&) noexcept = default;

    constexpr Index_proxy_variable(Matrix_unfixed<T>& matrix, const int row_number) noexcept
        :   matrix{matrix},
            row_number{row_number}
    {
        assert(this->matrix.container_size != 0);
        assert(this->row_number >= 0 && this->row_number < this->matrix.rows_amount);
    }

    ~Index_proxy_variable() noexcept = default;

    constexpr Index_proxy_variable&
		operator=(const Index_proxy_variable&) noexcept = default;

    constexpr Index_proxy_variable&
		operator=(Index_proxy_variable&&) noexcept = default;

    constexpr T& operator[](const int col_number) noexcept
    {
        assert(col_number >= 0 && col_number < this->matrix.cols_amount);
        const int index = this->matrix.get_converted_index(row_number, col_number);
        return this->matrix.matrix_data[index];
    }

    constexpr operator std::vector<T>() noexcept
    {
        std::vector<T> res(this->matrix.cols_amount);
        for (int n = 0; n < this->matrix.cols_amount; ++n) {
            const int index = this->matrix.get_converted_index(this->row_number, n);
            res[n] = this->matrix.matrix_data[index];
        }
        return res;
    }

private:

    Matrix_unfixed<T>& matrix;
    const int row_number = 0;

};

template<typename T>
class Matrix_unfixed<T>::Index_proxy_constant
{
public:

    constexpr Index_proxy_constant() noexcept = default;

    constexpr Index_proxy_constant(const Index_proxy_constant&) noexcept = default;
    constexpr Index_proxy_constant(Index_proxy_constant&&) noexcept = default;

    constexpr Index_proxy_constant(const Matrix_unfixed<T>& matrix, const int row_number) noexcept
        :   matrix{matrix},
            row_number{row_number}
    {
        assert(this->matrix.container_size != 0);
        assert(this->row_number >= 0 && this->row_number < this->matrix.rows_amount);
    }

    ~Index_proxy_constant() noexcept = default;

    constexpr Index_proxy_constant&
		operator=(const Index_proxy_constant&) noexcept = default;

    constexpr Index_proxy_constant&
		operator=(Index_proxy_constant&&) noexcept = default;

    constexpr const T& operator[](const int col_number) const noexcept
    {
        assert(col_number >= 0 && col_number < this->matrix.cols_amount);
        const int index = this->matrix.get_converted_index(this->row_number, col_number);
        return this->matrix.matrix_data[index];
    }

    constexpr operator std::vector<T>() const noexcept
    {
        std::vector<T> res(this->matrix.cols_amount);
        for (int n = 0; n < this->matrix.cols_amount; ++n) {
            const int index = this->matrix.get_converted_index(this->row_number, n);
            res[n] = this->matrix.matrix_data[index];
        }
        return res;
    }

private:

    const Matrix_unfixed<T>& matrix;
    const int row_number = 0;

};

template<typename T_lhs, typename T_rhs>
constexpr bool operator==(const Matrix_unfixed<T_lhs>& lhs,
                          const Matrix_unfixed<T_rhs>& rhs) noexcept
{
    assert(lhs.get_rows_amount() == rhs.get_rows_amount());
    assert(lhs.get_cols_amount() == rhs.get_cols_amount());

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

    const int rows_amount = lhs.get_rows_amount();
    const int cols_amount = lhs.get_cols_amount();

    for (int i = 0; i < rows_amount; ++i) {
        for (int j = 0; j < cols_amount; ++j) {
            if (!comp(lhs[i][j], rhs[i][j])) {
                return false;
            }
        }
    }

    return true;
}

template<typename T_lhs, typename T_rhs>
constexpr Matrix_unfixed<decltype(std::declval<T_lhs>() +
                                  std::declval<T_rhs>())>
    operator+(const Matrix_unfixed<T_lhs>& lhs,
              const Matrix_unfixed<T_rhs>& rhs) noexcept
{
    if (rhs.get_rows_amount() == 0 && rhs.get_cols_amount() == 0) {

        assert(lhs.get_rows_amount() != 0);
        assert(lhs.get_cols_amount() != 0);

        return lhs;

    } else {

        assert(lhs.get_rows_amount() == rhs.get_rows_amount());
        assert(lhs.get_cols_amount() == rhs.get_cols_amount());

        const int rows_amount = lhs.get_rows_amount();
        const int cols_amount = lhs.get_cols_amount();

        using res_type = decltype(std::declval<T_lhs>() +
                                  std::declval<T_rhs>());

        Matrix_unfixed<res_type> res(rows_amount, cols_amount);

        for (int i = 0; i < rows_amount; ++i) {
            for (int j = 0; j < cols_amount; ++j) {
                res[i][j] = lhs[i][j] + rhs[i][j];
            }
        }

        return res;

    }
}

template<typename T_lhs, typename T_rhs>
constexpr Matrix_unfixed<decltype(std::declval<T_lhs>() -
                                  std::declval<T_rhs>())>
    operator-(const Matrix_unfixed<T_lhs>& lhs,
              const Matrix_unfixed<T_rhs>& rhs) noexcept
{
    if (rhs.get_rows_amount() == 0 && rhs.get_cols_amount() == 0) {

        assert(lhs.get_rows_amount() != 0);
        assert(lhs.get_cols_amount() != 0);

        return lhs;

    } else {

        assert(lhs.get_rows_amount() == rhs.get_rows_amount());
        assert(lhs.get_cols_amount() == rhs.get_cols_amount());

        const int rows_amount = lhs.get_rows_amount();
        const int cols_amount = lhs.get_cols_amount();

        using res_type = decltype(std::declval<T_lhs>() -
                                  std::declval<T_rhs>());

        Matrix_unfixed<res_type> res(rows_amount, cols_amount);

        for (int i = 0; i < rows_amount; ++i) {
            for (int j = 0; j < cols_amount; ++j) {
                res[i][j] = lhs[i][j] - rhs[i][j];
            }
        }

        return res;

    }
}

template<typename T_lhs, typename T_rhs>
constexpr Matrix_unfixed<decltype(std::declval<T_lhs>() *
                                  std::declval<T_rhs>())>
    operator*(const Matrix_unfixed<T_lhs>& lhs, const T_rhs& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() *
                              std::declval<T_rhs>());

    const int rows_amount = lhs.get_rows_amount();
    const int cols_amount = lhs.get_cols_amount();

    Matrix_unfixed<res_type> res(rows_amount, cols_amount);

    for (int i = 0; i < rows_amount; ++i) {
        for (int j = 0; j < cols_amount; ++j) {
            res[i][j] = lhs[i][j] * rhs;
        }
    }

    return res;
}

template<typename T_lhs, typename T_rhs>
constexpr Matrix_unfixed<decltype(std::declval<T_lhs>() *
                                  std::declval<T_rhs>())>
    operator*(const T_lhs& lhs, const Matrix_unfixed<T_rhs>& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() *
                              std::declval<T_rhs>());

    const int rows_amount = rhs.get_rows_amount();
    const int cols_amount = rhs.get_cols_amount();

    Matrix_unfixed<res_type> res(rows_amount, cols_amount);

    for (int i = 0; i < rows_amount; ++i) {
        for (int j = 0; j < cols_amount; ++j) {
            res[i][j] = lhs * rhs[i][j];
        }
    }

    return res;
}

template<typename T_lhs, typename T_rhs>
constexpr Matrix_unfixed<decltype(std::declval<T_lhs>() *
                                  std::declval<T_rhs>())>
    operator*(const Matrix_unfixed<T_lhs>& lhs,
              const Matrix_unfixed<T_rhs>& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() *
                              std::declval<T_rhs>());

    const int rows_amount_lhs = lhs.get_rows_amount();
    const int cols_amount_lhs = lhs.get_cols_amount();

    const int rows_amount_rhs = rhs.get_rows_amount();
    const int cols_amount_rhs = rhs.get_cols_amount();

    assert(cols_amount_lhs == rows_amount_rhs);

    Matrix_unfixed<res_type> res(rows_amount_lhs,
                                 cols_amount_rhs);

    if constexpr (std::is_arithmetic_v<std::decay_t<res_type>>) {
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

template<typename T_lhs, typename T_rhs>
constexpr Vector_unfixed<decltype(std::declval<T_lhs>() *
                                  std::declval<T_rhs>())>
    operator*(const Matrix_unfixed<T_lhs>& lhs,
              const Vector_unfixed<T_rhs>& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() *
                              std::declval<T_rhs>());

    const int rows_amount_lhs = lhs.get_rows_amount();
    const int cols_amount_lhs = lhs.get_cols_amount();

    const int vector_size = rhs.size();

    assert(cols_amount_lhs == vector_size);

    Vector_unfixed<res_type> res(rows_amount_lhs);

    constexpr bool is_simple_standard_math_type =
                            std::is_arithmetic_v<std::decay_t<res_type>> ||
                            std::is_same_v<std::decay_t<res_type>, std::complex<float>> ||
                            std::is_same_v<std::decay_t<res_type>, std::complex<double>> ||
                            std::is_same_v<std::decay_t<res_type>, std::complex<long double>>;

    if constexpr (is_simple_standard_math_type) {
        for (int i = 0; i < rows_amount_lhs; ++i) {
            res[i] = res_type(0);
        }
    } else {
        res.set_null_value();
    }

    for (int i = 0; i < rows_amount_lhs; ++i) {
        for (int j = 0; j < vector_size; ++j) {
            res[i] += lhs[i][j] * rhs[j];
        }
    }

    return res;
}

template<typename T_lhs, typename T_rhs>
constexpr Matrix_unfixed<decltype(std::declval<T_lhs>() /
                                  std::declval<T_rhs>())>
    operator/(const Matrix_unfixed<T_lhs>& lhs, const T_rhs& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() /
                              std::declval<T_rhs>());

    const int rows_amount = lhs.get_rows_amount();
    const int cols_amount = lhs.get_cols_amount();

    Matrix_unfixed<res_type> res(rows_amount, cols_amount);

    for (int i = 0; i < rows_amount; ++i) {
        for (int j = 0; j < cols_amount; ++j) {
            res[i][j] = lhs[i][j] / rhs;
        }
    }

    return res;
}

template<typename T>
constexpr auto max_abs(const Matrix_unfixed<T>& m) noexcept
{
    constexpr bool is_simple_standard_math_type =
                            std::is_arithmetic_v<std::decay_t<T>> ||
                            std::is_same_v<std::decay_t<T>, std::complex<float>> ||
                            std::is_same_v<std::decay_t<T>, std::complex<double>> ||
                            std::is_same_v<std::decay_t<T>, std::complex<long double>>;

    if constexpr (is_simple_standard_math_type) {

        const int rows_amount = m.get_rows_amount();
        const int cols_amount = m.get_cols_amount();

        using max_value_t = decltype(std::abs(std::declval<decltype(m[0][0])>()));

        max_value_t max_value = max_value_t(0);

        for (int i = 0; i < rows_amount; ++i) {
            for (int j = 0; j < cols_amount; ++j) {
                max_value = std::max(max_value, std::abs(m[i][j]));
            }
        }

        return max_value;

    } else {

        const int rows_amount = m.get_rows_amount();
        const int cols_amount = m.get_cols_amount();

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

#endif // MATRIX_UNFIXED_HXX_INCLUDED
