#ifndef REGULAR_GRID_3D_HXX_INCLUDED
#define REGULAR_GRID_3D_HXX_INCLUDED

#include "basic_libs.hxx"
#include "basic_math_functions.hxx"

namespace Math_structures
{

template<typename T = double>
class Regular_grid_3d
{
private:

    class Index_proxy_variable;
    class Index_proxy_constant;

    friend class Index_proxy_variable;
    friend class Index_proxy_constant;

public:

    using value_type = T;

    constexpr Regular_grid_3d() noexcept = default;

    constexpr Regular_grid_3d(const Regular_grid_3d<T>& constr_grid) noexcept
        :   size_x{constr_grid.size_x},
            size_y{constr_grid.size_y},
            size_z{constr_grid.size_z},
            container_size{constr_grid.container_size},
            grid_data{constr_grid.grid_data}
    { }

    constexpr Regular_grid_3d(Regular_grid_3d<T>&& constr_grid) noexcept
        :   size_x{std::move(constr_grid.size_x)},
            size_y{std::move(constr_grid.size_y)},
            size_z{std::move(constr_grid.size_z)},
            container_size{std::move(constr_grid.container_size)},
            grid_data{std::move(constr_grid.grid_data)}
    { }

    constexpr Regular_grid_3d(int size_y, int size_x, int size_z) noexcept
        :   size_x{size_x},
            size_y{size_y},
            size_z{size_z},
            container_size{size_x * size_y * size_z},
            grid_data(container_size)
    { }

    constexpr Regular_grid_3d(int size_y, int size_x, int size_z, const T& initial_value) noexcept
        :   size_x{size_x},
            size_y{size_y},
            size_z{size_z},
            container_size{size_x * size_y * size_z},
            grid_data(container_size, initial_value)
    { }

    template<
        int size_x,
        int size_y,
        int size_z
        >
    constexpr Regular_grid_3d(const T (&data)[size_x][size_y][size_z]) noexcept
        :   size_x{size_x},
            size_y{size_y},
            size_z{size_z},
            container_size{size_x * size_y * size_z},
            grid_data(container_size)
    {
        for (int i = 0; i < size_y; ++i) {
            for (int j = 0; j < size_x; ++j) {
                for (int k = 0; k < size_z; ++k) {
                    const int index = this->get_converted_index(i, j, k);
                    this->grid_data[index] = data[i][j][k];
                }
            }
        }
    }

    constexpr Regular_grid_3d(std::initializer_list<
                                std::initializer_list<
                                    std::initializer_list<T>
                                    >
                                > data) noexcept
        :   size_x(data.size()),
            size_y(data.begin()->size()),
            size_z(data.begin()->begin()->size()),
            container_size(size_x * size_y * size_z),
            grid_data(container_size)
    {
        int n = 0;
        for (auto const& data_x: data) {
            assert(int(data_x.size()) == this->size_y);
            for (auto const& data_y: data_x) {
                assert(int(data_y.size()) == this->size_z);
                for (auto const& elem: data_y) {
                    this->grid_data[n++] = elem;
                }
            }
        }
    }

    ~Regular_grid_3d() noexcept = default;

    constexpr Regular_grid_3d<T>& operator=(const Regular_grid_3d<T>& assign_grid) noexcept
    {
        if (&assign_grid != this) {
            if (this->size_x == 0 &&
                this->size_y == 0 &&
                this->size_z == 0) {

                this->size_x = assign_grid.size_x;
                this->size_y = assign_grid.size_y;
                this->size_z = assign_grid.size_z;
                this->container_size = assign_grid.container_size;
                this->grid_data = assign_grid.grid_data;

            } else {

                assert(this->size_x == assign_grid.size_x);
                assert(this->size_y == assign_grid.size_y);
                assert(this->size_z == assign_grid.size_z);
                assert(this->container_size == assign_grid.container_size);

                this->grid_data = assign_grid.grid_data;

            }
        } else {}

        return *this;
    }

    constexpr Regular_grid_3d<T>& operator=(Regular_grid_3d<T>&& assign_grid) noexcept
    {
        if (&assign_grid != this) {
            if (this->size_x == 0 &&
                this->size_y == 0 &&
                this->size_z == 0) {

                this->size_x = std::move(assign_grid.size_x);
                this->size_y = std::move(assign_grid.size_y);
                this->size_z = std::move(assign_grid.size_z);
                this->container_size = std::move(assign_grid.container_size);
                this->grid_data = std::move(assign_grid.grid_data);

            } else {

                assert(this->size_x == assign_grid.size_x);
                assert(this->size_y == assign_grid.size_y);
                assert(this->size_z == assign_grid.size_z);
                assert(this->container_size == assign_grid.container_size);

                this->grid_data = std::move(assign_grid.grid_data);

            }
        } else {}

        return *this;
    }

    constexpr Index_proxy_variable operator[](const int x_node_number) noexcept
    {
        assert(x_node_number >= 0 && x_node_number < this->size_x);
        return Index_proxy_variable{*this, x_node_number};
    }

    constexpr Index_proxy_constant operator[](const int x_node_number) const noexcept
    {
        assert(x_node_number >= 0 && x_node_number < this->size_x);
        return Index_proxy_constant{*this, x_node_number};
    }

    constexpr int get_size_x() const noexcept
    {
        return this->size_x;
    }

    constexpr int get_size_y() const noexcept
    {
        return this->size_y;
    }

    constexpr int get_size_z() const noexcept
    {
        return this->size_z;
    }

    constexpr void set_null_value() noexcept
    {
        if (this->container_size != 0) {
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

    constexpr Regular_grid_3d<T>& operator+=(const Regular_grid_3d<T>& add_grid) noexcept
    {
        if (this->size_x == 0 &&
            this->size_y == 0 &&
            this->size_z == 0) {

            this->size_x = add_grid.size_x;
            this->size_y = add_grid.size_y;
            this->size_z = add_grid.size_z;
            this->container_size = add_grid.container_size;
            this->grid_data = add_grid.grid_data;

        } else {

            assert(this->size_x == add_grid.size_x);
            assert(this->size_y == add_grid.size_y);
            assert(this->size_z == add_grid.size_z);
            assert(this->container_size == add_grid.container_size);

            for (int i = 0; i < this->container_size; ++i) {
                this->grid_data[i] += add_grid.grid_data[i];
            }
        }

        return *this;
    }

    constexpr Regular_grid_3d<T>& operator-=(const Regular_grid_3d<T>& sub_grid) noexcept
    {
        if (this->size_x == 0 &&
            this->size_y == 0 &&
            this->size_z == 0) {

            this->size_x = sub_grid.size_x;
            this->size_y = sub_grid.size_y;
            this->size_z = sub_grid.size_z;
            this->container_size = sub_grid.container_size;

            this->set_null_value();

        } else {

            assert(this->size_x == sub_grid.size_x);
            assert(this->size_y == sub_grid.size_y);
            assert(this->size_z == sub_grid.size_z);
            assert(this->container_size == sub_grid.container_size);

        }

        for (int i = 0; i < this->container_size; ++i) {
            this->grid_data[i] -= sub_grid.grid_data[i];
        }

        return *this;
    }

    constexpr Regular_grid_3d<T>& operator*=(const T& multiplier) noexcept
    {
        for (int i = 0; i < this->container_size; ++i) {
            this->grid_data[i] *= multiplier;
        }

        return *this;
    }

    constexpr Regular_grid_3d<T>& operator/=(const T& denom) noexcept
    {
        for (int i = 0; i < this->container_size; ++i) {
            this->grid_data[i] /= denom;
        }

        return *this;
    }

private:

    constexpr int get_converted_index(const int x_node_number,
                                      const int y_node_number,
                                      const int z_node_number) const noexcept
    {
        assert(x_node_number >= 0 && x_node_number < this->size_x);
        assert(y_node_number >= 0 && y_node_number < this->size_y);
        assert(z_node_number >= 0 && z_node_number < this->size_z);

        const int index = x_node_number * this->size_y * this->size_z +
                          y_node_number * this->size_z + z_node_number;

        assert(index >= 0 && index < this->container_size);

        return index;
    }

    int size_x = 0;
    int size_y = 0;
    int size_z = 0;

    int container_size = 0;

    std::vector<T> grid_data{};

};

template<typename T>
class Regular_grid_3d<T>::Index_proxy_variable
{
private:

    class Internal_index_proxy_variable;
    friend class Internal_index_proxy_variable;

public:

    constexpr Index_proxy_variable() noexcept = default;

    constexpr Index_proxy_variable(const Index_proxy_variable&) noexcept = default;
    constexpr Index_proxy_variable(Index_proxy_variable&&) noexcept = default;

    constexpr Index_proxy_variable(Regular_grid_3d<T>& grid,
                                   const int x_node_number) noexcept
        :   grid{grid},
            x_node_number{x_node_number}
    {
        assert(this->grid.container_size != 0);
        assert(this->x_node_number >= 0 && this->x_node_number < this->grid.size_x);
    }

    ~Index_proxy_variable() noexcept = default;

    constexpr Internal_index_proxy_variable operator[](const int y_node_number) noexcept
    {
        assert(y_node_number >= 0 && y_node_number < this->grid.size_y);
        return Internal_index_proxy_variable{*this, y_node_number};
    }

private:

    Regular_grid_3d<T>& grid;
    const int x_node_number = 0;

};

template<typename T>
class Regular_grid_3d<T>::Index_proxy_variable::Internal_index_proxy_variable
{
private:

    using Index_proxy = typename Regular_grid_3d<T>::Index_proxy_variable;

public:

    constexpr Internal_index_proxy_variable() noexcept = default;

    constexpr Internal_index_proxy_variable(const Internal_index_proxy_variable&) noexcept = default;
    constexpr Internal_index_proxy_variable(Internal_index_proxy_variable&&) noexcept = default;

    constexpr Internal_index_proxy_variable(Index_proxy& index_proxy,
                                            const int y_node_number) noexcept
        :   index_proxy{index_proxy},
            y_node_number{y_node_number}
    {
        assert(index_proxy.grid.container_size != 0);
        assert(index_proxy.x_node_number >= 0 && index_proxy.x_node_number < index_proxy.grid.size_x);
        assert(this->y_node_number >= 0 && this->y_node_number < index_proxy.grid.size_y);
    }

    ~Internal_index_proxy_variable() noexcept = default;

    constexpr T& operator[](const int z_node_number) noexcept
    {
        assert(z_node_number >= 0 && z_node_number < index_proxy.grid.size_z);

        const int index = index_proxy.grid.get_converted_index(index_proxy.x_node_number,
                                                               y_node_number, z_node_number);

        return index_proxy.grid.grid_data[index];
    }

private:

    Index_proxy& index_proxy;
    const int y_node_number = 0;

};

template<typename T>
class Regular_grid_3d<T>::Index_proxy_constant
{
private:

    class Internal_index_proxy_constant;
    friend class Internal_index_proxy_constant;

public:

    constexpr Index_proxy_constant() noexcept = default;

    constexpr Index_proxy_constant(const Index_proxy_constant&) noexcept = default;
    constexpr Index_proxy_constant(Index_proxy_constant&&) noexcept = default;

    constexpr Index_proxy_constant(const Regular_grid_3d<T>& grid,
                                   const int x_node_number) noexcept
        :   grid{grid},
            x_node_number{x_node_number}
    {
        assert(this->grid.container_size != 0);
        assert(this->x_node_number >= 0 && this->x_node_number < this->grid.size_x);
    }

    ~Index_proxy_constant() noexcept = default;

    constexpr Internal_index_proxy_constant operator[](const int y_node_number) const noexcept
    {
        assert(y_node_number >= 0 && y_node_number < this->grid.size_y);
        return Internal_index_proxy_constant{*this, y_node_number};
    }

private:

    const Regular_grid_3d<T>& grid;
    const int x_node_number = 0;

};

template<typename T>
class Regular_grid_3d<T>::Index_proxy_constant::Internal_index_proxy_constant
{
private:

    using Index_proxy = typename Regular_grid_3d<T>::Index_proxy_constant;

public:

    constexpr Internal_index_proxy_constant() noexcept = default;

    constexpr Internal_index_proxy_constant(const Internal_index_proxy_constant&) noexcept = default;
    constexpr Internal_index_proxy_constant(Internal_index_proxy_constant&&) noexcept = default;

    constexpr Internal_index_proxy_constant(const Index_proxy& index_proxy,
                                            const int y_node_number) noexcept
        :   index_proxy{index_proxy},
            y_node_number{y_node_number}
    {
        assert(index_proxy.grid.container_size != 0);
        assert(index_proxy.x_node_number >= 0 && index_proxy.x_node_number < index_proxy.grid.size_x);
        assert(this->y_node_number >= 0 && this->y_node_number < index_proxy.grid.size_y);
    }

    ~Internal_index_proxy_constant() noexcept = default;

    constexpr const T& operator[](const int z_node_number) const noexcept
    {
        assert(z_node_number >= 0 && z_node_number < index_proxy.grid.size_z);

        const int index = index_proxy.grid.get_converted_index(index_proxy.x_node_number,
                                                               y_node_number, z_node_number);

        return index_proxy.grid.grid_data[index];
    }

private:

    const Index_proxy& index_proxy;
    const int y_node_number = 0;

};

template<typename T_lhs, typename T_rhs>
constexpr bool operator==(const Regular_grid_3d<T_lhs>& lhs,
                          const Regular_grid_3d<T_rhs>& rhs) noexcept
{
    assert(lhs.get_size_x() == rhs.get_size_x());
    assert(lhs.get_size_y() == rhs.get_size_y());
    assert(lhs.get_size_z() == rhs.get_size_z());

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

    const int size_x = lhs.get_size_x();
    const int size_y = lhs.get_size_y();
    const int size_z = lhs.get_size_z();

    for (int i = 0; i < size_x; ++i) {
        for (int j = 0; j < size_y; ++j) {
            for (int k = 0; k < size_z; ++k) {
                if (!comp(lhs[i][j][k], rhs[i][j][k])) {
                    return false;
                }
            }
        }
    }

    return true;
}

template<typename T_lhs, typename T_rhs>
constexpr Regular_grid_3d<decltype(std::declval<T_lhs>() +
                                   std::declval<T_rhs>())>
    operator+(const Regular_grid_3d<T_lhs>& lhs,
              const Regular_grid_3d<T_rhs>& rhs) noexcept
{
    if (rhs.get_size_x() == 0 &&
        rhs.get_size_y() == 0 &&
        rhs.get_size_z() == 0) {

        assert(lhs.get_size_x() != 0);
        assert(lhs.get_size_y() != 0);
        assert(lhs.get_size_z() != 0);

        return lhs;

    } else {

        assert(lhs.get_size_x() == rhs.get_size_x());
        assert(lhs.get_size_y() == rhs.get_size_y());
        assert(lhs.get_size_z() == rhs.get_size_z());

        const int size_x = lhs.get_size_x();
        const int size_y = lhs.get_size_y();
        const int size_z = lhs.get_size_z();

        using res_type = decltype(std::declval<T_lhs>() +
                                  std::declval<T_rhs>());

        Regular_grid_3d<res_type> res(size_x, size_y, size_z);

        for (int i = 0; i < size_x; ++i) {
            for (int j = 0; j < size_y; ++j) {
                for (int k = 0; k < size_z; ++k) {
                    res[i][j][k] = lhs[i][j][k] + rhs[i][j][k];
                }
            }
        }

        return res;

    }
}

template<typename T_lhs, typename T_rhs>
constexpr Regular_grid_3d<decltype(std::declval<T_lhs>() -
                                   std::declval<T_rhs>())>
    operator-(const Regular_grid_3d<T_lhs>& lhs,
              const Regular_grid_3d<T_rhs>& rhs) noexcept
{
    if (rhs.get_size_x() == 0 &&
        rhs.get_size_y() == 0 &&
        rhs.get_size_z() == 0) {

        assert(lhs.get_size_x() != 0);
        assert(lhs.get_size_y() != 0);
        assert(lhs.get_size_z() != 0);

        return lhs;

    } else {

        assert(lhs.get_size_x() == rhs.get_size_x());
        assert(lhs.get_size_y() == rhs.get_size_y());
        assert(lhs.get_size_z() == rhs.get_size_z());

        const int size_x = lhs.get_size_x();
        const int size_y = lhs.get_size_y();
        const int size_z = lhs.get_size_z();

        using res_type = decltype(std::declval<T_lhs>() -
                                  std::declval<T_rhs>());

        Regular_grid_3d<res_type> res(size_x, size_y, size_z);

        for (int i = 0; i < size_x; ++i) {
            for (int j = 0; j < size_y; ++j) {
                for (int k = 0; k < size_z; ++k) {
                    res[i][j][k] = lhs[i][j][k] - rhs[i][j][k];
                }
            }
        }

        return res;

    }
}

template<typename T_lhs, typename T_rhs>
constexpr Regular_grid_3d<decltype(std::declval<T_lhs>() *
                                   std::declval<T_rhs>())>
    operator*(const Regular_grid_3d<T_lhs>& lhs, const T_rhs& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() *
                              std::declval<T_rhs>());

    const int size_x = lhs.get_size_x();
    const int size_y = lhs.get_size_y();
    const int size_z = lhs.get_size_z();

    Regular_grid_3d<res_type> res(size_x, size_y, size_z);

    for (int i = 0; i < size_x; ++i) {
        for (int j = 0; j < size_y; ++j) {
            for (int k = 0; k < size_z; ++k) {
                res[i][j][k] = lhs[i][j][k] * rhs;
            }
        }
    }

    return res;
}

template<typename T_lhs, typename T_rhs>
constexpr Regular_grid_3d<decltype(std::declval<T_lhs>() *
                                   std::declval<T_rhs>())>
    operator*(const T_lhs& lhs, const Regular_grid_3d<T_rhs>& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() *
                              std::declval<T_rhs>());

    const int size_x = rhs.get_size_x();
    const int size_y = rhs.get_size_y();
    const int size_z = rhs.get_size_z();

    Regular_grid_3d<res_type> res(size_x, size_y, size_z);

    for (int i = 0; i < size_x; ++i) {
        for (int j = 0; j < size_y; ++j) {
            for (int k = 0; k < size_z; ++k) {
                res[i][j][k] = lhs * rhs[i][j][k];
            }
        }
    }

    return res;
}

template<typename T_lhs, typename T_rhs>
constexpr Regular_grid_3d<decltype(std::declval<T_lhs>() /
                                   std::declval<T_rhs>())>
    operator/(const Regular_grid_3d<T_lhs>& lhs, const T_rhs& rhs) noexcept
{
    using res_type = decltype(std::declval<T_lhs>() /
                              std::declval<T_rhs>());

    const int size_x = lhs.get_size_x();
    const int size_y = lhs.get_size_y();
    const int size_z = lhs.get_size_z();

    Regular_grid_3d<res_type> res(size_x, size_y, size_z);

    for (int i = 0; i < size_x; ++i) {
        for (int j = 0; j < size_y; ++j) {
            for (int k = 0; k < size_z; ++k) {
                res[i][j][k] = lhs[i][j][k] / rhs;
            }
        }
    }

    return res;
}

template<typename T>
constexpr auto max_abs(const Regular_grid_3d<T>& m) noexcept
{
    constexpr bool is_simple_standard_math_type =
                            std::is_arithmetic_v<std::decay_t<T>> ||
                            std::is_same_v<std::decay_t<T>, std::complex<float>> ||
                            std::is_same_v<std::decay_t<T>, std::complex<double>> ||
                            std::is_same_v<std::decay_t<T>, std::complex<long double>>;

    if constexpr (is_simple_standard_math_type) {

        const int size_x = m.get_size_x();
        const int size_y = m.get_size_y();
        const int size_z = m.get_size_z();

        using max_value_t = decltype(std::abs(std::declval<decltype(m[0][0][0])>()));

        max_value_t max_value = max_value_t(0);

        for (int i = 0; i < size_x; ++i) {
            for (int j = 0; j < size_y; ++j) {
                for (int k = 0; k < size_z; ++k) {
                    max_value = std::max(max_value, std::abs(m[i][j][k]));
                }
            }
        }

        return max_value;

    } else {

        const int size_x = m.get_size_x();
        const int size_y = m.get_size_y();
        const int size_z = m.get_size_z();

        using max_value_t = decltype(max_abs(std::declval<decltype(m[0][0][0])>()));

        max_value_t max_value = max_value_t();

        for (int i = 0; i < size_x; ++i) {
            for (int j = 0; j < size_y; ++j) {
                for (int k = 0; k < size_z; ++k) {
                    max_value = std::max(max_value, max_abs(m[i][j][k]));
                }
            }
        }

        return max_value;

    }
}

}

#endif // REGULAR_GRID_3D_HXX_INCLUDED
