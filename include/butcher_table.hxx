#ifndef BUTCHER_TABLE_HXX_INCLUDED
#define BUTCHER_TABLE_HXX_INCLUDED

#include "basic_libs.hxx"

namespace Math_structures
{

template<
    int stages_amount,
    typename T_elems = double
    >
class Butcher_table
{
private:

    template<typename T>
    using Table_1d = std::array<T, stages_amount>;

    template<typename T>
    using Table_2d = std::array<Table_1d<T>, stages_amount>;

public:

    static_assert(stages_amount > 0);

    constexpr Butcher_table(const Table_2d<T_elems>& a_table,
                            const Table_1d<T_elems>& c_table,
                            const Table_1d<T_elems>& b_table) noexcept
        :   a_table{a_table},
            c_table{c_table},
            b_table{b_table}
    {
        assert(this->make_checked_b_table());
        assert(this->make_checked_ac_tables());
    }

    constexpr Butcher_table(Table_2d<T_elems>&& a_table,
                            Table_1d<T_elems>&& c_table,
                            Table_1d<T_elems>&& b_table) noexcept
        :   a_table{std::move(a_table)},
            c_table{std::move(c_table)},
            b_table{std::move(b_table)}
    {
        assert(this->make_checked_b_table());
        assert(this->make_checked_ac_tables());
    }

    template<
        typename T_elems_copy,
        typename = std::enable_if_t<
                        !std::is_same_v<
                                std::decay_t<T_elems_copy>,
                                std::decay_t<T_elems>
                            > &&
                        std::is_convertible_v<
                                std::decay_t<T_elems_copy>,
                                std::decay_t<T_elems>
                            >
                        >
        >
    constexpr Butcher_table(const Table_2d<T_elems_copy>& a_table,
                            const Table_1d<T_elems_copy>& c_table,
                            const Table_1d<T_elems_copy>& b_table) noexcept
    {
        for (int i = 0; i < stages_amount; ++i) {
            for (int j = 0; j < stages_amount; ++j) {
                this->a_table[i][j] = static_cast<T_elems>(a_table[i][j]);
            }
            this->b_table[i] = static_cast<T_elems>(b_table[i]);
            this->c_table[i] = static_cast<T_elems>(c_table[i]);
        }

        assert(this->make_checked_b_table());
        assert(this->make_checked_ac_tables());
    }

    ~Butcher_table() noexcept = default;

    constexpr T_elems get_element_a_table(const int i, const int j) const noexcept
    {
        assert(i >= 0 && i < stages_amount);
        assert(j >= 0 && j < stages_amount);

        return this->a_table[i][j];
    }

    constexpr T_elems get_element_c_table(const int i) const noexcept
    {
        assert(i >= 0 && i < stages_amount);

        return this->c_table[i];
    }

    constexpr T_elems get_element_b_table(const int i) const noexcept
    {
        assert(i >= 0 && i < stages_amount);

        return this->b_table[i];
    }

private:

    const Table_2d<T_elems> a_table{};
    const Table_1d<T_elems> c_table{};
    const Table_1d<T_elems> b_table{};

    static constexpr const bool is_floating_point_math_type =
            std::is_floating_point_v<std::decay_t<T_elems>> ||
            std::is_same_v<std::decay_t<T_elems>, std::complex<float>> ||
            std::is_same_v<std::decay_t<T_elems>, std::complex<double>> ||
            std::is_same_v<std::decay_t<T_elems>, std::complex<long double>>;

    constexpr bool make_checked_b_table() const noexcept
    {
        T_elems sum = 0.0;

        for (int i = 0; i < stages_amount; ++i) {
            sum += this->b_table[i];
        }

        const T_elems exact_sum = T_elems(1.0);

        if constexpr (is_floating_point_math_type) {
            return compare_floating_numbers(sum, exact_sum);
        } else {
            return sum == exact_sum;
        }
    }

    constexpr bool make_checked_ac_tables() const noexcept
    {
        for (int i = 0; i < stages_amount; ++i) {

            T_elems sum = 0.0;

            for (int j = 0; j < stages_amount; ++j) {
                sum += this->a_table[i][j];
            }

            const T_elems exact_sum = this->c_table[i];

            if constexpr (is_floating_point_math_type) {
                const bool res = compare_floating_numbers(sum, exact_sum);
                if (!res) {
                    return false;
                }
            } else {
                const bool res = {sum == exact_sum};
                if (!res) {
                    return false;
                }
            }

        }

        return true;
    }
};

}

#endif // BUTCHER_TABLE_HXX_INCLUDED
