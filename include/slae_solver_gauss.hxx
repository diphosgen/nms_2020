#ifndef SLAE_SOLVER_GAUSS_HXX_INCLUDED
#define SLAE_SOLVER_GAUSS_HXX_INCLUDED

#include "slae_solver.hxx"

namespace Math_structures
{

template<
    typename Matrix_type = Matrix_unfixed<double>,
    typename Vector_type = Vector_unfixed<double>
    >
class SLAE_solver_Gauss
    :   public SLAE_solver<Matrix_type, Vector_type>
{
public:

    constexpr SLAE_solver_Gauss() noexcept = default;

    constexpr SLAE_solver_Gauss(const SLAE_solver_Gauss<Matrix_type, Vector_type>&) noexcept = default;
    constexpr SLAE_solver_Gauss(SLAE_solver_Gauss<Matrix_type, Vector_type>&&) noexcept = default;

    constexpr SLAE_solver_Gauss<Matrix_type, Vector_type>&
        operator=(const SLAE_solver_Gauss<Matrix_type, Vector_type>&) noexcept = default;

    constexpr SLAE_solver_Gauss<Matrix_type, Vector_type>&
        operator=(SLAE_solver_Gauss<Matrix_type, Vector_type>&&) noexcept = default;

    virtual ~SLAE_solver_Gauss() noexcept = default;

    virtual Vector_type get_solution(const Matrix_type& matrix_slae,
                                     const Vector_type& vector_slae) const noexcept override
    {
        assert(matrix_slae.get_cols_amount() == matrix_slae.get_rows_amount());
        assert(matrix_slae.get_cols_amount() == vector_slae.size());

        Matrix_type matrix_system = matrix_slae;
        Vector_type vector_system = vector_slae;

        const int system_size = vector_slae.size();

        using vector_elem_t = std::decay_t<typename Vector_type::value_type>;
        using matrix_elem_t = std::decay_t<typename Matrix_type::value_type>;

        using vector_elem_ptr_t = std::add_pointer_t<vector_elem_t>;

        std::vector<vector_elem_ptr_t> vector_system_ptr(system_size);
        std::vector<int> index_ptr(system_size);

        for (int row = 0; row < system_size; ++row) {
            index_ptr[row] = row;
            vector_system_ptr[row] = &vector_system[row];
        }

        for (int row = 0; row < system_size - 1; ++row) {

            using max_value_t = decltype(max_abs(std::declval<decltype(matrix_system[0][0])>()));

            max_value_t max_value{};

            constexpr bool is_simple_standard_math_type =
                                std::is_arithmetic_v<std::decay_t<matrix_elem_t>> ||
                                std::is_same_v<std::decay_t<matrix_elem_t>, std::complex<float>> ||
                                std::is_same_v<std::decay_t<matrix_elem_t>, std::complex<double>> ||
                                std::is_same_v<std::decay_t<matrix_elem_t>, std::complex<long double>>;

            if constexpr (is_simple_standard_math_type) {
                max_value = max_value_t(0);
            } else {
                max_value.set_null_value();
            }

            int row_max = row;
            int col_max = row;

            for (int row_local = row; row_local < system_size; ++row_local) {
                for (int col_local = row; col_local < system_size; ++col_local) {

                    const max_value_t abs_of_element = max_abs(matrix_system[row_local][col_local]);

                    if (abs_of_element > max_value) {
                        max_value = abs_of_element;
                        row_max = row_local;
                        col_max = col_local;
                    }
                }
            }

            using std::swap;

            if (row != row_max) {
                for (int col_local = row; col_local < system_size; ++col_local) {
                    swap(matrix_system[row][col_local],matrix_system[row_max][col_local]);
                }

                swap(vector_system_ptr[row],vector_system_ptr[row_max]);
            }

            if (row != col_max) {
                for (int row_local = 0; row_local < system_size; ++row_local) {
                    swap(matrix_system[row_local][row],matrix_system[row_local][col_max]);
                }

                swap(index_ptr[row],index_ptr[col_max]);
            }

            for (int row_local = row + 1; row_local < system_size; ++row_local) {

                const matrix_elem_t mult_temp = matrix_system[row_local][row] / matrix_system[row][row];

                for (int col_local = row; col_local < system_size; ++col_local) {
                    matrix_system[row_local][col_local] -= mult_temp * matrix_system[row][col_local];
                }

                (*vector_system_ptr[row_local]) -= mult_temp * (*vector_system_ptr[row]);
            }
        }

        Vector_type vector_res_temp(system_size);

        for (int step = 0; step < system_size; ++step) {

            const int row = system_size - step - 1;

            vector_elem_t sum{};

            constexpr bool is_simple_standard_math_type =
                                std::is_arithmetic_v<std::decay_t<vector_elem_t>> ||
                                std::is_same_v<std::decay_t<vector_elem_t>, std::complex<float>> ||
                                std::is_same_v<std::decay_t<vector_elem_t>, std::complex<double>> ||
                                std::is_same_v<std::decay_t<vector_elem_t>, std::complex<long double>>;

            if constexpr (is_simple_standard_math_type) {
                sum = vector_elem_t(0);
            } else {
                sum.set_null_value();
            }

            for (int col = row + 1; col < system_size; ++col) {
                sum += matrix_system[row][col] * vector_res_temp[col];
            }

            vector_res_temp[row] = ((*vector_system_ptr[row]) - sum) / matrix_system[row][row];
        }

        Vector_type vector_res(system_size);

        for (int row = 0; row < system_size; ++row) {
            vector_res[index_ptr[row]] = vector_res_temp[row];
        }

        return vector_res;
    }

};

}

#endif // SLAE_SOLVER_GAUSS_HXX_INCLUDED
