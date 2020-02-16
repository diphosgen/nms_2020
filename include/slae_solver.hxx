#ifndef SLAE_SOLVER_HXX_INCLUDED
#define SLAE_SOLVER_HXX_INCLUDED

#include "basic_libs.hxx"
#include "basic_math_functions.hxx"
#include "vector_unfixed.hxx"
#include "matrix_unfixed.hxx"

namespace Math_structures
{

template<
    typename Matrix_type = Matrix_unfixed<double>,
    typename Vector_type = Vector_unfixed<double>
    >
class SLAE_solver
{
public:

    constexpr SLAE_solver() noexcept = default;

    constexpr SLAE_solver(const SLAE_solver<Matrix_type, Vector_type>&) noexcept = default;
    constexpr SLAE_solver(SLAE_solver<Matrix_type, Vector_type>&&) noexcept = default;

    virtual ~SLAE_solver() noexcept = default;

    constexpr SLAE_solver<Matrix_type, Vector_type>&
        operator=(const SLAE_solver<Matrix_type, Vector_type>&) noexcept = default;

    constexpr SLAE_solver<Matrix_type, Vector_type>&
        operator=(SLAE_solver<Matrix_type, Vector_type>&&) noexcept = default;

    virtual Vector_type get_solution(const Matrix_type& matrix_slae,
                                     const Vector_type& vector_slae) const = 0;

};

}

#endif // SLAE_SOLVER_HXX_INCLUDED
