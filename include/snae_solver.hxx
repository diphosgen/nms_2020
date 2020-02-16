#ifndef SNAE_SOLVER_HXX_INCLUDED
#define SNAE_SOLVER_HXX_INCLUDED

#include "basic_libs.hxx"
#include "basic_math_functions.hxx"

namespace Math_structures
{

template<typename T_args>
class SNAE_solver
{
public:

    constexpr SNAE_solver() noexcept = default;

    constexpr SNAE_solver(const SNAE_solver<T_args>&) noexcept = default;
    constexpr SNAE_solver(SNAE_solver<T_args>&&) noexcept = default;

    virtual ~SNAE_solver() noexcept = default;

    constexpr SNAE_solver<T_args>& operator=(const SNAE_solver<T_args>&) noexcept = default;
    constexpr SNAE_solver<T_args>& operator=(SNAE_solver<T_args>&&) noexcept = default;

    virtual T_args get_solution(const std::function<T_args(T_args)>& solve_function,
                                const T_args& initial_value) const = 0;

};

class SNAE_solution_not_obtained
    :   public std::runtime_error
{
public:

    SNAE_solution_not_obtained() noexcept
        :   std::runtime_error{"SNAE solution not obtained."}
    { }

    SNAE_solution_not_obtained(const SNAE_solution_not_obtained&) noexcept = default;
    SNAE_solution_not_obtained(SNAE_solution_not_obtained&&) noexcept = default;

    virtual ~SNAE_solution_not_obtained() noexcept = default;

    SNAE_solution_not_obtained& operator=(const SNAE_solution_not_obtained&) noexcept = default;
    SNAE_solution_not_obtained& operator=(SNAE_solution_not_obtained&&) noexcept = default;

    virtual const char* what() const noexcept
    {
        return "A solution to the system of nonlinear algebraic equations "
               "could not be obtained: try a otherwise method or a another "
               "parameters of this method.";
    }
};

}

#endif // SNAE_SOLVER_HXX_INCLUDED
