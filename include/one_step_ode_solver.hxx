#ifndef ONE_STEP_ODE_SOLVER_HXX_INCLUDED
#define ONE_STEP_ODE_SOLVER_HXX_INCLUDED

#include "basic_libs.hxx"
#include "basic_math_functions.hxx"

namespace Math_structures
{

template<
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class One_step_ODE_solver
{
public:

    constexpr One_step_ODE_solver() noexcept = default;

    constexpr One_step_ODE_solver(const One_step_ODE_solver<T_arg_u, T_arg_t>&) noexcept = default;
    constexpr One_step_ODE_solver(One_step_ODE_solver<T_arg_u, T_arg_t>&&) noexcept = default;

    constexpr One_step_ODE_solver(const int accuracy_order)
        :   accuracy_order{accuracy_order}
    { }

    virtual ~One_step_ODE_solver() noexcept = default;

    constexpr One_step_ODE_solver<T_arg_u, T_arg_t>&
        operator=(const One_step_ODE_solver<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr One_step_ODE_solver<T_arg_u, T_arg_t>&
        operator=(One_step_ODE_solver<T_arg_u, T_arg_t>&&) noexcept = default;

    virtual int get_accuracy_order() const noexcept
    {
        return this->accuracy_order;
    }

    virtual T_arg_u get_solution(const std::function<T_arg_u(T_arg_u, T_arg_t)>& f_rhs,
                                 const T_arg_u& u_prev, const T_arg_t& t, const T_arg_t& dt) = 0;

private:

    const int accuracy_order = 0;

};

class One_step_ODE_solver_does_not_obtained_solution
    :   public std::runtime_error
{
public:

    One_step_ODE_solver_does_not_obtained_solution() noexcept
        :   std::runtime_error{"One-step ODE solver does not obtained a solution."}
    { }

    One_step_ODE_solver_does_not_obtained_solution
        (const One_step_ODE_solver_does_not_obtained_solution&) noexcept = default;

    One_step_ODE_solver_does_not_obtained_solution
        (One_step_ODE_solver_does_not_obtained_solution&&) noexcept = default;

    virtual ~One_step_ODE_solver_does_not_obtained_solution() noexcept = default;

    One_step_ODE_solver_does_not_obtained_solution&
        operator=(const One_step_ODE_solver_does_not_obtained_solution&) noexcept = default;

    One_step_ODE_solver_does_not_obtained_solution&
        operator=(One_step_ODE_solver_does_not_obtained_solution&&) noexcept = default;

    virtual const char* what() const noexcept
    {
        return "A solution to the system of differential equations by "
               "this one-step ODE solver could not be obtained: try a otherwise "
               "method or a another parameters of this method.";
    }
};

}

#endif // ONE_STEP_ODE_SOLVER_HXX_INCLUDED
