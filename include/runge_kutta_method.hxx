#ifndef RUNGE_KUTTA_METHOD_HXX_INCLUDED
#define RUNGE_KUTTA_METHOD_HXX_INCLUDED

#include "vector_fixed.hxx"
#include "one_step_ode_solver.hxx"
#include "butcher_table.hxx"

namespace Math_structures
{

template<
    int stages_amount,
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class Runge_Kutta_method
    :   public One_step_ODE_solver<T_arg_u, T_arg_t>
{
protected:

    Vector_fixed<T_arg_u, stages_amount> rkm_increments_k =
        Vector_fixed<T_arg_u, stages_amount>{};

    const Butcher_table<stages_amount, T_arg_t> butcher_table{};

    virtual void make_increments_calculation(const std::function<T_arg_u(T_arg_u, T_arg_t)>& f_rhs,
                                             const T_arg_u& u_prev, const T_arg_t& t, const T_arg_t& dt) = 0;

    virtual T_arg_u get_increment_step(const T_arg_u& u_prev, const T_arg_t& dt) const noexcept
    {
        T_arg_u u_next = u_prev;

        for (int i = 0; i < stages_amount; ++i) {
            const T_arg_u weighted_increment = this->butcher_table.get_element_b_table(i) *
                                               this->rkm_increments_k[i];
            u_next += dt * weighted_increment;
        }

        return u_next;
    }

    virtual T_arg_u get_u_next(const std::function<T_arg_u(T_arg_u, T_arg_t)>& f_rhs,
                               const T_arg_u& u_prev, const T_arg_t& t, const T_arg_t& dt)
    {
        this->make_increments_calculation(f_rhs, u_prev, t, dt);
        const T_arg_u u_next = this->get_increment_step(u_prev, dt);

        return u_next;
    }

public:

    constexpr Runge_Kutta_method() noexcept = delete;

    constexpr Runge_Kutta_method(const Runge_Kutta_method<stages_amount, T_arg_u, T_arg_t>&) noexcept = default;
    constexpr Runge_Kutta_method(Runge_Kutta_method<stages_amount, T_arg_u, T_arg_t>&&) noexcept = default;

    constexpr Runge_Kutta_method(const int accuracy_order,
                                 const Butcher_table<stages_amount, T_arg_t>& butcher_table) noexcept
        :   One_step_ODE_solver<T_arg_u, T_arg_t>{accuracy_order},
            butcher_table{butcher_table}
    { }

    virtual ~Runge_Kutta_method() noexcept = default;

    constexpr Runge_Kutta_method<stages_amount, T_arg_u, T_arg_t>&
        operator=(const Runge_Kutta_method<stages_amount, T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Runge_Kutta_method<stages_amount, T_arg_u, T_arg_t>&
        operator=(Runge_Kutta_method<stages_amount, T_arg_u, T_arg_t>&&) noexcept = default;

    virtual int get_accuracy_order() const noexcept override
    {
        return One_step_ODE_solver<T_arg_u, T_arg_t>
               ::get_accuracy_order();
    }

    virtual T_arg_u get_solution(const std::function<T_arg_u(T_arg_u, T_arg_t)>& f_rhs,
                                 const T_arg_u& u_prev, const T_arg_t& t, const T_arg_t& dt) override
    {
        return this->get_u_next(f_rhs, u_prev, t, dt);
    }
};

class Runge_Kutta_method_solution_not_obtained
    :   public One_step_ODE_solver_does_not_obtained_solution
{
public:

    Runge_Kutta_method_solution_not_obtained() noexcept
        :   One_step_ODE_solver_does_not_obtained_solution{}
    { }

    Runge_Kutta_method_solution_not_obtained(const Runge_Kutta_method_solution_not_obtained&) noexcept = default;
    Runge_Kutta_method_solution_not_obtained(Runge_Kutta_method_solution_not_obtained&&) noexcept = default;

    virtual ~Runge_Kutta_method_solution_not_obtained() noexcept = default;

    Runge_Kutta_method_solution_not_obtained&
        operator=(const Runge_Kutta_method_solution_not_obtained&) noexcept = default;

    Runge_Kutta_method_solution_not_obtained&
        operator=(Runge_Kutta_method_solution_not_obtained&&) noexcept = default;

    virtual const char* what() const noexcept
    {
        return "A solution to the system of differential equations by "
               "Runge Kutta method could not be obtained: try a otherwise "
               "method or a another parameters of this method.";
    }
};

}

#endif // RUNGE_KUTTA_METHOD_HXX_INCLUDED
