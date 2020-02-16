#ifndef RUNGE_KUTTA_METHOD_EXPL_HXX_INCLUDED
#define RUNGE_KUTTA_METHOD_EXPL_HXX_INCLUDED

#include "runge_kutta_method.hxx"

namespace Math_structures
{

template<
    int stages_amount,
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class Runge_Kutta_method_expl
    :   public Runge_Kutta_method<stages_amount, T_arg_u, T_arg_t>
{
private:

    virtual void make_increments_calculation(const std::function<T_arg_u(T_arg_u, T_arg_t)>& f_rhs,
                                             const T_arg_u& u_prev, const T_arg_t& t, const T_arg_t& dt) override
    {
        for (int stage = 0; stage < stages_amount; ++stage) {

            T_arg_u temp_value = u_prev;

            for (int i = 0; i < stage; ++i) {
                const T_arg_u temp_increment_value = this->butcher_table.get_element_a_table(stage, i) *
                                                     this->rkm_increments_k[i];
                temp_value += dt * temp_increment_value;
            }

            const T_arg_t temp_time = t + this->butcher_table.get_element_c_table(stage) * dt;
            this->rkm_increments_k[stage] = f_rhs(temp_value, temp_time);
        }
    }

public:

    constexpr Runge_Kutta_method_expl() noexcept = delete;

    constexpr Runge_Kutta_method_expl
        (const Runge_Kutta_method_expl<stages_amount, T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Runge_Kutta_method_expl
        (Runge_Kutta_method_expl<stages_amount, T_arg_u, T_arg_t>&&) noexcept = default;

    constexpr Runge_Kutta_method_expl(const int accuracy_order,
                                      const Butcher_table<stages_amount, T_arg_t>& butcher_table)
        :   Runge_Kutta_method<stages_amount, T_arg_u, T_arg_t>{accuracy_order, butcher_table}
    { }

    constexpr Runge_Kutta_method_expl<stages_amount, T_arg_u, T_arg_t>&
        operator=(const Runge_Kutta_method_expl<stages_amount, T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Runge_Kutta_method_expl<stages_amount, T_arg_u, T_arg_t>&
        operator=(Runge_Kutta_method_expl<stages_amount, T_arg_u, T_arg_t>&&) noexcept = default;

    virtual ~Runge_Kutta_method_expl() noexcept = default;

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

}

#endif // RUNGE_KUTTA_METHOD_EXPL_HXX_INCLUDED
