#ifndef EXPLICIT_RALSTON_FOURTH_ORDER_METHOD_HXX_INCLUDED
#define EXPLICIT_RALSTON_FOURTH_ORDER_METHOD_HXX_INCLUDED

#include "runge_kutta_method.hxx"
#include "runge_kutta_method_expl.hxx"

namespace Math_structures
{

template<
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class Explicit_Ralston_fourth_order_method
    :   public Runge_Kutta_method_expl<4, T_arg_u, T_arg_t>
{
public:

    constexpr Explicit_Ralston_fourth_order_method() noexcept
        :   Runge_Kutta_method_expl<4, T_arg_u, T_arg_t>{
                    {4}, {
                        {{
                            {{0.0,        0.0,         0.0,        0.0}},
                            {{0.4,        0.0,         0.0,        0.0}},
                            {{0.29697761, 0.15875964,  0.0,        0.0}},
                            {{0.21810040, -3.05096516, 3.83286476, 0.0}}
                        }},
                        {{0.0,        0.4,         0.45573725, 1.0}},
                        {{0.17476028, -0.55148066, 1.20553560, 0.17118478}}
                    }}
    { }

    constexpr Explicit_Ralston_fourth_order_method
        (const Explicit_Ralston_fourth_order_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Explicit_Ralston_fourth_order_method
        (Explicit_Ralston_fourth_order_method<T_arg_u, T_arg_t>&&) noexcept = default;

    virtual ~Explicit_Ralston_fourth_order_method() noexcept = default;

    constexpr Explicit_Ralston_fourth_order_method<T_arg_u, T_arg_t>&
        operator=(const Explicit_Ralston_fourth_order_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Explicit_Ralston_fourth_order_method<T_arg_u, T_arg_t>&
        operator=(Explicit_Ralston_fourth_order_method<T_arg_u, T_arg_t>&&) noexcept = default;

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

#endif // EXPLICIT_RALSTON_FOURTH_ORDER_METHOD_HXX_INCLUDED
