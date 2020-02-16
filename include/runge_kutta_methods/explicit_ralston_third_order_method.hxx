#ifndef EXPLICIT_RALSTON_THIRD_ORDER_METHOD_HXX_INCLUDED
#define EXPLICIT_RALSTON_THIRD_ORDER_METHOD_HXX_INCLUDED

#include "runge_kutta_method.hxx"
#include "runge_kutta_method_expl.hxx"

namespace Math_structures
{

template<
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class Explicit_Ralston_third_order_method
    :   public Runge_Kutta_method_expl<3, T_arg_u, T_arg_t>
{
public:

    constexpr Explicit_Ralston_third_order_method() noexcept
        :   Runge_Kutta_method_expl<3, T_arg_u, T_arg_t>{
                    {3}, {
                        {{
                            {{0.0, 0.0,  0.0}},
                            {{0.5, 0.0,  0.0}},
                            {{0.0, 0.75, 0.0}}
                        }},
                        {{0.0,     0.5,     0.75}},
                        {{2.0/9.0, 1.0/3.0, 4.0/9.0}}
                    }}
    { }

    constexpr Explicit_Ralston_third_order_method
        (const Explicit_Ralston_third_order_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Explicit_Ralston_third_order_method
        (Explicit_Ralston_third_order_method<T_arg_u, T_arg_t>&&) noexcept = default;

    virtual ~Explicit_Ralston_third_order_method() noexcept = default;

    constexpr Explicit_Ralston_third_order_method<T_arg_u, T_arg_t>&
        operator=(const Explicit_Ralston_third_order_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Explicit_Ralston_third_order_method<T_arg_u, T_arg_t>&
        operator=(Explicit_Ralston_third_order_method<T_arg_u, T_arg_t>&&) noexcept = default;

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

#endif // EXPLICIT_RALSTON_THIRD_ORDER_METHOD_HXX_INCLUDED
