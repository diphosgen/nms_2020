#ifndef EXPLICIT_RUNGE_KUTTA_FEHLBERG_FIFTH_ORDER_METHOD_HXX_INCLUDED
#define EXPLICIT_RUNGE_KUTTA_FEHLBERG_FIFTH_ORDER_METHOD_HXX_INCLUDED

#include "runge_kutta_method.hxx"
#include "runge_kutta_method_expl.hxx"

namespace Math_structures
{

template<
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class Explicit_Runge_Kutta_Fehlberg_fifth_order_method
    :   public Runge_Kutta_method_expl<6, T_arg_u, T_arg_t>
{
public:

    constexpr Explicit_Runge_Kutta_Fehlberg_fifth_order_method() noexcept
        :   Runge_Kutta_method_expl<6, T_arg_u, T_arg_t>{
                    {5}, {
                        {{
                            {{0.0,              0.0,            0.0,            0.0,            0.0,            0.0}},
                            {{1.0/4.0,          0.0,            0.0,            0.0,            0.0,            0.0}},
                            {{3.0/32.0,         9.0/32.0,       0.0,            0.0,            0.0,            0.0}},
                            {{1932.0/2197.0,    -7200.0/2197.0, 7296.0/2197.0,  0.0,            0.0,            0.0}},
                            {{439.0/216.0,      -8.0,           3680.0/513.0,   -845.0/4104.0,  0.0,            0.0}},
                            {{-8.0/27.0,        2.0,            -3544.0/2565.0, 1859.0/4104.0,  -11.0/40.0,     0.0}}
                        }},
                        {{0.0,          1.0/4.0,    3.0/8.0,        12.0/13.0,          1.0,        1.0/2.0}},
                        {{16.0/135.0,   0.0,        6656.0/12825.0, 28561.0/56430.0,    -9.0/50.0,  2.0/55.0}}
                    }}
    { }

    constexpr Explicit_Runge_Kutta_Fehlberg_fifth_order_method
        (const Explicit_Runge_Kutta_Fehlberg_fifth_order_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Explicit_Runge_Kutta_Fehlberg_fifth_order_method
        (Explicit_Runge_Kutta_Fehlberg_fifth_order_method<T_arg_u, T_arg_t>&&) noexcept = default;

    virtual ~Explicit_Runge_Kutta_Fehlberg_fifth_order_method() noexcept = default;

    constexpr Explicit_Runge_Kutta_Fehlberg_fifth_order_method<T_arg_u, T_arg_t>&
        operator=(const Explicit_Runge_Kutta_Fehlberg_fifth_order_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Explicit_Runge_Kutta_Fehlberg_fifth_order_method<T_arg_u, T_arg_t>&
        operator=(Explicit_Runge_Kutta_Fehlberg_fifth_order_method<T_arg_u, T_arg_t>&&) noexcept = default;

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

#endif // EXPLICIT_RUNGE_KUTTA_FEHLBERG_FIFTH_ORDER_METHOD_HXX_INCLUDED
