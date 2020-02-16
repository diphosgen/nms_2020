#ifndef EXPLICIT_CASH_KARP_FIFTH_ORDER_METHOD_HXX_INCLUDED
#define EXPLICIT_CASH_KARP_FIFTH_ORDER_METHOD_HXX_INCLUDED

#include "runge_kutta_method.hxx"
#include "runge_kutta_method_expl.hxx"

namespace Math_structures
{

template<
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class Explicit_Cash_Karp_fifth_order_method
    :   public Runge_Kutta_method_expl<6, T_arg_u, T_arg_t>
{
public:

    constexpr Explicit_Cash_Karp_fifth_order_method() noexcept
        :   Runge_Kutta_method_expl<6, T_arg_u, T_arg_t>{
                    {5}, {
                        {{
                            {{0.0,              0.0,            0.0,            0.0,                0.0,            0.0}},
                            {{1.0/5.0,          0.0,            0.0,            0.0,                0.0,            0.0}},
                            {{3.0/40.0,         9.0/40.0,       0.0,            0.0,                0.0,            0.0}},
                            {{3.0/10.0,         -9.0/10.0,      6.0/5.0,        0.0,                0.0,            0.0}},
                            {{-11.0/54.0,       5.0/2.0,        -70.0/27.0,     35.0/27.0,          0.0,            0.0}},
                            {{1631.0/55296.0,   175.0/512.0,    575.0/13824.0,  44275.0/110592.0,   253.0/4096.0,   0.0}}
                        }},
                        {{0.0,          1.0/5.0,    3.0/10.0,       3.0/5.0,            1.0,        7.0/8.0}},
                        {{37.0/378.0,   0.0,        250.0/621.0,    125.0/594.0,        0.0,        512.0/1771.0}}
                    }}
    { }

    constexpr Explicit_Cash_Karp_fifth_order_method
        (const Explicit_Cash_Karp_fifth_order_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Explicit_Cash_Karp_fifth_order_method
        (Explicit_Cash_Karp_fifth_order_method<T_arg_u, T_arg_t>&&) noexcept = default;

    virtual ~Explicit_Cash_Karp_fifth_order_method() noexcept = default;

    constexpr Explicit_Cash_Karp_fifth_order_method<T_arg_u, T_arg_t>&
        operator=(const Explicit_Cash_Karp_fifth_order_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Explicit_Cash_Karp_fifth_order_method<T_arg_u, T_arg_t>&
        operator=(Explicit_Cash_Karp_fifth_order_method<T_arg_u, T_arg_t>&&) noexcept = default;

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

#endif // EXPLICIT_CASH_KARP_FIFTH_ORDER_METHOD_HXX_INCLUDED
