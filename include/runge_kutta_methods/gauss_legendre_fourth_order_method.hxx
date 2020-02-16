#ifndef GAUSS_LEGENDRE_FOURTH_ORDER_METHOD_HXX_INCLUDED
#define GAUSS_LEGENDRE_FOURTH_ORDER_METHOD_HXX_INCLUDED

#include "runge_kutta_method.hxx"
#include "runge_kutta_method_impl.hxx"

namespace Math_structures
{

template<
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class Gauss_Legendre_fourth_order_method
    :   public Runge_Kutta_method_impl<2, T_arg_u, T_arg_t>
{
public:

    constexpr Gauss_Legendre_fourth_order_method() noexcept
        :   Runge_Kutta_method_impl<2, T_arg_u, T_arg_t>{
                    {4}, {
                        {{
                            {{
                                0.25,
                                0.25 - std::sqrt(3.0)/6.0
                            }},
                            {{
                                0.25 + std::sqrt(3.0)/6.0,
                                0.25
                            }}
                        }},
                        {{
                            0.5 - std::sqrt(3.0)/6.0,
                            0.5 + std::sqrt(3.0)/6.0
                        }},
                        {{0.5, 0.5}}
                    }}
    { }

    constexpr Gauss_Legendre_fourth_order_method
        (const Gauss_Legendre_fourth_order_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Gauss_Legendre_fourth_order_method
        (Gauss_Legendre_fourth_order_method<T_arg_u, T_arg_t>&&) noexcept = default;

    virtual ~Gauss_Legendre_fourth_order_method() noexcept = default;

    constexpr Gauss_Legendre_fourth_order_method<T_arg_u, T_arg_t>&
        operator=(const Gauss_Legendre_fourth_order_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Gauss_Legendre_fourth_order_method<T_arg_u, T_arg_t>&
        operator=(Gauss_Legendre_fourth_order_method<T_arg_u, T_arg_t>&&) noexcept = default;

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

#endif // GAUSS_LEGENDRE_FOURTH_ORDER_METHOD_HXX_INCLUDED
