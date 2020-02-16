#ifndef RADAU_IIA_FIFTH_ORDER_METHOD_HXX_INCLUDED
#define RADAU_IIA_FIFTH_ORDER_METHOD_HXX_INCLUDED

#include "runge_kutta_method.hxx"
#include "runge_kutta_method_impl.hxx"

namespace Math_structures
{

template<
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class Radau_IIA_fifth_order_method
    :   public Runge_Kutta_method_impl<3, T_arg_u, T_arg_t>
{
public:

    constexpr Radau_IIA_fifth_order_method() noexcept
        :   Runge_Kutta_method_impl<3, T_arg_u, T_arg_t>{
                    {5}, {
                        {{
                            {{
                                11.0/45.0 - 7.0 * std::sqrt(6.0)/360.0,
                                37.0/225.0 - 169.0 * std::sqrt(6.0)/1800.0,
                                -2.0/225.0 + std::sqrt(6.0)/75.0
                            }},
                            {{
                                37.0/225.0 + 169.0 * std::sqrt(6.0)/1800.0,
                                11.0/45.0 + 7.0 * std::sqrt(6.0)/360.0,
                                -2.0/225.0 - std::sqrt(6.0)/75.0
                            }},
                            {{
                                4.0/9.0 - std::sqrt(6.0)/36.0,
                                4.0/9.0 + std::sqrt(6.0)/36.0,
                                1.0/9.0
                            }}
                        }},
                        {{
                            2.0/5.0 - std::sqrt(6.0)/10.0,
                            2.0/5.0 + std::sqrt(6.0)/10.0,
                            1.0
                        }},
                        {{
                            4.0/9.0 - std::sqrt(6.0)/36.0,
                            4.0/9.0 + std::sqrt(6.0)/36.0,
                            1.0/9.0
                        }}
                    }}
    { }

    constexpr Radau_IIA_fifth_order_method
        (const Radau_IIA_fifth_order_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Radau_IIA_fifth_order_method
        (Radau_IIA_fifth_order_method<T_arg_u, T_arg_t>&&) noexcept = default;

    virtual ~Radau_IIA_fifth_order_method() noexcept = default;

    constexpr Radau_IIA_fifth_order_method<T_arg_u, T_arg_t>&
        operator=(const Radau_IIA_fifth_order_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Radau_IIA_fifth_order_method<T_arg_u, T_arg_t>&
        operator=(Radau_IIA_fifth_order_method<T_arg_u, T_arg_t>&&) noexcept = default;

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

#endif // RADAU_IIA_FIFTH_ORDER_METHOD_HXX_INCLUDED
