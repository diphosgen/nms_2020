#ifndef LOBATTO_IIIC_EIGHTH_ORDER_METHOD_HXX_INCLUDED
#define LOBATTO_IIIC_EIGHTH_ORDER_METHOD_HXX_INCLUDED

#include "runge_kutta_method.hxx"
#include "runge_kutta_method_impl.hxx"

namespace Math_structures
{

template<
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class Lobatto_IIIC_eighth_order_method
    :   public Runge_Kutta_method_impl<5, T_arg_u, T_arg_t>
{
public:

    constexpr Lobatto_IIIC_eighth_order_method() noexcept
        :   Runge_Kutta_method_impl<5, T_arg_u, T_arg_t>{
                    {8}, {
                        {{
                            {{
                                1.0/20.0,
                                -7.0/60.0,
                                2.0/15.0,
                                -7.0/60.0,
                                1.0/20.0
                            }},
                            {{
                                1.0/20.0,
                                29.0/180.0,
                                (47.0 - 15.0 * std::sqrt(21.0))/315.0,
                                (203.0 - 30.0 * std::sqrt(21.0))/1260.0,
                                -3.0/140.0
                            }},
                            {{
                                1.0/20.0,
                                (329.0 + 105.0 * std::sqrt(21.0))/2880.0,
                                73.0/360.0,
                                (329.0 - 105.0 * std::sqrt(21.0))/2880.0,
                                3.0/160.0
                            }},
                            {{
                                1.0/20.0,
                                (203.0 + 30.0 * std::sqrt(21.0))/1260.0,
                                (47.0 + 15.0 * std::sqrt(21.0))/315.0,
                                29.0/180.0,
                                -3.0/140.0
                            }},
                            {{1.0/20.0, 49.0/180.0, 16.0/45.0, 49.0/180.0, 1.0/20.0}}
                        }},
                        {{
                            0.0,
                            1.0/2.0 - std::sqrt(21.0)/14.0,
                            1.0/2.0,
                            1.0/2.0 + std::sqrt(21.0)/14.0,
                            1.0
                        }},
                        {{1.0/20.0, 49.0/180.0, 16.0/45.0, 49.0/180.0, 1.0/20.0}}
                    }}
    { }

    constexpr Lobatto_IIIC_eighth_order_method
        (const Lobatto_IIIC_eighth_order_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Lobatto_IIIC_eighth_order_method
        (Lobatto_IIIC_eighth_order_method<T_arg_u, T_arg_t>&&) noexcept = default;

    virtual ~Lobatto_IIIC_eighth_order_method() noexcept = default;

    constexpr Lobatto_IIIC_eighth_order_method<T_arg_u, T_arg_t>&
        operator=(const Lobatto_IIIC_eighth_order_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Lobatto_IIIC_eighth_order_method<T_arg_u, T_arg_t>&
        operator=(Lobatto_IIIC_eighth_order_method<T_arg_u, T_arg_t>&&) noexcept = default;

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

#endif // LOBATTO_IIIC_EIGHTH_ORDER_METHOD_HXX_INCLUDED
