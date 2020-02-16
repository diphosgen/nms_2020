#ifndef LOBATTO_IIIA_EIGHTH_ORDER_METHOD_HXX_INCLUDED
#define LOBATTO_IIIA_EIGHTH_ORDER_METHOD_HXX_INCLUDED

#include "runge_kutta_method.hxx"
#include "runge_kutta_method_impl.hxx"

namespace Math_structures
{

template<
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class Lobatto_IIIA_eighth_order_method
    :   public Runge_Kutta_method_impl<5, T_arg_u, T_arg_t>
{
public:

    constexpr Lobatto_IIIA_eighth_order_method() noexcept
        :   Runge_Kutta_method_impl<5, T_arg_u, T_arg_t>{
                    {8}, {
                        {{
                            {{0.0, 0.0, 0.0, 0.0, 0.0}},
                            {{
                                (119.0 + 3.0 * std::sqrt(21.0))/1960.0,
                                (343.0 - 9.0 * std::sqrt(21.0))/2520.0,
                                (392.0 - 96.0 * std::sqrt(21.0))/2205.0,
                                (343.0 - 69.0 * std::sqrt(21.0))/2520.0,
                                (-21.0 + 3.0 * std::sqrt(21.0))/1960.0,
                            }},
                            {{
                                13.0/320.0,
                                (392.0 + 105.0 * std::sqrt(21.0))/2880.0,
                                8.0/45.0,
                                (392.0 - 105.0 * std::sqrt(21.0))/2880.0,
                                3.0/320.0
                            }},
                            {{
                                (119.0 - 3.0 * std::sqrt(21.0))/1960.0,
                                (343.0 + 69.0 * std::sqrt(21.0))/2520.0,
                                (392.0 + 96.0 * std::sqrt(21.0))/2205.0,
                                (343.0 + 9.0 * std::sqrt(21.0))/2520.0,
                                (-21.0 - 3.0 * std::sqrt(21.0))/1960.0,
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

    constexpr Lobatto_IIIA_eighth_order_method
        (const Lobatto_IIIA_eighth_order_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Lobatto_IIIA_eighth_order_method
        (Lobatto_IIIA_eighth_order_method<T_arg_u, T_arg_t>&&) noexcept = default;

    virtual ~Lobatto_IIIA_eighth_order_method() noexcept = default;

    constexpr Lobatto_IIIA_eighth_order_method<T_arg_u, T_arg_t>&
        operator=(const Lobatto_IIIA_eighth_order_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Lobatto_IIIA_eighth_order_method<T_arg_u, T_arg_t>&
        operator=(Lobatto_IIIA_eighth_order_method<T_arg_u, T_arg_t>&&) noexcept = default;

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

#endif // LOBATTO_IIIA_EIGHTH_ORDER_METHOD_HXX_INCLUDED
