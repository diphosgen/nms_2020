#ifndef NORSETT_THIRD_ORDER_DIRK_HXX_INCLUDED
#define NORSETT_THIRD_ORDER_DIRK_HXX_INCLUDED

#include "runge_kutta_method.hxx"
#include "runge_kutta_method_impl.hxx"

namespace Math_structures
{

template<
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class Norsett_third_order_DIRK
    :   public Runge_Kutta_method_impl<4, T_arg_u, T_arg_t>
{
public:

    constexpr Norsett_third_order_DIRK() noexcept
        :   Runge_Kutta_method_impl<4, T_arg_u, T_arg_t>{
                    {3}, {
                        {{
                            {{1.0/2.0,  0.0,      0.0,     0.0}},
                            {{1.0/6.0,  1.0/2.0,  0.0,     0.0}},
                            {{-1.0/2.0, 1.0/2.0,  1.0/2.0, 0.0}},
                            {{3.0/2.0,  -3.0/2.0, 1.0/2,   1.0/2.0}}
                        }},
                        {{1.0/2.0, 2.0/3.0,  1.0/2.0, 1.0}},
                        {{3.0/2.0, -3.0/2.0, 1.0/2.0, 1.0/2.0}}
                    }}
    { }

    constexpr Norsett_third_order_DIRK
        (const Norsett_third_order_DIRK<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Norsett_third_order_DIRK
        (Norsett_third_order_DIRK<T_arg_u, T_arg_t>&&) noexcept = default;

    virtual ~Norsett_third_order_DIRK() noexcept = default;

    constexpr Norsett_third_order_DIRK<T_arg_u, T_arg_t>&
        operator=(const Norsett_third_order_DIRK<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Norsett_third_order_DIRK<T_arg_u, T_arg_t>&
        operator=(Norsett_third_order_DIRK<T_arg_u, T_arg_t>&&) noexcept = default;

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

#endif // NORSETT_THIRD_ORDER_DIRK_HXX_INCLUDED
