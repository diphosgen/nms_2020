#ifndef PARESCHI_RUSSO_SECOND_ORDER_DIRK_HXX_INCLUDED
#define PARESCHI_RUSSO_SECOND_ORDER_DIRK_HXX_INCLUDED

#include "runge_kutta_method.hxx"
#include "runge_kutta_method_impl.hxx"

namespace Math_structures
{

template<
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class Pareschi_Russo_second_order_DIRK
    :   public Runge_Kutta_method_impl<2, T_arg_u, T_arg_t>
{
public:

    constexpr Pareschi_Russo_second_order_DIRK
        (const T_arg_t& x = 1.0 - std::sqrt(2.0)/2.0) noexcept
        :   Runge_Kutta_method_impl<2, T_arg_u, T_arg_t>{
                    {2}, {
                        {{
                            {{x,             0.0}},
                            {{1.0 - 2.0 * x, x}}
                        }},
                        {{x,    1.0 - x}},
                        {{0.5,  0.5}}
                    }}
    { }

    constexpr Pareschi_Russo_second_order_DIRK
        (const Pareschi_Russo_second_order_DIRK<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Pareschi_Russo_second_order_DIRK
        (Pareschi_Russo_second_order_DIRK<T_arg_u, T_arg_t>&&) noexcept = default;

    virtual ~Pareschi_Russo_second_order_DIRK() noexcept = default;

    constexpr Pareschi_Russo_second_order_DIRK<T_arg_u, T_arg_t>&
        operator=(const Pareschi_Russo_second_order_DIRK<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Pareschi_Russo_second_order_DIRK<T_arg_u, T_arg_t>&
        operator=(Pareschi_Russo_second_order_DIRK<T_arg_u, T_arg_t>&&) noexcept = default;

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

#endif // PARESCHI_RUSSO_SECOND_ORDER_DIRK_HXX_INCLUDED
