#ifndef CESCHINO_KUNZMANN_THIRD_ORDER_DIRK_HXX_INCLUDED
#define CESCHINO_KUNZMANN_THIRD_ORDER_DIRK_HXX_INCLUDED

#include "runge_kutta_method.hxx"
#include "runge_kutta_method_impl.hxx"

namespace Math_structures
{

template<
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class Ceschino_Kunzmann_third_order_DIRK
    :   public Runge_Kutta_method_impl<3, T_arg_u, T_arg_t>
{
public:

    constexpr Ceschino_Kunzmann_third_order_DIRK() noexcept
        :   Runge_Kutta_method_impl<3, T_arg_u, T_arg_t>{
                    {3}, {
                        {{
                            {{0.0,     0.0,        0.0}},
                            {{1.0/4.0, 1.0/4.0,    0.0}},
                            {{1.0/6.0, 4.0/6.0,    1.0/6.0}}
                        }},
                        {{0.0,     1.0/2.0, 1.0}},
                        {{1.0/6.0, 4.0/6.0, 1.0/6.0}}
                    }}
    { }

    constexpr Ceschino_Kunzmann_third_order_DIRK
        (const Ceschino_Kunzmann_third_order_DIRK<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Ceschino_Kunzmann_third_order_DIRK
        (Ceschino_Kunzmann_third_order_DIRK<T_arg_u, T_arg_t>&&) noexcept = default;

    virtual ~Ceschino_Kunzmann_third_order_DIRK() noexcept = default;

    constexpr Ceschino_Kunzmann_third_order_DIRK<T_arg_u, T_arg_t>&
        operator=(const Ceschino_Kunzmann_third_order_DIRK<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Ceschino_Kunzmann_third_order_DIRK<T_arg_u, T_arg_t>&
        operator=(Ceschino_Kunzmann_third_order_DIRK<T_arg_u, T_arg_t>&&) noexcept = default;

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

#endif // CESCHINO_KUNZMANN_THIRD_ORDER_DIRK_HXX_INCLUDED
