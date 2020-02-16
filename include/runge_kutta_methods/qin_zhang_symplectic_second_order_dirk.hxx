#ifndef QIN_ZHANG_SYMPLECTIC_SECOND_ORDER_DIRK_HXX_INCLUDED
#define QIN_ZHANG_SYMPLECTIC_SECOND_ORDER_DIRK_HXX_INCLUDED

#include "runge_kutta_method.hxx"
#include "runge_kutta_method_impl.hxx"

namespace Math_structures
{

template<
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class Qin_Zhang_symplectic_second_order_DIRK
    :   public Runge_Kutta_method_impl<2, T_arg_u, T_arg_t>
{
public:

    constexpr Qin_Zhang_symplectic_second_order_DIRK() noexcept
        :   Runge_Kutta_method_impl<2, T_arg_u, T_arg_t>{
                    {2}, {
                        {{
                            {{0.25, 0.00}},
                            {{0.50, 0.25}}
                        }},
                        {{0.25, 0.75}},
                        {{0.50, 0.50}}
                    }}
    { }

    constexpr Qin_Zhang_symplectic_second_order_DIRK
        (const Qin_Zhang_symplectic_second_order_DIRK<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Qin_Zhang_symplectic_second_order_DIRK
        (Qin_Zhang_symplectic_second_order_DIRK<T_arg_u, T_arg_t>&&) noexcept = default;

    virtual ~Qin_Zhang_symplectic_second_order_DIRK() noexcept = default;

    constexpr Qin_Zhang_symplectic_second_order_DIRK<T_arg_u, T_arg_t>&
        operator=(const Qin_Zhang_symplectic_second_order_DIRK<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Qin_Zhang_symplectic_second_order_DIRK<T_arg_u, T_arg_t>&
        operator=(Qin_Zhang_symplectic_second_order_DIRK<T_arg_u, T_arg_t>&&) noexcept = default;

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

#endif // QIN_ZHANG_SYMPLECTIC_SECOND_ORDER_DIRK_HXX_INCLUDED
