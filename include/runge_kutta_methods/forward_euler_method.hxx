#ifndef FORWARD_EULER_METHOD_HXX_INCLUDED
#define FORWARD_EULER_METHOD_HXX_INCLUDED

#include "runge_kutta_method.hxx"
#include "runge_kutta_method_expl.hxx"

namespace Math_structures
{

template<
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class Forward_Euler_method
    :   public Runge_Kutta_method_expl<1, T_arg_u, T_arg_t>
{
public:

    constexpr Forward_Euler_method() noexcept
        :   Runge_Kutta_method_expl<1, T_arg_u, T_arg_t>{
                    {1}, {
                        {{
                            {{0.0}}
                        }},
                        {{0.0}},
                        {{1.0}}
                    }}
    { }

    constexpr Forward_Euler_method
        (const Forward_Euler_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Forward_Euler_method
        (Forward_Euler_method<T_arg_u, T_arg_t>&&) noexcept = default;

    virtual ~Forward_Euler_method() noexcept = default;

    constexpr Forward_Euler_method<T_arg_u, T_arg_t>&
        operator=(const Forward_Euler_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Forward_Euler_method<T_arg_u, T_arg_t>&
        operator=(Forward_Euler_method<T_arg_u, T_arg_t>&&) noexcept = default;

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

#endif // FORWARD_EULER_METHOD_HXX_INCLUDED
