#ifndef EXPLICIT_SSPRK3_METHOD_HXX_INCLUDED
#define EXPLICIT_SSPRK3_METHOD_HXX_INCLUDED

#include "runge_kutta_method.hxx"
#include "runge_kutta_method_expl.hxx"
#include "runge_kutta_method_impl.hxx"

namespace Math_structures
{

template<
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class Explicit_SSPRK3_method
    :   public Runge_Kutta_method_expl<3, T_arg_u, T_arg_t>
{
public:

    constexpr Explicit_SSPRK3_method() noexcept
        :   Runge_Kutta_method_expl<3, T_arg_u, T_arg_t>{
                    {3}, {
                        {{
                            {{0.0,  0.0,  0.0}},
                            {{1.0,  0.0,  0.0}},
                            {{0.25, 0.25, 0.0}}
                        }},
                        {{0.0,     1.0,     0.5}},
                        {{1.0/6.0, 1.0/6.0, 2.0/3.0}}
                    }}
    { }

    constexpr Explicit_SSPRK3_method
        (const Explicit_SSPRK3_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Explicit_SSPRK3_method
        (Explicit_SSPRK3_method<T_arg_u, T_arg_t>&&) noexcept = default;

    virtual ~Explicit_SSPRK3_method() noexcept = default;

    constexpr Explicit_SSPRK3_method<T_arg_u, T_arg_t>&
        operator=(const Explicit_SSPRK3_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Explicit_SSPRK3_method<T_arg_u, T_arg_t>&
        operator=(Explicit_SSPRK3_method<T_arg_u, T_arg_t>&&) noexcept = default;

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

#endif // EXPLICIT_SSPRK3_METHOD_HXX_INCLUDED
