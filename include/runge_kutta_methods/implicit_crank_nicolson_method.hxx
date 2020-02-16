#ifndef IMPLICIT_CRANK_NICOLSON_METHOD_HXX_INCLUDED
#define IMPLICIT_CRANK_NICOLSON_METHOD_HXX_INCLUDED

#include "runge_kutta_method.hxx"
#include "runge_kutta_method_impl.hxx"

namespace Math_structures
{

template<
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class Implicit_Crank_Nicolson_method
    :   public Runge_Kutta_method_impl<2, T_arg_u, T_arg_t>
{
public:

    constexpr Implicit_Crank_Nicolson_method() noexcept
        :   Runge_Kutta_method_impl<2, T_arg_u, T_arg_t>{
                    {2}, {
                        {{
                            {{0.0, 0.0}},
                            {{0.5, 0.5}}
                        }},
                        {{0.0, 1.0}},
                        {{0.5, 0.5}}
                    }}
    { }

    constexpr Implicit_Crank_Nicolson_method
        (const Implicit_Crank_Nicolson_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Implicit_Crank_Nicolson_method
        (Implicit_Crank_Nicolson_method<T_arg_u, T_arg_t>&&) noexcept = default;

    virtual ~Implicit_Crank_Nicolson_method() noexcept = default;

    constexpr Implicit_Crank_Nicolson_method<T_arg_u, T_arg_t>&
        operator=(const Implicit_Crank_Nicolson_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Implicit_Crank_Nicolson_method<T_arg_u, T_arg_t>&
        operator=(Implicit_Crank_Nicolson_method<T_arg_u, T_arg_t>&&) noexcept = default;

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

#endif // IMPLICIT_CRANK_NICOLSON_METHOD_HXX_INCLUDED
