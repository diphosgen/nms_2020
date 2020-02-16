#ifndef DORMAND_PRINCE_METHOD_HXX_INCLUDED
#define DORMAND_PRINCE_METHOD_HXX_INCLUDED

#include "runge_kutta_method.hxx"
#include "runge_kutta_method_expl.hxx"

namespace Math_structures
{

template<
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class Dormand_Prince_method
    :   public Runge_Kutta_method_expl<7, T_arg_u, T_arg_t>
{
public:

    constexpr Dormand_Prince_method() noexcept
        :   Runge_Kutta_method_expl<7, T_arg_u, T_arg_t>{
                    {5}, {
                        {{
                            {{0.0,            0.0,             0.0,            0.0,          0.0,             0.0,       0.0}},
                            {{0.2,            0.0,             0.0,            0.0,          0.0,             0.0,       0.0}},
                            {{3.0/40.0,       9.0/40.0,        0.0,            0.0,          0.0,             0.0,       0.0}},
                            {{44.0/45.0,      -56.0/15.0,      32.0/9.0,       0.0,          0.0,             0.0,       0.0}},
                            {{19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0, 0.0,             0.0,       0.0}},
                            {{9017.0/3168.0,  -355.0/33.0,     46732.0/5247.0, 49.0/176.0,   -5103.0/18656.0, 0.0,       0.0}},
                            {{35.0/384.0,     0.0,             500.0/1113.0,   125.0/192.0,  -2187.0/6784.0,  11.0/84.0, 0.0}}
                        }},
                        {{0.0,        0.2, 0.3,          0.8,         8.0/9.0,        1.0,       1.0}},
                        {{35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0, 0.0}}
                    }}
    { }

    constexpr Dormand_Prince_method
        (const Dormand_Prince_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Dormand_Prince_method
        (Dormand_Prince_method<T_arg_u, T_arg_t>&&) noexcept = default;

    virtual ~Dormand_Prince_method() noexcept = default;

    constexpr Dormand_Prince_method<T_arg_u, T_arg_t>&
        operator=(const Dormand_Prince_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Dormand_Prince_method<T_arg_u, T_arg_t>&
        operator=(Dormand_Prince_method<T_arg_u, T_arg_t>&&) noexcept = default;

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

#endif // DORMAND_PRINCE_METHOD_HXX_INCLUDED
