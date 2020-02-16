#ifndef CLASSIC_THREE_EIGHTHS_RULE_METHOD_HXX_INCLUDED
#define CLASSIC_THREE_EIGHTHS_RULE_METHOD_HXX_INCLUDED

#include "runge_kutta_method.hxx"
#include "runge_kutta_method_expl.hxx"

namespace Math_structures
{

template<
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class Classic_three_eighths_rule_method
    :   public Runge_Kutta_method_expl<4, T_arg_u, T_arg_t>
{
public:

    constexpr Classic_three_eighths_rule_method() noexcept
        :   Runge_Kutta_method_expl<4, T_arg_u, T_arg_t>{
                    {4}, {
                        {{
                            {{0.0,      0.0,  0.0, 0.0}},
                            {{1.0/3.0,  0.0,  0.0, 0.0}},
                            {{-1.0/3.0, 1.0,  0.0, 0.0}},
                            {{1.0,      -1.0, 1.0, 0.0}}
                        }},
                        {{0.0,     1.0/3.0, 2.0/3.0, 1.0}},
                        {{1.0/8.0, 3.0/8.0, 3.0/8.0, 1.0/8.0}}
                    }}
    { }

    constexpr Classic_three_eighths_rule_method
        (const Classic_three_eighths_rule_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Classic_three_eighths_rule_method
        (Classic_three_eighths_rule_method<T_arg_u, T_arg_t>&&) noexcept = default;

    virtual ~Classic_three_eighths_rule_method() noexcept = default;

    constexpr Classic_three_eighths_rule_method<T_arg_u, T_arg_t>&
        operator=(const Classic_three_eighths_rule_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Classic_three_eighths_rule_method<T_arg_u, T_arg_t>&
        operator=(Classic_three_eighths_rule_method<T_arg_u, T_arg_t>&&) noexcept = default;

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

#endif // CLASSIC_THREE_EIGHTHS_RULE_METHOD_HXX_INCLUDED
