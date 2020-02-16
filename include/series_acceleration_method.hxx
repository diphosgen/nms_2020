#ifndef SERIES_ACCELERATION_METHOD_HXX_INCLUDED
#define SERIES_ACCELERATION_METHOD_HXX_INCLUDED

#include "one_step_ode_solver.hxx"

namespace Math_structures
{

template<
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class Series_acceleration_method
    :   public One_step_ODE_solver<T_arg_u, T_arg_t>
{
public:

    constexpr Series_acceleration_method() noexcept = default;

    constexpr Series_acceleration_method(const Series_acceleration_method<T_arg_u, T_arg_t>&) noexcept = default;
    constexpr Series_acceleration_method(Series_acceleration_method<T_arg_u, T_arg_t>&&) noexcept = default;

    virtual ~Series_acceleration_method() noexcept = default;

    constexpr Series_acceleration_method<T_arg_u, T_arg_t>&
        operator=(const Series_acceleration_method<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Series_acceleration_method<T_arg_u, T_arg_t>&
        operator=(Series_acceleration_method<T_arg_u, T_arg_t>&&) noexcept = default;

    virtual int get_accuracy_order() const noexcept override = 0;

    virtual T_arg_u get_solution(const std::function<T_arg_u(T_arg_u, T_arg_t)>& f_rhs,
                                 const T_arg_u& u_prev, const T_arg_t& t, const T_arg_t& dt) override = 0;
};

}

#endif // SERIES_ACCELERATION_METHOD_HXX_INCLUDED
