#ifndef RUNGE_KUTTA_METHOD_IMPL_HXX_INCLUDED
#define RUNGE_KUTTA_METHOD_IMPL_HXX_INCLUDED

#include "runge_kutta_method.hxx"
#include "snae_solver_factory.hxx"

namespace Math_structures
{

template<
    int stages_amount,
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class Runge_Kutta_method_impl
    :   public Runge_Kutta_method<stages_amount, T_arg_u, T_arg_t>
{
private:

    const std::experimental::propagate_const<
        const std::unique_ptr<SNAE_solver<Vector_fixed<T_arg_u, stages_amount>>>
        > snae_solver_ptr = nullptr;

    virtual void make_increments_calculation(const std::function<T_arg_u(T_arg_u, T_arg_t)>& f_rhs,
                                             const T_arg_u& u_prev, const T_arg_t& t, const T_arg_t& dt) override
    {
        const T_arg_u initial_value_k = f_rhs(u_prev, t);

        for (int stage = 0; stage < stages_amount; ++stage) {
            this->rkm_increments_k[stage] = initial_value_k;
        }

        auto rkm_impl_snae = [this, &f_rhs, &u_prev, &t, &dt]
                             (const Vector_fixed<T_arg_u, stages_amount>& internal_rkm_increments)
                                constexpr -> Vector_fixed<T_arg_u, stages_amount>
        {
            Vector_fixed<T_arg_u, stages_amount> res{};

            for (int stage = 0; stage < stages_amount; ++stage) {

                T_arg_u temp_value = u_prev;

                for (int i = 0; i < stages_amount; ++i) {
                    const T_arg_u temp_increment_value = this->butcher_table.get_element_a_table(stage, i) *
                                                         this->rkm_increments_k[i];
                    temp_value += dt * temp_increment_value;
                }

                const T_arg_t temp_time = t + this->butcher_table.get_element_c_table(stage) * dt;
                res[stage] = f_rhs(temp_value, temp_time) - internal_rkm_increments[stage];
            }

            return res;
        };

        try {
            this->rkm_increments_k = this->snae_solver_ptr->
                                           get_solution(rkm_impl_snae,
                                                        this->rkm_increments_k);
        } catch (const SNAE_solution_not_obtained& ex) {
            throw Runge_Kutta_method_solution_not_obtained{};
        }
    }

public:

    constexpr Runge_Kutta_method_impl() noexcept = delete;

    constexpr Runge_Kutta_method_impl
        (const Runge_Kutta_method_impl<stages_amount, T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Runge_Kutta_method_impl
        (Runge_Kutta_method_impl<stages_amount, T_arg_u, T_arg_t>&&) noexcept = default;

    constexpr Runge_Kutta_method_impl
        (const int accuracy_order, const Butcher_table<stages_amount, T_arg_t>& butcher_table,
         std::unique_ptr<SNAE_solver<Vector_fixed<T_arg_u, stages_amount>>> snae_solver_ptr =
                SNAE_solver_default_factory<Vector_fixed<T_arg_u, stages_amount>>{}.create()) noexcept
        :   Runge_Kutta_method<stages_amount, T_arg_u, T_arg_t>{accuracy_order, butcher_table},
            snae_solver_ptr{std::move(snae_solver_ptr)}
    { }

    constexpr Runge_Kutta_method_impl<stages_amount, T_arg_u, T_arg_t>&
        operator=(const Runge_Kutta_method_impl<stages_amount, T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Runge_Kutta_method_impl<stages_amount, T_arg_u, T_arg_t>&
        operator=(Runge_Kutta_method_impl<stages_amount, T_arg_u, T_arg_t>&&) noexcept = default;

    virtual ~Runge_Kutta_method_impl() noexcept = default;

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

#endif // RUNGE_KUTTA_METHOD_IMPL_HXX_INCLUDED
