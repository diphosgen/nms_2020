#ifndef RICHARDSON_EXTRAPOLATION_HXX_INCLUDED
#define RICHARDSON_EXTRAPOLATION_HXX_INCLUDED

#include "matrix_unfixed.hxx"
#include "slae_solver_factory.hxx"
#include "accuracy_sequence_factory.hxx"
#include "series_acceleration_method.hxx"

namespace Math_structures
{

template<
    typename T_arg_u = double,
    typename T_arg_t = double
    >
class Richardson_extrapolation
    :   public Series_acceleration_method<T_arg_u, T_arg_t>
{
private:

    using Ode_solver_ptr_t =
        std::experimental::propagate_const<
            std::unique_ptr<One_step_ODE_solver<T_arg_u, T_arg_t>>
            >;

    using Accuracy_sequence_t =
        std::experimental::propagate_const<
            std::unique_ptr<Accuracy_sequence>
            >;

    using Slae_solver_ptr_t =
        std::experimental::propagate_const<
            std::unique_ptr<SLAE_solver<
                Matrix_unfixed<T_arg_t>,
                Vector_unfixed<T_arg_u>>
                >>;

private:

    Ode_solver_ptr_t ode_solver_ptr = nullptr;
    const int accuracy_order_boost = 0;
    const Accuracy_sequence_t accuracy_sequence = nullptr;
    const Slae_solver_ptr_t slae_solver_ptr = nullptr;

public:

    constexpr Richardson_extrapolation() noexcept = default;

    constexpr Richardson_extrapolation(const Richardson_extrapolation<T_arg_u, T_arg_t>&) noexcept = default;
    constexpr Richardson_extrapolation(Richardson_extrapolation<T_arg_u, T_arg_t>&&) noexcept = default;

    constexpr Richardson_extrapolation(Ode_solver_ptr_t ode_solver_ptr,
                                       const int accuracy_order_boost,
                                       Accuracy_sequence_t accuracy_sequence =
                                            Accuracy_sequence_default_factory{}.create(),
                                       Slae_solver_ptr_t slae_solver_ptr =
                                            SLAE_solver_default_factory<
                                                        Matrix_unfixed<T_arg_t>,
                                                        Vector_unfixed<T_arg_u>
                                                        >{}.create()) noexcept
        :   ode_solver_ptr{std::move(ode_solver_ptr)},
            accuracy_order_boost{accuracy_order_boost},
            accuracy_sequence{std::move(accuracy_sequence)},
            slae_solver_ptr{std::move(slae_solver_ptr)}
    {
		assert(this->accuracy_order_boost >= 0);
	}

    virtual ~Richardson_extrapolation() noexcept = default;

    constexpr Richardson_extrapolation<T_arg_u, T_arg_t>&
        operator=(const Richardson_extrapolation<T_arg_u, T_arg_t>&) noexcept = default;

    constexpr Richardson_extrapolation<T_arg_u, T_arg_t>&
        operator=(Richardson_extrapolation<T_arg_u, T_arg_t>&&) noexcept = default;

    virtual int get_accuracy_order() const noexcept override
    {
        return ode_solver_ptr->get_accuracy_order() + this->accuracy_order_boost;
    }

    virtual T_arg_u get_solution(const std::function<T_arg_u(T_arg_u, T_arg_t)>& f_rhs,
                                 const T_arg_u& u_prev, const T_arg_t& t, const T_arg_t& dt) override
    {
        if (this->accuracy_order_boost == 0) {

            return this->ode_solver_ptr->get_solution(f_rhs, u_prev, t, dt);

        } else {

            const int system_size = this->accuracy_order_boost + 1;

            Matrix_unfixed<T_arg_t> matrix_slae(system_size, system_size);

            for (int i = 0; i < system_size; ++i) {
                for (int j = 0; j < system_size; ++j) {
                    if (j == 0) {
                        matrix_slae[i][j] = T_arg_t(1);
                    } else {
                        const int denominator = this->accuracy_sequence->get_sequence_element(i);
                        const int exponent = this->ode_solver_ptr->get_accuracy_order() + j - 1;
                        const T_arg_t dt_divided = dt/denominator;
                        matrix_slae[i][j] = std::pow(dt_divided, exponent);
                    }
                }
            }

            Vector_unfixed<T_arg_u> vector_slae(system_size, u_prev);

            for (int i = 0; i < system_size; ++i) {
                const int updates_amount = this->accuracy_sequence->get_sequence_element(i);
                const T_arg_t dt_divided = dt/updates_amount;
                T_arg_t t_internal = t;
                for (int update_step = 0; update_step < updates_amount; ++update_step) {
                    vector_slae[i] = this->ode_solver_ptr->
                                     get_solution(f_rhs, vector_slae[i],
                                                  t_internal, dt_divided);
                    t_internal += dt_divided;
                }
            }

            const Vector_unfixed<T_arg_u> vector_res =
                    this->slae_solver_ptr->get_solution(matrix_slae, vector_slae);

            constexpr int solution_component_number = 0;
            return vector_res[solution_component_number];
        }
    }
};

}

#endif // RICHARDSON_EXTRAPOLATION_HXX_INCLUDED
