#ifndef SNAE_SOLVER_FPI_HXX_INCLUDED
#define SNAE_SOLVER_FPI_HXX_INCLUDED

#include "snae_solver.hxx"

namespace Math_structures
{

template<typename T_args>
class SNAE_solver_FPI
    :   public SNAE_solver<T_args>
{
private:

    using Invoke_res_type = T_args;
    using Invoke_arg_type = decltype(std::declval<T_args>());
    using Invoke_abs_type = decltype(max_abs(std::declval<Invoke_res_type>()));

public:

    using T_norm = Invoke_abs_type;

    constexpr SNAE_solver_FPI(const T_norm& precision = T_norm(1.0e-8),
                              const int max_amount_iterations = 1'000'000) noexcept
        :   precision{precision},
            max_amount_iterations{max_amount_iterations}
    { }

    virtual ~SNAE_solver_FPI() noexcept = default;

    virtual T_args get_solution(const std::function<T_args(T_args)>& solve_function,
                                const T_args& initial_value) const override final
    {
        int amount_iterations = 0;

        Invoke_res_type result = initial_value;
        Invoke_res_type iteration_error = solve_function(result);

        while(max_abs(iteration_error) > this->precision) {
            result += iteration_error;
            iteration_error = solve_function(result);
            amount_iterations++;

            if (amount_iterations >= this->max_amount_iterations) {
                throw SNAE_solution_not_obtained{};
            }
        }

        return result;
    }

private:

    T_norm precision = T_norm(0.0);
    int max_amount_iterations = 0;

};

}

#endif // SNAE_SOLVER_FPI_HXX_INCLUDED
