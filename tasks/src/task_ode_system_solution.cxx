#include "task_ode_system_solution.hxx"

namespace Tasks
{

void task_ode_system_solution() noexcept
{
    std::cout << "|---- task_ode_system_solution()" << std::endl;

    harmonic_oscillator_solution();
    individual_ode_system_solution();
}

void harmonic_oscillator_solution() noexcept
{
    std::cout << "      |---- harmonic_oscillator_solution()" << std::endl;

    using namespace Math_structures;

    constexpr int equations_amount = 2;
    using T_arg = Vector_fixed<double, equations_amount>;

    using RK_method = Classic_fourth_order_method<T_arg, double>;
    using Refined_method = Richardson_extrapolation<T_arg, double>;

    constexpr double frequency = 1.0;
    constexpr double t_end = 100.0;
    constexpr int time_steps_amount = 10'000;
    constexpr int time_steps_between_outputs = 1;

    constexpr int accuracy_order_boost = 5;

    auto ode_system =
        [frequency](const T_arg& u, const double t)
        constexpr noexcept -> T_arg
    {
        constexpr double frequency_sqr = frequency * frequency;
        return {{u[1], -frequency_sqr * u[0]}};
    };

    auto initial_conditions = [frequency]()
        constexpr noexcept -> T_arg
    {
        return {{1.0, 0.0}};
    };

    auto exact_solution = [frequency](const double t)
        constexpr noexcept -> T_arg
    {
        return {{std::cos(frequency * t), -std::sin(frequency * t)}};
    };

    RK_method rkm_solver = RK_method{};
    Refined_method refined_solver = Refined_method{
                                        std::make_unique<RK_method>(),
                                        accuracy_order_boost
                                    };

    std::ofstream output_data_stream{"output_data/task_1/harmonic_oscillator_solution.dat"};

    T_arg rkm_solution = initial_conditions();
    T_arg refined_solution = initial_conditions();

    constexpr double dt = t_end / time_steps_amount;

    for (int time_step = 0; time_step <= time_steps_amount; ++time_step) {

        const double t = dt * time_step;

        if (time_step % time_steps_between_outputs == 0) {
            const T_arg res_exact_solution = exact_solution(t);
            output_data_stream << std::fixed << std::setprecision(15) << t << "\t";
            for (int i = 0; i < equations_amount; ++i) {
                output_data_stream << "\t" << rkm_solution[i];
            }
            output_data_stream << "\t";
            for (int i = 0; i < equations_amount; ++i) {
                output_data_stream << "\t" << refined_solution[i];
            }
            output_data_stream << "\t";
            for (int i = 0; i < equations_amount; ++i) {
                output_data_stream << "\t" << res_exact_solution[i];
            }
            const double rkm_solution_error = max_abs(rkm_solution - res_exact_solution);
            const double refined_solution_error = max_abs(refined_solution - res_exact_solution);
            output_data_stream << std::scientific << "\t"
                               << rkm_solution_error << "\t"
                               << refined_solution_error << std::endl;
        }

        rkm_solution = rkm_solver.get_solution(ode_system, rkm_solution, t, dt);
        refined_solution = refined_solver.get_solution(ode_system, refined_solution, t, dt);
    }

    output_data_stream.flush();
    output_data_stream.close();
}

void individual_ode_system_solution() noexcept
{
    std::cout << "      |---- individual_ode_system_solution()" << std::endl;

    using namespace Math_structures;

    constexpr int equations_amount = 2;
    using T_arg = Vector_fixed<double, equations_amount>;

    using RK_method = Classic_fourth_order_method<T_arg, double>;
    using Refined_method = Richardson_extrapolation<T_arg, double>;

    constexpr double frequency = 1.0;
    constexpr double t_end = 100.0;
    constexpr int time_steps_amount = 10'000;
    constexpr int time_steps_between_outputs = 1;

    constexpr int accuracy_order_boost = 5;

    auto individual_ode_system =
        [frequency](const T_arg& u, const double t)
        constexpr noexcept -> T_arg
    {
        constexpr double frequency_sqr = frequency * frequency;
        return {{u[1], -frequency_sqr * u[0]}};
    };

    auto initial_conditions = [frequency]()
        constexpr noexcept -> T_arg
    {
        return {{1.0, 0.0}};
    };

    RK_method rkm_solver = RK_method{};
    Refined_method refined_solver = Refined_method{
                                        std::make_unique<RK_method>(),
                                        accuracy_order_boost
                                    };

    std::ofstream output_data_stream{"output_data/task_1/individual_ode_system_solution.dat"};

    T_arg rkm_solution = initial_conditions();
    T_arg refined_solution = initial_conditions();

    constexpr double dt = t_end / time_steps_amount;

    for (int time_step = 0; time_step <= time_steps_amount; ++time_step) {

        const double t = dt * time_step;

        if (time_step % time_steps_between_outputs == 0) {
            output_data_stream << std::scientific << std::setprecision(15) << t << "\t";
            for (int i = 0; i < equations_amount; ++i) {
                output_data_stream << "\t" << rkm_solution[i];
            }
            output_data_stream << "\t";
            for (int i = 0; i < equations_amount; ++i) {
                output_data_stream << "\t" << refined_solution[i];
            }
            output_data_stream << std::endl;
        }

        rkm_solution = rkm_solver.get_solution(individual_ode_system, rkm_solution, t, dt);
        refined_solution = refined_solver.get_solution(individual_ode_system, refined_solution, t, dt);
    }

    output_data_stream.flush();
    output_data_stream.close();
}

}
