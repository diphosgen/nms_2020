#include "task_poisson_equation.hxx"

namespace Tasks
{

void task_poisson_equation() noexcept
{
    std::cout << "|---- task_poisson_equation()" << std::endl;

    constexpr int task_option = 2;

    if constexpr (task_option == 1) {

        poisson_equation_1D_solution();

    } else if constexpr (task_option == 2) {

        poisson_equation_2D_solution();

    } else {

        std::cout << "task_option error" << std::endl;

    }
}

void poisson_equation_1D_solution() noexcept
{
    std::cout << "      |---- poisson_equation_1D_solution()" << std::endl;

    using namespace Math_structures;

    constexpr double lbc_coordinate = 0.0;
    constexpr double rbc_coordinate = 1.0;

    constexpr double coordinate_length = rbc_coordinate - lbc_coordinate;

    constexpr double thermal_diffusivity = 1.0;

    constexpr int grid_steps_amount = 100;

    constexpr double time_of_compute = 5.0;
    constexpr double CFL = 0.2;

    constexpr double sufficient_convergence_accuracy = 1.0e-5;

    constexpr int accuracy_order_boost = 0;

    constexpr Boundary_conditions lbc = Boundary_conditions::Dirichlet;
    constexpr Boundary_conditions rbc = Boundary_conditions::Dirichlet;

    using ODE_solver = One_step_ODE_solver<Regular_grid_1d<double>, double>;
    using RK_method = Classic_fourth_order_method<Regular_grid_1d<double>, double>;
    using Refined_method = Richardson_extrapolation<Regular_grid_1d<double>, double>;

    auto init_conditions = [](const double x)
        noexcept -> double
    {
        return 1.0 + 4.0 * x + std::sin(M_PI * x / coordinate_length);
    };

    auto function_lbc = [](const double t) noexcept -> double
    {
        return 1.0;
    };

    auto function_rbc = [](const double t) noexcept -> double
    {
        return 4.0;
    };

    auto source = [](const double x, const double t) noexcept -> double
    {
        return 10.0 * x * (1.0 - x);
    };

    Heat_conduction_equation_1d<lbc, rbc>
        hce{
            init_conditions,
            function_lbc,
            function_rbc,
            source,
            grid_steps_amount,
            lbc_coordinate,
            rbc_coordinate,
            thermal_diffusivity
        };

    Regular_grid_1d<double> solution_prev_step = hce.get_initial_conditions_grid();
    Regular_grid_1d<double> solution_next_step = hce.get_initial_conditions_grid();

    constexpr double spartial_step = (rbc_coordinate - lbc_coordinate)/grid_steps_amount;
    constexpr double delta_t = CFL * spartial_step * spartial_step / thermal_diffusivity;

    constexpr int time_steps_between_checks = 1000;

    int time_step = 0;
    int output_file_number = 0;
    double t = 0.0;

    const std::unique_ptr<ODE_solver> solver_ptr =
        std::make_unique<Refined_method>(std::make_unique<RK_method>(), accuracy_order_boost);

    double convergence_accuracy = 0.0;

    do {

        if (time_step % time_steps_between_checks == 0) {

            std::stringstream str_stream_file_name{};

            str_stream_file_name << std::fixed << std::setprecision(10)
                                 << "output_data/task_4/"
                                 << output_file_number << "_[t=" << t << "].dat"
                                 << std::ends;

            std::ofstream output_data_stream{str_stream_file_name.str(),
                                             std::ios_base::trunc};

            for (int i = 0; i < solution_prev_step.size(); ++i) {
                output_data_stream << std::scientific << std::setprecision(15)
                                   << lbc_coordinate + i * spartial_step << "\t"
                                   << solution_prev_step[i] << std::endl;
            }

            output_data_stream.flush();
            output_data_stream.close();

            output_file_number++;

        }

        solution_next_step = hce.get_solution(solution_prev_step, solver_ptr, t, delta_t);

        convergence_accuracy = max_abs(solution_next_step - solution_prev_step);

        std::swap(solution_next_step, solution_prev_step);

        time_step++;
        t += delta_t;

        if (t > time_of_compute) {
            std::cout << "A sufficient accuracy of convergence "
                         "for solving the Poisson equation could "
                         "not be obtained for the selected time period "
                         "time_of_compute = " << time_of_compute << "." << std::endl;
        }

    } while (convergence_accuracy > sufficient_convergence_accuracy);

    std::ofstream output_data_stream{"output_data/task_4/convergence_solution.dat",
                                     std::ios_base::trunc};

    for (int i = 0; i < solution_next_step.size(); ++i) {
        output_data_stream << std::scientific << std::setprecision(15)
                           << lbc_coordinate + i * spartial_step << "\t"
                           << solution_next_step[i] << std::endl;
    }

    output_data_stream.flush();
    output_data_stream.close();
}

void poisson_equation_2D_solution() noexcept
{
    std::cout << "      |---- poisson_equation_2D_solution()" << std::endl;

    using namespace Math_structures;

    constexpr double lbc_coordinate_x = 0.0;
    constexpr double rbc_coordinate_x = 1.0;
    constexpr double lbc_coordinate_y = 0.0;
    constexpr double rbc_coordinate_y = 0.5;

    constexpr double coordinate_length_x = rbc_coordinate_x - lbc_coordinate_x;
    constexpr double coordinate_length_y = rbc_coordinate_y - lbc_coordinate_y;

    constexpr double thermal_diffusivity = 1.0;

    constexpr int grid_steps_amount_x = 40;
    constexpr int grid_steps_amount_y = 40;

    constexpr double time_of_compute = 5.0;
    constexpr double CFL = 0.2;

    constexpr double sufficient_convergence_accuracy = 1.0e-5;

    constexpr int accuracy_order_boost = 0;

    constexpr Boundary_conditions lbc_x = Boundary_conditions::Dirichlet;
    constexpr Boundary_conditions rbc_x = Boundary_conditions::Dirichlet;
    constexpr Boundary_conditions lbc_y = Boundary_conditions::Dirichlet;
    constexpr Boundary_conditions rbc_y = Boundary_conditions::Dirichlet;

    using ODE_solver = One_step_ODE_solver<Regular_grid_2d<double>, double>;
    using RK_method = Classic_fourth_order_method<Regular_grid_2d<double>, double>;
    using Refined_method = Richardson_extrapolation<Regular_grid_2d<double>, double>;

    auto init_conditions = [](const double x, const double y)
        noexcept -> double
    {
        return std::sin(M_PI * x / coordinate_length_x) *
               std::sin(M_PI * y / coordinate_length_y);
    };

    auto function_lbc_x = [](const double y, const double t) noexcept -> double
    {
        return 0.0;
    };

    auto function_rbc_x = [](const double y, const double t) noexcept -> double
    {
        return 0.0;
    };

    auto function_lbc_y = [](const double x, const double t) noexcept -> double
    {
        return 0.0;
    };

    auto function_rbc_y = [](const double x, const double t) noexcept -> double
    {
        return 0.0;
    };

    auto source = [](const double x, const double y, const double t) noexcept -> double
    {
        constexpr double q = 10.0;

        constexpr double x_c = 0.5 * (rbc_coordinate_x - lbc_coordinate_x);
        constexpr double y_c = 0.5 * (rbc_coordinate_y - lbc_coordinate_y);

        const double dist_x = x - x_c;
        const double dist_y = y - y_c;

        constexpr double eps = 1.0e-10;

        const double denom = dist_x * dist_x + dist_y * dist_y + eps;

        return q/denom;
    };

    Heat_conduction_equation_2d<lbc_x, rbc_x, lbc_y, rbc_y>
        hce{
            init_conditions,
            function_lbc_x,
            function_rbc_x,
            function_lbc_y,
            function_rbc_y,
            source,
            grid_steps_amount_x,
            grid_steps_amount_y,
            lbc_coordinate_x,
            rbc_coordinate_x,
            lbc_coordinate_y,
            rbc_coordinate_y,
            thermal_diffusivity
        };

    Regular_grid_2d<double> solution_prev_step = hce.get_initial_conditions_grid();
    Regular_grid_2d<double> solution_next_step = hce.get_initial_conditions_grid();

    constexpr double spartial_step_x = (rbc_coordinate_x - lbc_coordinate_x)/grid_steps_amount_x;
    constexpr double spartial_step_y = (rbc_coordinate_y - lbc_coordinate_y)/grid_steps_amount_y;

    constexpr double spartial_step_min = std::min(spartial_step_x, spartial_step_y);
    constexpr double delta_t = CFL * spartial_step_min * spartial_step_min / thermal_diffusivity;

    constexpr int time_steps_between_checks = 1000;

    int time_step = 0;
    int output_file_number = 0;
    double t = 0.0;

    const std::unique_ptr<ODE_solver> solver_ptr =
        std::make_unique<Refined_method>(std::make_unique<RK_method>(), accuracy_order_boost);

    double convergence_accuracy = 0.0;

    do {

        if (time_step % time_steps_between_checks == 0) {

            std::stringstream str_stream_file_name{};

            str_stream_file_name << std::fixed << std::setprecision(10)
                                 << "output_data/task_4/"
                                 << output_file_number << "_[t=" << t << "].dat"
                                 << std::ends;

            std::ofstream output_data_stream{str_stream_file_name.str().c_str(),
                                             std::ios_base::trunc};

            for (int i = 0; i < solution_prev_step.get_size_x(); ++i) {
                for (int j = 0; j < solution_prev_step.get_size_y(); ++j) {
                    output_data_stream << std::scientific << std::setprecision(15)
                                       << lbc_coordinate_x + i * spartial_step_x << "\t"
                                       << lbc_coordinate_y + j * spartial_step_y << "\t"
                                       << solution_prev_step[i][j] << std::endl;
                }
            }

            output_data_stream.flush();
            output_data_stream.close();

            output_file_number++;

        }

        solution_next_step = hce.get_solution(solution_prev_step, solver_ptr, t, delta_t);

        convergence_accuracy = max_abs(solution_next_step - solution_prev_step);

        std::swap(solution_next_step, solution_prev_step);

        time_step++;
        t += delta_t;

        if (t > time_of_compute) {
            std::cout << "A sufficient accuracy of convergence "
                         "for solving the Poisson equation could "
                         "not be obtained for the selected time period "
                         "time_of_compute = " << time_of_compute << "." << std::endl;
        }

    } while (convergence_accuracy > sufficient_convergence_accuracy);

    std::ofstream output_data_stream{"output_data/task_4/convergence_solution.dat",
                                     std::ios_base::trunc};

    for (int i = 0; i < solution_next_step.get_size_x(); ++i) {
        for (int j = 0; j < solution_next_step.get_size_y(); ++j) {
            output_data_stream << std::scientific << std::setprecision(15)
                               << lbc_coordinate_x + i * spartial_step_x << "\t"
                               << lbc_coordinate_y + j * spartial_step_y << "\t"
                               << solution_next_step[i][j] << std::endl;
        }
    }

    output_data_stream.flush();
    output_data_stream.close();
}

}
