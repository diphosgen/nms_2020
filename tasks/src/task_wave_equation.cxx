#include "task_wave_equation.hxx"

namespace Tasks
{

void task_wave_equation() noexcept
{
    std::cout << "|---- task_wave_equation()" << std::endl;

    constexpr int task_option = 2;

    if constexpr (task_option == 1) {

        wave_equation_1D_solution();

    } else if constexpr (task_option == 2) {

        wave_equation_2D_solution();

    } else {

        std::cout << "task_option error" << std::endl;

    }
}

void wave_equation_1D_solution() noexcept
{
    std::cout << "      |---- wave_equation_1D_solution()" << std::endl;

    using namespace Math_structures;

    constexpr double lbc_coordinate = 0.0;
    constexpr double rbc_coordinate = 1.0;

    constexpr double coordinate_length = rbc_coordinate - lbc_coordinate;

    constexpr double thermal_diffusivity = 1.0;

    constexpr int grid_steps_amount = 100;

    constexpr double time_of_compute = 5.0;
    constexpr double CFL = 0.2;

    constexpr int accuracy_order_boost = 0;

    constexpr Boundary_conditions lbc = Boundary_conditions::Dirichlet;
    constexpr Boundary_conditions rbc = Boundary_conditions::Dirichlet;

    using Grid_variable = typename Wave_equation_1d<lbc, rbc>::Grid_variable;

    using ODE_solver = One_step_ODE_solver<Regular_grid_1d<Grid_variable>, double>;
    using RK_method = Classic_fourth_order_method<Regular_grid_1d<Grid_variable>, double>;
    using Refined_method = Richardson_extrapolation<Regular_grid_1d<Grid_variable>, double>;

    auto init_conditions = [](const double x)
        noexcept -> Grid_variable
    {
        return {{std::sin(M_PI * x / coordinate_length), 0.0}};
    };

    auto function_lbc = [](const double t) noexcept -> double
    {
        return 0.0;
    };

    auto function_rbc = [](const double t) noexcept -> double
    {
        return 0.0;
    };

    auto source = [](const double x, const double t) noexcept -> double
    {
        return 0.0;
    };

    Wave_equation_1d<lbc, rbc>
        wave_equation{
            init_conditions,
            function_lbc,
            function_rbc,
            source,
            grid_steps_amount,
            lbc_coordinate,
            rbc_coordinate,
            thermal_diffusivity
        };

    Regular_grid_1d<Grid_variable> solution = wave_equation.get_initial_conditions_grid();

    constexpr int grid_variable_displacement_number =
        Wave_equation_1d<lbc, rbc>::grid_variable_displacement_number;

    constexpr int grid_variable_velocity_number =
        Wave_equation_1d<lbc, rbc>::grid_variable_velocity_number;

    constexpr double spartial_step = (rbc_coordinate - lbc_coordinate)/grid_steps_amount;
    constexpr double delta_t = CFL * spartial_step * spartial_step / thermal_diffusivity;

    constexpr int time_steps_between_checks = 1'000;

    int time_step = 0;
    int output_file_number = 0;
    double t = 0.0;

    const std::unique_ptr<ODE_solver> solver_ptr =
        std::make_unique<Refined_method>(std::make_unique<RK_method>(), accuracy_order_boost);

    while (t <= time_of_compute) {

        if (time_step % time_steps_between_checks == 0) {

            std::stringstream str_stream_file_name{};

            str_stream_file_name << std::fixed << std::setprecision(10)
                                 << "output_data/task_3/"
                                 << output_file_number << "_[t=" << t << "].dat"
                                 << std::ends;

            std::ofstream output_data_stream{str_stream_file_name.str().c_str(),
                                             std::ios_base::trunc};

            for (int i = 0; i < solution.size(); ++i) {
                output_data_stream << std::scientific << std::setprecision(15)
                                   << lbc_coordinate + i * spartial_step << "\t"
                                   << solution[i][grid_variable_displacement_number] << "\t"
                                   << solution[i][grid_variable_velocity_number] << std::endl;
            }

            output_data_stream.flush();
            output_data_stream.close();

            output_file_number++;

        }

        solution = wave_equation.get_solution(solution, solver_ptr, t, delta_t);

        time_step++;
        t += delta_t;

    }
}

void wave_equation_2D_solution() noexcept
{
    std::cout << "      |---- wave_equation_2D_solution()" << std::endl;

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

    constexpr int accuracy_order_boost = 0;

    constexpr Boundary_conditions lbc_x = Boundary_conditions::Dirichlet;
    constexpr Boundary_conditions rbc_x = Boundary_conditions::Dirichlet;
    constexpr Boundary_conditions lbc_y = Boundary_conditions::Dirichlet;
    constexpr Boundary_conditions rbc_y = Boundary_conditions::Dirichlet;

    using Grid_variable = typename Wave_equation_2d<lbc_x, rbc_x, lbc_y, rbc_y>::Grid_variable;

    using ODE_solver = One_step_ODE_solver<Regular_grid_2d<Grid_variable>, double>;
    using RK_method = Classic_fourth_order_method<Regular_grid_2d<Grid_variable>, double>;
    using Refined_method = Richardson_extrapolation<Regular_grid_2d<Grid_variable>, double>;

    auto init_conditions = [](const double x, const double y)
        noexcept -> Grid_variable
    {
        return {{
                    std::sin(M_PI * x / coordinate_length_x) *
                    std::sin(M_PI * y / coordinate_length_y),
                    0.0
               }};
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
        return 0.0;
    };

    Wave_equation_2d<lbc_x, rbc_x, lbc_y, rbc_y>
        wave_equation{
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

    Regular_grid_2d<Grid_variable> solution = wave_equation.get_initial_conditions_grid();

    constexpr int grid_variable_displacement_number =
        Wave_equation_2d<lbc_x, rbc_x, lbc_y, rbc_y>::grid_variable_displacement_number;

    constexpr int grid_variable_velocity_number =
        Wave_equation_2d<lbc_x, rbc_x, lbc_y, rbc_y>::grid_variable_velocity_number;

    constexpr double spartial_step_x = (rbc_coordinate_x - lbc_coordinate_x)/grid_steps_amount_x;
    constexpr double spartial_step_y = (rbc_coordinate_y - lbc_coordinate_y)/grid_steps_amount_y;

    constexpr double spartial_step_min = std::min(spartial_step_x, spartial_step_y);
    constexpr double delta_t = CFL * spartial_step_min * spartial_step_min / thermal_diffusivity;

    constexpr int time_steps_between_checks = 10'000;

    int time_step = 0;
    int output_file_number = 0;
    double t = 0.0;

    const std::unique_ptr<ODE_solver> solver_ptr =
        std::make_unique<Refined_method>(std::make_unique<RK_method>(), accuracy_order_boost);

    while (t <= time_of_compute) {

        if (time_step % time_steps_between_checks == 0) {

            std::stringstream str_stream_file_name{};

            str_stream_file_name << std::fixed << std::setprecision(10)
                                 << "output_data/task_3/"
                                 << output_file_number << "_[t=" << t << "].dat"
                                 << std::ends;

            std::ofstream output_data_stream{str_stream_file_name.str().c_str(),
                                             std::ios_base::trunc};

            for (int i = 0; i < solution.get_size_x(); ++i) {
                for (int j = 0; j < solution.get_size_y(); ++j) {
                    output_data_stream << std::scientific << std::setprecision(15)
                                       << lbc_coordinate_x + i * spartial_step_x << "\t"
                                       << lbc_coordinate_y + j * spartial_step_y << "\t"
                                       << solution[i][j][grid_variable_displacement_number] << "\t"
                                       << solution[i][j][grid_variable_velocity_number] << std::endl;
                }
            }

            output_data_stream.flush();
            output_data_stream.close();

            output_file_number++;

        }

        solution = wave_equation.get_solution(solution, solver_ptr, t, delta_t);

        time_step++;
        t += delta_t;

    }
}

}
