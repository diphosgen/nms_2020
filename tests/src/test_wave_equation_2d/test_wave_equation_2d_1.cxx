#include "test_wave_equation_2d.hxx"

namespace Tests
{

void test_wave_equation_2D_1() noexcept
{
    std::cout << "      |---- test_wave_equation_2D_1()" << std::endl;

    using namespace Math_structures;

    constexpr int max_harmonic_x = 2;
    constexpr int max_harmonic_y = 2;

    constexpr int harmonic_amplitudes_container_size_x = max_harmonic_x;
    constexpr int harmonic_amplitudes_container_size_y = max_harmonic_y;

    constexpr double min_coefficient_value = -10.0;
    constexpr double max_coefficient_value = +10.0;

    constexpr double lbc_coordinate_x = 0.0;
    constexpr double rbc_coordinate_x = 1.0;
    constexpr double lbc_coordinate_y = 0.0;
    constexpr double rbc_coordinate_y = 0.5;

    static_assert(rbc_coordinate_x > lbc_coordinate_x);
    static_assert(rbc_coordinate_y > lbc_coordinate_y);

    constexpr double coordinate_length_x = rbc_coordinate_x - lbc_coordinate_x;
    constexpr double coordinate_length_y = rbc_coordinate_y - lbc_coordinate_y;

    constexpr double wave_speed = 1.0;

    constexpr int grid_steps_amount_x = 100;
    constexpr int grid_steps_amount_y = 80;

    constexpr double time_of_compute = 0.5;
    constexpr double CFL = 0.2;
    constexpr double spartial_step_x = (rbc_coordinate_x - lbc_coordinate_x)/grid_steps_amount_x;
    constexpr double spartial_step_y = (rbc_coordinate_y - lbc_coordinate_y)/grid_steps_amount_y;
    constexpr double spartial_step_min = std::min(spartial_step_x, spartial_step_y);
    constexpr double delta_t = CFL * spartial_step_min / wave_speed;

    constexpr double max_error = 1.0e-1;

    constexpr Boundary_conditions lbc_x = Boundary_conditions::Dirichlet;
    constexpr Boundary_conditions rbc_x = Boundary_conditions::Dirichlet;
    constexpr Boundary_conditions lbc_y = Boundary_conditions::Dirichlet;
    constexpr Boundary_conditions rbc_y = Boundary_conditions::Dirichlet;

    using Grid_variable = typename Wave_equation_2d<lbc_x, rbc_x, lbc_y, rbc_y>::Grid_variable;

    constexpr int grid_variable_displacement_number =
        Wave_equation_2d<lbc_x, rbc_x, lbc_y, rbc_y>::grid_variable_displacement_number;

    constexpr int grid_variable_velocity_number =
        Wave_equation_2d<lbc_x, rbc_x, lbc_y, rbc_y>::grid_variable_velocity_number;

    using RKM = Classic_fourth_order_method<Regular_grid_2d<Grid_variable>, double>;

    auto eigenvalue_harmonic_factor_x = [](const int harmonic_number)
        noexcept -> double
    {
        return double(harmonic_number);
    };

    auto eigenvalue_harmonic_factor_y = [](const int harmonic_number)
        noexcept -> double
    {
        return double(harmonic_number);
    };

    auto basis_function_x = [coordinate_length_x,&eigenvalue_harmonic_factor_x]
                            (const double x, const int harmonic_number)
        noexcept -> double
    {
        const double eigenvalue_factor = eigenvalue_harmonic_factor_x(harmonic_number);
        return std::sin(M_PI * eigenvalue_factor * x / coordinate_length_x);
    };

    auto basis_function_y = [coordinate_length_y,&eigenvalue_harmonic_factor_y]
                            (const double y, const int harmonic_number)
        noexcept -> double
    {
        const double eigenvalue_factor = eigenvalue_harmonic_factor_y(harmonic_number);
        return std::sin(M_PI * eigenvalue_factor * y / coordinate_length_y);
    };

    std::mt19937 gen_rand(0);
    std::uniform_real_distribution<double> distr_amplitudes(min_coefficient_value,
                                                            max_coefficient_value);

    auto rand_amplitude = [&gen_rand, &distr_amplitudes]() noexcept -> double
    {
        return distr_amplitudes(gen_rand);
    };

    std::array<Grid_variable, harmonic_amplitudes_container_size_x> amplitudes_x{};
    std::array<Grid_variable, harmonic_amplitudes_container_size_y> amplitudes_y{};

    for (int harmonic = 1; harmonic <= max_harmonic_x; ++harmonic) {
        amplitudes_x[harmonic-1][grid_variable_displacement_number] = rand_amplitude();
        amplitudes_x[harmonic-1][grid_variable_velocity_number] = rand_amplitude();
    }

    for (int harmonic = 1; harmonic <= max_harmonic_y; ++harmonic) {
        amplitudes_y[harmonic-1][grid_variable_displacement_number] = rand_amplitude();
        amplitudes_y[harmonic-1][grid_variable_velocity_number] = rand_amplitude();
    }

    auto init_conditions = [max_harmonic_x, max_harmonic_y,
                            &amplitudes_x, &amplitudes_y,
                            &basis_function_x, &basis_function_y]
                           (const double x, const double y)
        noexcept -> Grid_variable
    {
        Grid_variable res_x = {{0.0, 0.0}};
        Grid_variable res_y = {{0.0, 0.0}};

        for (int harmonic = 1; harmonic <= max_harmonic_x; ++harmonic) {
            res_x[grid_variable_displacement_number] +=
                amplitudes_x[harmonic-1][grid_variable_displacement_number] *
                basis_function_x(x, harmonic);
            res_x[grid_variable_velocity_number] +=
                amplitudes_x[harmonic-1][grid_variable_velocity_number] *
                basis_function_x(x, harmonic);
        }

        for (int harmonic = 1; harmonic <= max_harmonic_y; ++harmonic) {
            res_y[grid_variable_displacement_number] +=
                amplitudes_y[harmonic-1][grid_variable_displacement_number] *
                basis_function_y(y, harmonic);
            res_y[grid_variable_velocity_number] +=
                amplitudes_y[harmonic-1][grid_variable_velocity_number] *
                basis_function_y(y, harmonic);
        }

        const Grid_variable res = {{
                res_x[grid_variable_displacement_number] *
                res_y[grid_variable_displacement_number],
                res_x[grid_variable_velocity_number] *
                res_y[grid_variable_velocity_number]
            }};

        return res;
    };

    auto function_lbc_x = []([[maybe_unused]] const double y, [[maybe_unused]] const double t) noexcept -> double
    {
        return 0.0;
    };

    auto function_rbc_x = []([[maybe_unused]] const double y, [[maybe_unused]] const double t) noexcept -> double
    {
        return 0.0;
    };

    auto function_lbc_y = []([[maybe_unused]] const double x, [[maybe_unused]] const double t) noexcept -> double
    {
        return 0.0;
    };

    auto function_rbc_y = []([[maybe_unused]] const double x, [[maybe_unused]] const double t) noexcept -> double
    {
        return 0.0;
    };

    auto source = []([[maybe_unused]] const double x, [[maybe_unused]] const double y, [[maybe_unused]] const double t) noexcept -> double
    {
        return 0.0;
    };

    auto exact_solution = [wave_speed, max_harmonic_x, max_harmonic_y,
                           coordinate_length_x, coordinate_length_y,
                           &eigenvalue_harmonic_factor_x, &eigenvalue_harmonic_factor_y,
                           &amplitudes_x, &amplitudes_y, &basis_function_x, &basis_function_y]
                          (const double x, const double y, const double t)
        noexcept -> Grid_variable
    {
        Grid_variable res = {{0.0, 0.0}};

        for (int harmonic_x = 1; harmonic_x <= max_harmonic_x; ++harmonic_x) {
            for (int harmonic_y = 1; harmonic_y <= max_harmonic_y; ++harmonic_y) {
                const double harmonic_factor_x =
                    M_PI * eigenvalue_harmonic_factor_x(harmonic_x) / coordinate_length_x;
                const double harmonic_factor_y =
                    M_PI * eigenvalue_harmonic_factor_y(harmonic_y) / coordinate_length_y;
                const double omega_xy =
                    wave_speed * std::sqrt(harmonic_factor_x * harmonic_factor_x +
                                           harmonic_factor_y * harmonic_factor_y);

                const double basis_functions_tensor_product_component =
                    basis_function_x(x, harmonic_x) *
                    basis_function_y(y, harmonic_y);

                const double amplitudes_displacement_tensor_product_component =
                    amplitudes_x[harmonic_x-1][grid_variable_displacement_number] *
                    amplitudes_y[harmonic_y-1][grid_variable_displacement_number];

                const double amplitudes_velocity_tensor_product_component =
                    amplitudes_x[harmonic_x-1][grid_variable_velocity_number] *
                    amplitudes_y[harmonic_y-1][grid_variable_velocity_number];

                const Grid_variable res_a_xy =
                    amplitudes_displacement_tensor_product_component *
                    basis_functions_tensor_product_component *
                    Grid_variable{{
                        std::cos(omega_xy * t),
                        -omega_xy * std::sin(omega_xy * t)
                    }};

                const Grid_variable res_b_xy =
                    amplitudes_velocity_tensor_product_component *
                    basis_functions_tensor_product_component *
                    Grid_variable{{
                        std::sin(omega_xy * t) / omega_xy,
                        std::cos(omega_xy * t)
                    }};

                res += res_a_xy + res_b_xy;
            }
        }

        return res;
    };

    auto grid_exact_solution = [rbc_coordinate_x, lbc_coordinate_x, grid_steps_amount_x,
                                rbc_coordinate_y, lbc_coordinate_y, grid_steps_amount_y,
                                spartial_step_x, spartial_step_y, &exact_solution](const double t)
        noexcept -> Regular_grid_2d<Grid_variable>
    {
        Regular_grid_2d<Grid_variable> grid_solution(grid_steps_amount_x + 1,
                                                     grid_steps_amount_y + 1);

        for (int i = 0; i <= grid_steps_amount_x; ++i) {
            for (int j = 0; j <= grid_steps_amount_y; ++j) {
                const double x = lbc_coordinate_x + i * spartial_step_x;
                const double y = lbc_coordinate_y + j * spartial_step_y;
                grid_solution[i][j] = exact_solution(x, y, t);
            }
        }

        return grid_solution;
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
            wave_speed
        };

    Regular_grid_2d<Grid_variable> solution = wave_equation.get_initial_conditions_grid();

    constexpr int time_steps_between_checks = 1;

    int time_step = 0;
    double t = 0.0;

    while (t <= time_of_compute) {

        solution = wave_equation.get_solution(solution, std::make_unique<RKM>(), t, delta_t);

        if (time_step % time_steps_between_checks == 0) {
            const auto exact_solution = grid_exact_solution(t + delta_t);
            const double eps = 1.0e-12;
            const double abs_norm = max_abs(solution + exact_solution) + eps;
            const double error = max_abs(solution - exact_solution)/abs_norm;
            assert(error < max_error);
        }

        time_step++;
        t += delta_t;

    }
}

}
