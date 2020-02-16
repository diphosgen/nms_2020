#include "test_wave_equation_1d.hxx"

namespace Tests
{

void test_wave_equation_1D_6() noexcept
{
    std::cout << "      |---- test_wave_equation_1D_6()" << std::endl;

    using namespace Math_structures;

    constexpr int max_harmonic = 4;
    constexpr int harmonic_amplitudes_container_size = max_harmonic;

    constexpr double min_coefficient_value = -10.0;
    constexpr double max_coefficient_value = +10.0;

    constexpr double lbc_coordinate = 0.0;
    constexpr double rbc_coordinate = 1.0;

    static_assert(rbc_coordinate > lbc_coordinate);

    constexpr double coordinate_length = rbc_coordinate - lbc_coordinate;

    constexpr double wave_speed = 1.0;

    constexpr int grid_steps_amount = 500;

    constexpr double time_of_compute = 2.0;
    constexpr double CFL = 0.2;
    constexpr double spartial_step = (rbc_coordinate - lbc_coordinate)/grid_steps_amount;
    constexpr double delta_t = CFL * spartial_step / wave_speed;

    constexpr double max_error = 1.0e-2;

    constexpr Boundary_conditions lbc = Boundary_conditions::Dirichlet;
    constexpr Boundary_conditions rbc = Boundary_conditions::Robin;

    auto eigenvalue_harmonic_factor = [](const int harmonic_number)
        constexpr noexcept -> double
    {
        return double(harmonic_number) + 0.5;
    };

    auto basis_function = [coordinate_length,&eigenvalue_harmonic_factor]
                          (const double x, const int harmonic_number)
        constexpr noexcept -> double
    {
        return std::sin(M_PI * eigenvalue_harmonic_factor(harmonic_number) * x / coordinate_length);
    };

    using Grid_variable = typename Wave_equation_1d<lbc, rbc>::Grid_variable;

    constexpr int grid_variable_displacement_number =
        Wave_equation_1d<lbc, rbc>::grid_variable_displacement_number;

    constexpr int grid_variable_velocity_number =
        Wave_equation_1d<lbc, rbc>::grid_variable_velocity_number;

    using RKM = Classic_fourth_order_method<Regular_grid_1d<Grid_variable>, double>;

    std::mt19937 gen_rand(0);
    std::uniform_real_distribution<double> distr_amplitudes(min_coefficient_value,
                                                            max_coefficient_value);

    auto rand_amplitude = [&gen_rand, &distr_amplitudes]() noexcept -> double
    {
        return distr_amplitudes(gen_rand);
    };

    std::array<Grid_variable, harmonic_amplitudes_container_size> amplitudes{};

    for (int harmonic = 1; harmonic <= max_harmonic; ++harmonic) {
        amplitudes[harmonic-1][grid_variable_displacement_number] = rand_amplitude();
        amplitudes[harmonic-1][grid_variable_velocity_number] = rand_amplitude();
    }

    auto init_conditions = [max_harmonic, &amplitudes, &basis_function](const double x)
        noexcept -> Grid_variable
    {
        Grid_variable res = {{0.0, 0.0}};

        for (int harmonic = 1; harmonic <= max_harmonic; ++harmonic) {
            res[grid_variable_displacement_number] +=
                amplitudes[harmonic-1][grid_variable_displacement_number] * basis_function(x, harmonic);
            res[grid_variable_velocity_number] +=
                amplitudes[harmonic-1][grid_variable_velocity_number] * basis_function(x, harmonic);
        }

        return res;
    };

    auto function_lbc = [](const double t) noexcept -> double
    {
        return 0.0;
    };

    auto function_rbc_m = [](const double t) noexcept -> double
    {
        return 0.0;
    };

    auto function_rbc_c = [](const double t) noexcept -> double
    {
        return 0.0;
    };

    auto source = [](const double x, const double t) noexcept -> double
    {
        return 0.0;
    };

    auto exact_solution = [wave_speed, max_harmonic, coordinate_length,
                           &amplitudes, &basis_function, &eigenvalue_harmonic_factor]
                          (const double x, const double t) noexcept -> Grid_variable
    {
        Grid_variable res = {{0.0, 0.0}};

        for (int harmonic = 1; harmonic <= max_harmonic; ++harmonic) {

            const double eigenvalue_factor = eigenvalue_harmonic_factor(harmonic);

            const double harmonic_factor =
                M_PI * eigenvalue_factor / coordinate_length;

            const double velocity_factor = harmonic_factor * wave_speed;
            const double time_phase = velocity_factor * t;

            const double sin_fourier_series_factor =
                coordinate_length / (M_PI * eigenvalue_factor * wave_speed);

            const Grid_variable amplitude = amplitudes[harmonic-1];

            const Grid_variable time_factor_component_cos =
                amplitude[grid_variable_displacement_number] *
                Grid_variable{{
                    std::cos(time_phase),
                    -velocity_factor * std::sin(time_phase)
                }};

            const Grid_variable time_factor_component_sin =
                sin_fourier_series_factor *
                amplitude[grid_variable_velocity_number] *
                Grid_variable{{
                    std::sin(time_phase),
                    velocity_factor * std::cos(time_phase)
                }};

            const Grid_variable time_factor = time_factor_component_cos +
                                              time_factor_component_sin;

            const double basis_function_value = basis_function(x, harmonic);

            res += time_factor * basis_function_value;

        }

        return res;
    };

    auto grid_exact_solution = [rbc_coordinate, lbc_coordinate, grid_steps_amount,
                                spartial_step, &exact_solution](const double t)
        noexcept -> Regular_grid_1d<Grid_variable>
    {
        Regular_grid_1d<Grid_variable> grid_solution(grid_steps_amount + 1);

        for (int n = 0; n <= grid_steps_amount; ++n) {
            const double x = lbc_coordinate + n * spartial_step;
            grid_solution[n] = exact_solution(x, t);
        }

        return grid_solution;
    };

    Wave_equation_1d<lbc, rbc>
        wave_equation{
            init_conditions,
            function_lbc,
            std::make_pair(function_rbc_c, function_rbc_m),
            source,
            grid_steps_amount,
            lbc_coordinate,
            rbc_coordinate,
            wave_speed
        };

    Regular_grid_1d<Grid_variable> solution = wave_equation.get_initial_conditions_grid();

    constexpr int time_steps_between_checks = 1000;

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
