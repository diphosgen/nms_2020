#include "test_heat_conduction_equation_1d.hxx"

namespace Tests
{

void test_heat_conduction_equation_1D_2() noexcept
{
    std::cout << "      |---- test_heat_conduction_equation_1D_2()" << std::endl;

    using namespace Math_structures;

    constexpr int max_harmonic = 8;
    constexpr int harmonic_amplitudes_container_size = max_harmonic;

    constexpr double min_coefficient_value = -10.0;
    constexpr double max_coefficient_value = +10.0;

    constexpr double lbc_coordinate = 0.0;
    constexpr double rbc_coordinate = 1.0;

    static_assert(rbc_coordinate > lbc_coordinate);

    constexpr double coordinate_length = rbc_coordinate - lbc_coordinate;

    constexpr double thermal_diffusivity = 1.0;

    constexpr int grid_steps_amount = 100;

    constexpr double time_of_compute = 1.0;
    constexpr double CFL = 0.2;
    constexpr double spartial_step = (rbc_coordinate - lbc_coordinate)/grid_steps_amount;
    constexpr double delta_t = CFL * spartial_step * spartial_step / thermal_diffusivity;

    constexpr double max_error = 1.0e-2;

    using RKM = Classic_fourth_order_method<Regular_grid_1d<double>, double>;

    constexpr Boundary_conditions lbc = Boundary_conditions::Neumann;
    constexpr Boundary_conditions rbc = Boundary_conditions::Neumann;

    auto eigenvalue_harmonic_factor = [](const int harmonic_number)
        noexcept -> double
    {
        return double(harmonic_number);
    };

    auto basis_function = [coordinate_length,&eigenvalue_harmonic_factor]
                          (const double x, const int harmonic_number)
        noexcept -> double
    {
        return std::cos(M_PI * eigenvalue_harmonic_factor(harmonic_number) * x / coordinate_length);
    };

    std::mt19937 gen_rand(0);
    std::uniform_real_distribution<double> distr_amplitudes(min_coefficient_value,
                                                            max_coefficient_value);

    auto rand_amplitude = [&gen_rand, &distr_amplitudes]() noexcept -> double
    {
        return distr_amplitudes(gen_rand);
    };

    std::array<double, harmonic_amplitudes_container_size> amplitudes{};

    for (int harmonic = 1; harmonic <= max_harmonic; ++harmonic) {
        amplitudes[harmonic-1] = rand_amplitude();
    }

    auto init_conditions = [max_harmonic, &amplitudes, &basis_function](const double x)
        noexcept -> double
    {
        double res = 0.0;

        for (int harmonic = 1; harmonic <= max_harmonic; ++harmonic) {
            res += amplitudes[harmonic-1] * basis_function(x, harmonic);
        }

        return res;
    };

    auto function_lbc = []([[maybe_unused]] const double t) noexcept -> double
    {
        return 0.0;
    };

    auto function_rbc = []([[maybe_unused]] const double t) noexcept -> double
    {
        return 0.0;
    };

    auto source = []([[maybe_unused]] const double x, [[maybe_unused]] const double t) noexcept -> double
    {
        return 0.0;
    };

    auto exact_solution = [thermal_diffusivity, max_harmonic, coordinate_length,
                           &amplitudes, &basis_function, &eigenvalue_harmonic_factor]
                          (const double x, const double t) noexcept -> double
    {
        double res = 0.0;

        for (int harmonic = 1; harmonic <= max_harmonic; ++harmonic) {
            const double harmonic_factor_sqrt =
                M_PI * eigenvalue_harmonic_factor(harmonic) / coordinate_length;
            const double harmonic_factor = harmonic_factor_sqrt * harmonic_factor_sqrt;
            const double time_factor = std::exp(-thermal_diffusivity * harmonic_factor * t);
            res += amplitudes[harmonic-1] * basis_function(x, harmonic) * time_factor;
        }

        return res;
    };

    auto grid_exact_solution = [rbc_coordinate, lbc_coordinate, grid_steps_amount,
                                spartial_step, &exact_solution](const double t)
        noexcept -> Regular_grid_1d<double>
    {
        Regular_grid_1d<double> grid_solution(grid_steps_amount + 1);

        for (int n = 0; n <= grid_steps_amount; ++n) {
            const double x = lbc_coordinate + n * spartial_step;
            grid_solution[n] = exact_solution(x, t);
        }

        return grid_solution;
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

    Regular_grid_1d<double> solution = hce.get_initial_conditions_grid();

    constexpr int time_steps_between_checks = 100;

    int time_step = 0;
    double t = 0.0;

    while (t <= time_of_compute) {

        solution = hce.get_solution(solution, std::make_unique<RKM>(), t, delta_t);

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
