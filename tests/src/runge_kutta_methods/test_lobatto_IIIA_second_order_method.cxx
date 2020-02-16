#include "test_runge_kutta_method.hxx"

namespace Tests
{

void test_Lobatto_IIIA_second_order_method() noexcept
{
    std::cout << "      |---- test_Lobatto_IIIA_second_order_method()" << std::endl;

    using namespace Math_structures;

    constexpr int equations_amount = 2;
    constexpr int frequencies_amount = 10;

    using T_arg_t = double;
    using T_arg_u = Vector_unfixed<Vector_fixed<T_arg_t, equations_amount>>;

    using ODE_solver = One_step_ODE_solver<T_arg_u, T_arg_t>;
    using RKM = Lobatto_IIIA_second_order_method<T_arg_u, T_arg_t>;

    constexpr T_arg_t min_frequency = 0.01;
    constexpr T_arg_t max_frequency = 2.00;

    std::mt19937 gen_rand(0);
    std::uniform_real_distribution<double> distr_frequency(min_frequency,
                                                           max_frequency);

    auto rand_frequency = [&gen_rand, &distr_frequency]() noexcept -> double
    {
        return distr_frequency(gen_rand);
    };

    std::array<T_arg_t, frequencies_amount> frequencies{};

    for (auto& frequency: frequencies) {
        frequency = rand_frequency();
    }

    auto f = [frequencies_amount, &frequencies](const T_arg_u& u, const T_arg_t& t)
        noexcept -> T_arg_u
    {
        T_arg_u f_rhs(frequencies_amount, {{0.0, 0.0}});

        for (int i = 0; i < frequencies_amount; ++i) {
            const T_arg_t frequency = frequencies[i];
            const T_arg_t frequency_sqr = frequency * frequency;
            f_rhs[i] = {{u[i][1], (-frequency_sqr) * u[i][0]}};
        }

        return f_rhs;
    };

    auto exact_solution = [frequencies_amount, &frequencies](const T_arg_t& t)
        noexcept -> T_arg_u
    {
        T_arg_u res(frequencies_amount, {{0.0, 0.0}});

        for (int i = 0; i < frequencies_amount; ++i) {
            const T_arg_t frequency = frequencies[i];
            res[i] = {{
                    std::cos(frequency * t),
                    std::sin(frequency * t) * (-frequency)
                }};
        }

        return res;
    };

    auto initial_conditions = [frequencies_amount, &frequencies]()
        noexcept -> T_arg_u
    {
        return T_arg_u(frequencies_amount, {{1.0, 0.0}});
    };

    constexpr T_arg_t max_error = 1.0e-4;

    constexpr T_arg_t t_end = 100.0;
    constexpr int steps_amount = 1'000'000;
    constexpr T_arg_t dt = t_end / steps_amount;

    constexpr int time_steps_between_checks = 100;

    T_arg_u solution = initial_conditions();

    std::unique_ptr<ODE_solver> solver = std::make_unique<RKM>();

    for (int time_step = 0; time_step <= steps_amount; ++time_step) {

        const T_arg_t t = dt * time_step;
        solution = solver->get_solution(f, solution, t, dt);

        if (time_step % time_steps_between_checks == 0) {
            const auto res_exact_solution = exact_solution(t + dt);
            const double eps = 1.0e-12;
            const double abs_norm = max_abs(solution + res_exact_solution) + eps;
            const double error = max_abs(solution - res_exact_solution)/abs_norm;
            assert(error < max_error);
        }

    }
}

}
