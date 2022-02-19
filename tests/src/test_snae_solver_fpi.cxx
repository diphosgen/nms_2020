#include "test_snae_solver_fpi.hxx"

namespace Tests
{

static void test_SNAE_solver_FPI_1() noexcept;
static void test_SNAE_solver_FPI_2() noexcept;
static void test_SNAE_solver_FPI_3() noexcept;
static void test_SNAE_solver_FPI_4() noexcept;

void test_SNAE_solver_FPI() noexcept
{
    std::cout << "|---- test_SNAE_solver_FPI()" << std::endl;

    test_SNAE_solver_FPI_1();
    test_SNAE_solver_FPI_2();
    test_SNAE_solver_FPI_3();
    test_SNAE_solver_FPI_4();
}

static void test_SNAE_solver_FPI_1() noexcept
{
    std::cout << "      |---- test_SNAE_solver_FPI_1()" << std::endl;

    using namespace Math_structures;

    auto snae_func = [](const double& x)  -> double {
        return sin(x) - 0.5;
    };

    const double solution = SNAE_solver_FPI<double>{}.get_solution(snae_func, 0.0);
    constexpr double error_max = 1.0e-5;
    assert(max_abs(snae_func(solution)) < error_max);
}

static void test_SNAE_solver_FPI_2() noexcept
{
    std::cout << "      |---- test_SNAE_solver_FPI_2()" << std::endl;

    using namespace Math_structures;

    auto snae_func = [](const Vector_unfixed<double>& x) -> Vector_unfixed<double> {
        return Vector_unfixed<double>{cos(x[0])};
    };

    const Vector_unfixed<double> solution = SNAE_solver_FPI<Vector_unfixed<double>>{}.
                                                get_solution(snae_func, {1.0});

    constexpr double error_max = 1.0e-5;
    assert(max_abs(snae_func(solution)) < error_max);
}

static void test_SNAE_solver_FPI_3() noexcept
{
    std::cout << "      |---- test_SNAE_solver_FPI_3()" << std::endl;

    using namespace Math_structures;

    auto snae_func = [](const Vector_unfixed<Vector_unfixed<double>>& x)
        -> Vector_unfixed<Vector_unfixed<double>> {
        return {{cos(x[0][0])}, {sin(x[1][0])},
                {cos(x[2][0])}, {sin(x[3][0])}};
    };

    const Vector_unfixed<Vector_unfixed<double>> solution =
            SNAE_solver_FPI<Vector_unfixed<Vector_unfixed<double>>>{}.
                get_solution(snae_func, {{1.0}, {0.1}, {0.9}, {1.99}});

    constexpr double error_max = 1.0e-5;
    assert(max_abs(snae_func(solution)) < error_max);
}

static void test_SNAE_solver_FPI_4() noexcept
{
    std::cout << "      |---- test_SNAE_solver_FPI_4()" << std::endl;

    using namespace Math_structures;

    auto snae_func = [](const Matrix_unfixed<double>& x)
        -> Matrix_unfixed<double> {
        return {{cos(x[0][0])},
                {sin(x[1][0])},
                {cos(x[2][0])},
                {sin(x[3][0])}};
    };

    const Matrix_unfixed<double> solution =
            SNAE_solver_FPI<Matrix_unfixed<double>>{}.
                get_solution(snae_func, {{1.0}, {0.1}, {0.9}, {1.99}});

    constexpr double error_max = 1.0e-5;
    assert(max_abs(snae_func(solution)) < error_max);
}

}
