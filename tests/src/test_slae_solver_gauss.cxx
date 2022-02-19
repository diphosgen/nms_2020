#include "test_slae_solver_gauss.hxx"

namespace Tests
{

static void test_SLAE_solver_Gauss_1() noexcept;
static void test_SLAE_solver_Gauss_2() noexcept;
static void test_SLAE_solver_Gauss_3() noexcept;
static void test_SLAE_solver_Gauss_4() noexcept;
static void test_SLAE_solver_Gauss_5() noexcept;
static void test_SLAE_solver_Gauss_6() noexcept;
static void test_SLAE_solver_Gauss_7() noexcept;

void test_SLAE_solver_Gauss() noexcept
{
    std::cout << "|---- test_SLAE_solver_Gauss()" << std::endl;

    test_SLAE_solver_Gauss_1();
    test_SLAE_solver_Gauss_2();
    test_SLAE_solver_Gauss_3();
    test_SLAE_solver_Gauss_4();
    test_SLAE_solver_Gauss_5();
    test_SLAE_solver_Gauss_6();
    test_SLAE_solver_Gauss_7();
}

static void test_SLAE_solver_Gauss_1() noexcept
{
    std::cout << "      |---- test_SLAE_solver_Gauss_1()" << std::endl;

    using namespace Math_structures;

    const int system_size = 200;

    Matrix_unfixed<double> m(system_size, system_size);
    Vector_unfixed<double> v_exact(system_size);

    std::mt19937 gen_rand{0};
    std::uniform_real_distribution<double> distr(-1.0e+5, +1.0e+5);

    for (int i = 0; i < system_size; ++i) {
        for (int j = 0; j < system_size; ++j) {
            if (i == j) {
                m[i][j] = 0.0;
            } else {
                m[i][j] = std::pow(i - j, 5) * distr(gen_rand);
            }
        }
        v_exact[i] = distr(gen_rand);
    }

    Vector_unfixed<double> v_rhs = m * v_exact;

    std::unique_ptr<SLAE_solver<Matrix_unfixed<double>, Vector_unfixed<double>>> solver =
        std::make_unique<SLAE_solver_Gauss<Matrix_unfixed<double>, Vector_unfixed<double>>>();

    Vector_unfixed<double> v = solver->get_solution(m, v_rhs);

    Vector_unfixed<double> v_error = v_exact - v;

    constexpr double error_max = 1.0e-5;
    assert(max_abs(v_exact - v) < error_max);
}

static void test_SLAE_solver_Gauss_2() noexcept
{
    std::cout << "      |---- test_SLAE_solver_Gauss_2()" << std::endl;

    using namespace Math_structures;

    const int system_size = 200;

    Matrix_unfixed<std::complex<double>> m(system_size, system_size);
    Vector_unfixed<std::complex<double>> v_exact(system_size);

    std::mt19937 gen_rand{0};
    std::uniform_real_distribution<double> distr(-1.0e+5, +1.0e+5);

    for (int i = 0; i < system_size; ++i) {
        for (int j = 0; j < system_size; ++j) {
            if (i == j) {
                m[i][j] = {0.0, 0.0};
            } else {
                m[i][j] = {std::pow(i - j, 5) * distr(gen_rand),
                           std::pow(i - j, 5) * distr(gen_rand)};
            }
        }
        v_exact[i] = {distr(gen_rand), distr(gen_rand)};
    }

    Vector_unfixed<std::complex<double>> v_rhs = m * v_exact;

    std::unique_ptr<SLAE_solver<Matrix_unfixed<std::complex<double>>,
                                Vector_unfixed<std::complex<double>>>> solver =
        std::make_unique<SLAE_solver_Gauss<Matrix_unfixed<std::complex<double>>,
                                           Vector_unfixed<std::complex<double>>>>();

    Vector_unfixed<std::complex<double>> v = solver->get_solution(m, v_rhs);

    Vector_unfixed<std::complex<double>> v_error = v_exact - v;

    constexpr double error_max = 1.0e-5;
    assert(max_abs(v_exact - v) < error_max);
}

static void test_SLAE_solver_Gauss_3() noexcept
{
    std::cout << "      |---- test_SLAE_solver_Gauss_3()" << std::endl;

    using namespace Math_structures;

    const int system_size = 50;
    const int vector_size = 20'000;

    Matrix_unfixed<double> m(system_size, system_size);
    Vector_unfixed<Vector_unfixed<double>> v_exact(system_size);

    std::mt19937 gen_rand{0};
    std::uniform_real_distribution<double> distr(-1.0e+5, +1.0e+5);

    for (int i = 0; i < system_size; ++i) {
        for (int j = 0; j < system_size; ++j) {
            if (i == j) {
                m[i][j] = 0.0;
            } else {
                m[i][j] = std::pow(i - j, 5) * distr(gen_rand);
            }
        }

        v_exact[i] = Vector_unfixed<double>(vector_size);
        for (int j = 0; j < vector_size; ++j) {
            v_exact[i][j] = distr(gen_rand);
        }
    }

    Vector_unfixed<Vector_unfixed<double>> v_rhs = m * v_exact;

    std::unique_ptr<SLAE_solver<Matrix_unfixed<double>, Vector_unfixed<Vector_unfixed<double>>>> solver =
        std::make_unique<SLAE_solver_Gauss<Matrix_unfixed<double>, Vector_unfixed<Vector_unfixed<double>>>>();

    Vector_unfixed<Vector_unfixed<double>> v = solver->get_solution(m, v_rhs);

    Vector_unfixed<Vector_unfixed<double>> v_error = v_exact - v;

    constexpr double error_max = 1.0e-5;
    assert(max_abs(v_exact - v) < error_max);
}

static void test_SLAE_solver_Gauss_4() noexcept
{
    std::cout << "      |---- test_SLAE_solver_Gauss_4()" << std::endl;

    using namespace Math_structures;

    const int system_size = 50;
    const int matrix_size = 100;

    Matrix_unfixed<double> m(system_size, system_size);
    Vector_unfixed<Matrix_unfixed<double>> v_exact(system_size,
                                                   Matrix_unfixed<double>(matrix_size,
                                                                          matrix_size));

    std::mt19937 gen_rand{0};
    std::uniform_real_distribution<double> distr(-1.0e+5, +1.0e+5);

    for (int i = 0; i < system_size; ++i) {
        for (int j = 0; j < system_size; ++j) {
            if (i == j) {
                m[i][j] = 0.0;
            } else {
                m[i][j] = std::pow(i - j, 5) * distr(gen_rand);
            }
        }

        v_exact[i] = Matrix_unfixed<double>(matrix_size, matrix_size);
        for (int j = 0; j < matrix_size; ++j) {
            for (int k = 0; k < matrix_size; ++k) {
                v_exact[i][j][k] = distr(gen_rand);
            }
        }
    }

    Vector_unfixed<Matrix_unfixed<double>> v_rhs = m * v_exact;

    std::unique_ptr<SLAE_solver<Matrix_unfixed<double>,
                                Vector_unfixed<Matrix_unfixed<double>>>> solver =
        std::make_unique<SLAE_solver_Gauss<Matrix_unfixed<double>,
                                           Vector_unfixed<Matrix_unfixed<double>>>>();

    Vector_unfixed<Matrix_unfixed<double>> v = solver->get_solution(m, v_rhs);

    Vector_unfixed<Matrix_unfixed<double>> v_error = v_exact - v;

    constexpr double error_max = 1.0e-5;
    assert(max_abs(v_exact - v) < error_max);
}

static void test_SLAE_solver_Gauss_5() noexcept
{
    std::cout << "      |---- test_SLAE_solver_Gauss_5()" << std::endl;

    using namespace Math_structures;

    constexpr int system_size = 50;
    constexpr int matrix_size = 10;

    Matrix_fixed<double, system_size, system_size> m;
    Vector_fixed<Matrix_fixed<double, matrix_size, matrix_size>, system_size> v_exact{};

    std::mt19937 gen_rand{0};
    std::uniform_real_distribution<double> distr(-1.0e+5, +1.0e+5);

    for (int i = 0; i < system_size; ++i) {
        for (int j = 0; j < system_size; ++j) {
            if (i == j) {
                m[i][j] = 0.0;
            } else {
                m[i][j] = std::pow(i - j, 5) * distr(gen_rand);
            }
        }

        for (int j = 0; j < matrix_size; ++j) {
            for (int k = 0; k < matrix_size; ++k) {
                v_exact[i][j][k] = distr(gen_rand);
            }
        }
    }

    Vector_fixed<Matrix_fixed<double, matrix_size, matrix_size>, system_size> v_rhs = m * v_exact;

    std::unique_ptr<SLAE_solver<Matrix_fixed<double, system_size, system_size>,
                                Vector_fixed<Matrix_fixed<double, matrix_size, matrix_size>,
                                system_size>>> solver =
        std::make_unique<SLAE_solver_Gauss<Matrix_fixed<double, system_size, system_size>,
                                           Vector_fixed<Matrix_fixed<double, matrix_size, matrix_size>,
                                           system_size>>>();

    Vector_fixed<Matrix_fixed<double, matrix_size, matrix_size>, system_size> v = solver->get_solution(m, v_rhs);

    constexpr double error_max = 1.0e-5;
    assert(max_abs(v_exact - v) < error_max);
}

static void test_SLAE_solver_Gauss_6() noexcept
{
    std::cout << "      |---- test_SLAE_solver_Gauss_6()" << std::endl;

    using namespace Math_structures;

    constexpr int system_size = 50;
    constexpr int matrix_size = 10;
    constexpr int matrix_internal_size = 10;

    Matrix_fixed<double, system_size, system_size> m;
    Vector_fixed<Matrix_fixed<Matrix_unfixed<double>, matrix_size, matrix_size>, system_size> v_exact{};

    std::mt19937 gen_rand{0};
    std::uniform_real_distribution<double> distr(-1.0e+5, +1.0e+5);

    for (int i = 0; i < system_size; ++i) {
        for (int j = 0; j < system_size; ++j) {
            if (i == j) {
                m[i][j] = 0.0;
            } else {
                m[i][j] = std::pow(i - j, 5) * distr(gen_rand);
            }
        }

        for (int j = 0; j < matrix_size; ++j) {
            for (int k = 0; k < matrix_size; ++k) {
                Matrix_unfixed<double> m_temp(matrix_internal_size, matrix_internal_size);

                for (int _i = 0; _i < matrix_internal_size; ++_i) {
                    for (int _j = 0; _j < matrix_internal_size; ++_j) {
                        m_temp[_i][_j] = distr(gen_rand);
                    }
                }

                v_exact[i][j][k] = m_temp;
            }
        }
    }

    Vector_fixed<Matrix_fixed<Matrix_unfixed<double>, matrix_size, matrix_size>, system_size> v_rhs = m * v_exact;

    std::unique_ptr<SLAE_solver<Matrix_fixed<double, system_size, system_size>,
                                Vector_fixed<Matrix_fixed<Matrix_unfixed<double>, matrix_size, matrix_size>,
                                system_size>>> solver =
        std::make_unique<SLAE_solver_Gauss<Matrix_fixed<double, system_size, system_size>,
                                           Vector_fixed<Matrix_fixed<Matrix_unfixed<double>, matrix_size, matrix_size>,
                                           system_size>>>();

    Vector_fixed<Matrix_fixed<Matrix_unfixed<double>,
                              matrix_size, matrix_size>,
                              system_size> v = solver->get_solution(m, v_rhs);

    Vector_fixed<Matrix_fixed<Matrix_unfixed<double>, matrix_size, matrix_size>, system_size> v_error = v_exact - v;

    constexpr double error_max = 1.0e-5;
    assert(max_abs(v_exact - v) < error_max);
}

static void test_SLAE_solver_Gauss_7() noexcept
{
    std::cout << "      |---- test_SLAE_solver_Gauss_7()" << std::endl;

    using namespace Math_structures;

    constexpr int system_size = 50;
    constexpr int matrix_size = 10;
    constexpr int matrix_internal_size = 10;

    Matrix_fixed<std::complex<double>, system_size, system_size> m;
    Vector_fixed<Matrix_fixed<Matrix_unfixed<std::complex<double>>, matrix_size, matrix_size>, system_size> v_exact{};

    std::mt19937 gen_rand{0};
    std::uniform_real_distribution<double> distr(-1.0e+5, +1.0e+5);

    for (int i = 0; i < system_size; ++i) {
        for (int j = 0; j < system_size; ++j) {
            if (i == j) {
                m[i][j] = {0.0, 0.0};
            } else {
                m[i][j] = std::pow(i - j, 5) * std::complex<double>{distr(gen_rand), distr(gen_rand)};
            }
        }

        for (int j = 0; j < matrix_size; ++j) {
            for (int k = 0; k < matrix_size; ++k) {
                Matrix_unfixed<std::complex<double>> m_temp(matrix_internal_size, matrix_internal_size);

                for (int _i = 0; _i < matrix_internal_size; ++_i) {
                    for (int _j = 0; _j < matrix_internal_size; ++_j) {
                        m_temp[_i][_j] = {distr(gen_rand), distr(gen_rand)};
                    }
                }

                v_exact[i][j][k] = m_temp;
            }
        }
    }

    Vector_fixed<Matrix_fixed<Matrix_unfixed<std::complex<double>>,
                              matrix_size, matrix_size>,
                 system_size> v_rhs = m * v_exact;

    std::unique_ptr<SLAE_solver<Matrix_fixed<std::complex<double>,
                                             system_size, system_size>,
                                Vector_fixed<Matrix_fixed<Matrix_unfixed<std::complex<double>>,
                                                          matrix_size, matrix_size>,
                                             system_size>>> solver =
        std::make_unique<SLAE_solver_Gauss<Matrix_fixed<std::complex<double>,
                                                        system_size, system_size>,
                                           Vector_fixed<Matrix_fixed<Matrix_unfixed<std::complex<double>>,
                                                                     matrix_size, matrix_size>,
                                                        system_size>>>();

    Vector_fixed<Matrix_fixed<Matrix_unfixed<std::complex<double>>,
                              matrix_size, matrix_size>,
                 system_size> v = solver->get_solution(m, v_rhs);

    Vector_fixed<Matrix_fixed<Matrix_unfixed<std::complex<double>>,
                              matrix_size, matrix_size>,
                 system_size> v_error = v_exact - v;

    constexpr double error_max = 1.0e-5;
    assert(max_abs(v_exact - v) < error_max);
}

}
