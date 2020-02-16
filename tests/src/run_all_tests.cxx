#include "run_all_tests.hxx"

namespace Tests
{

void run_all_tests() noexcept
{
    std::cout << "run_all_tests()" << std::endl;

    test_vector_unfixed();
    test_vector_fixed();
    test_matrix_unfixed();
    test_matrix_fixed();
    test_SLAE_solver_Gauss();
    test_SNAE_solver_FPI();
    test_butcher_table();
    test_runge_kutta_method();
    test_richardson_extrapolation();
    test_regular_grid_1d();
    test_regular_grid_2d();
    test_regular_grid_3d();
    test_heat_conduction_equation_1D();
    test_heat_conduction_equation_2D();
    test_wave_equation_1D();
    test_wave_equation_2D();
}

}
