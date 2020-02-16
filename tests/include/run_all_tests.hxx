#ifndef RUN_ALL_TESTS_HXX_INCLUDED
#define RUN_ALL_TESTS_HXX_INCLUDED

#include "test_vector_unfixed.hxx"
#include "test_vector_fixed.hxx"
#include "test_matrix_unfixed.hxx"
#include "test_matrix_fixed.hxx"
#include "test_slae_solver_gauss.hxx"
#include "test_snae_solver_fpi.hxx"
#include "test_butcher_table.hxx"
#include "test_runge_kutta_method.hxx"
#include "test_richardson_extrapolation.hxx"
#include "test_regular_grid_1d.hxx"
#include "test_regular_grid_2d.hxx"
#include "test_regular_grid_3d.hxx"
#include "test_heat_conduction_equation_1d.hxx"
#include "test_heat_conduction_equation_2d.hxx"
#include "test_wave_equation_1d.hxx"
#include "test_wave_equation_2d.hxx"

namespace Tests
{

void run_all_tests() noexcept;

}

#endif // RUN_ALL_TESTS_HXX_INCLUDED
