#ifndef TEST_RUNGE_KUTTA_METHOD_HXX_INCLUDED
#define TEST_RUNGE_KUTTA_METHOD_HXX_INCLUDED

#include "one_step_ode_solver.hxx"
#include "runge_kutta_method.hxx"
#include "runge_kutta_method_expl.hxx"
#include "runge_kutta_method_impl.hxx"
#include "runge_kutta_method_lib.hxx"

#include "vector_unfixed.hxx"
#include "vector_fixed.hxx"
#include "matrix_unfixed.hxx"
#include "matrix_fixed.hxx"

namespace Tests
{

void test_runge_kutta_method() noexcept;

void test_forward_Euler_method() noexcept;
void test_explicit_midpoint_method() noexcept;
void test_explicit_Heun_second_order_method() noexcept;
void test_explicit_Ralston_second_order_method() noexcept;
void test_explicit_generic_second_order_method() noexcept;
void test_explicit_Kutta_third_order_method() noexcept;
void test_explicit_generic_third_order_method() noexcept;
void test_explicit_Heun_third_order_method() noexcept;
void test_explicit_Ralston_third_order_method() noexcept;
void test_explicit_SSPRK3_method() noexcept;
void test_classic_fourth_order_method() noexcept;
void test_explicit_Ralston_fourth_order_method() noexcept;
void test_classic_three_eighths_rule_method() noexcept;
void test_explicit_Fehlberg_second_order_method() noexcept;
void test_explicit_Bogacki_Shampine_third_order_method() noexcept;
void test_explicit_Runge_Kutta_Fehlberg_fifth_order_method() noexcept;
void test_explicit_Cash_Karp_fifth_order_method() noexcept;
void test_Dormand_Prince_method() noexcept;
void test_backward_Euler_method() noexcept;
void test_implicit_midpoint_method() noexcept;
void test_implicit_Crank_Nicolson_method() noexcept;
void test_Gauss_Legendre_fourth_order_method() noexcept;
void test_Gauss_Legendre_sixth_order_method() noexcept;
void test_Qin_Zhang_symplectic_second_order_DIRK() noexcept;
void test_Pareschi_Russo_second_order_DIRK() noexcept;
void test_Crouzeix_third_order_DIRK() noexcept;
void test_Norsett_fourth_order_DIRK() noexcept;
void test_Norsett_third_order_DIRK() noexcept;
void test_Lobatto_IIIA_second_order_method() noexcept;
void test_Lobatto_IIIA_fourth_order_method() noexcept;
void test_Lobatto_IIIA_sixth_order_method() noexcept;
void test_Lobatto_IIIA_eighth_order_method() noexcept;
void test_Lobatto_IIIB_second_order_method() noexcept;
void test_Lobatto_IIIB_fourth_order_method() noexcept;
void test_Lobatto_IIIB_sixth_order_method() noexcept;
void test_Lobatto_IIIB_eighth_order_method() noexcept;
void test_Lobatto_IIIC_second_order_method() noexcept;
void test_Lobatto_IIIC_fourth_order_method() noexcept;
void test_Lobatto_IIIC_sixth_order_method() noexcept;
void test_Lobatto_IIIC_eighth_order_method() noexcept;
void test_Lobatto_IIINW_second_order_method() noexcept;
void test_Lobatto_IIINW_fourth_order_method() noexcept;
void test_Radau_IA_third_order_method() noexcept;
void test_Radau_IA_fifth_order_method() noexcept;
void test_Radau_IIA_third_order_method() noexcept;
void test_Radau_IIA_fifth_order_method() noexcept;
void test_Ceschino_Kunzmann_third_order_DIRK() noexcept;
void test_Ceschino_Kunzmann_fifth_order_DIRK() noexcept;

}

#endif // TEST_RUNGE_KUTTA_METHOD_HXX_INCLUDED
