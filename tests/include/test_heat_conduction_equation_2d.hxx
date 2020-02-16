#ifndef TEST_HEAT_CONDUCTION_EQUATION_2D_HXX_INCLUDED
#define TEST_HEAT_CONDUCTION_EQUATION_2D_HXX_INCLUDED

#include "heat_conduction_equation_2d.hxx"
#include "runge_kutta_method_lib.hxx"
#include "richardson_extrapolation.hxx"

namespace Tests
{

void test_heat_conduction_equation_2D() noexcept;

void test_heat_conduction_equation_2D_1() noexcept;
void test_heat_conduction_equation_2D_2() noexcept;

}

#endif // TEST_HEAT_CONDUCTION_EQUATION_2D_HXX_INCLUDED
