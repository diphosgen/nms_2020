#ifndef TASK_HEAT_CONDUCTION_EQUATION_HXX_INCLUDED
#define TASK_HEAT_CONDUCTION_EQUATION_HXX_INCLUDED

#include "basic_libs.hxx"
#include "basic_math_functions.hxx"
#include "heat_conduction_equation_1d.hxx"
#include "heat_conduction_equation_2d.hxx"
#include "runge_kutta_method_lib.hxx"
#include "richardson_extrapolation.hxx"

namespace Tasks
{

void heat_conduction_equation_1D_solution() noexcept;
void heat_conduction_equation_2D_solution() noexcept;

void task_heat_conduction_equation() noexcept;

}

#endif // TASK_HEAT_CONDUCTION_EQUATION_HXX_INCLUDED
