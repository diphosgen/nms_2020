#ifndef TASK_ODE_SYSTEM_SOLUTION_HXX_INCLUDED
#define TASK_ODE_SYSTEM_SOLUTION_HXX_INCLUDED

#include "basic_libs.hxx"
#include "basic_math_functions.hxx"
#include "vector_fixed.hxx"
#include "runge_kutta_method_lib.hxx"
#include "richardson_extrapolation.hxx"

namespace Tasks
{

void harmonic_oscillator_solution() noexcept;
void individual_ode_system_solution() noexcept;

void task_ode_system_solution() noexcept;

}

#endif // TASK_ODE_SYSTEM_SOLUTION_HXX_INCLUDED
