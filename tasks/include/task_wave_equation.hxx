#ifndef TASK_WAVE_EQUATION_HXX_INCLUDED
#define TASK_WAVE_EQUATION_HXX_INCLUDED

#include "basic_libs.hxx"
#include "basic_math_functions.hxx"
#include "wave_equation_1d.hxx"
#include "wave_equation_2d.hxx"
#include "runge_kutta_method_lib.hxx"
#include "richardson_extrapolation.hxx"

namespace Tasks
{

void wave_equation_1D_solution() noexcept;
void wave_equation_2D_solution() noexcept;

void task_wave_equation() noexcept;

}

#endif // TASK_WAVE_EQUATION_HXX_INCLUDED
