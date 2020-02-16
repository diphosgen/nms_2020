#ifndef TASKS_HXX_INCLUDED
#define TASKS_HXX_INCLUDED

#include "task_ode_system_solution.hxx"
#include "task_heat_conduction_equation.hxx"
#include "task_wave_equation.hxx"
#include "task_poisson_equation.hxx"

namespace Tasks
{

void run_task(const int task_number) noexcept;

}

#endif // TASKS_HXX_INCLUDED
