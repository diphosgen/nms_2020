#include "tasks.hxx"

namespace Tasks
{

void run_task(const int task_number) noexcept
{
    switch(task_number) {
        case 1  : { task_ode_system_solution(); }       break;
        case 2  : { task_heat_conduction_equation(); }  break;
        case 3  : { task_wave_equation(); }             break;
        case 4  : { task_poisson_equation(); }          break;
        default : { std::cout << "task_number error" << std::endl; }
    }
}

}
