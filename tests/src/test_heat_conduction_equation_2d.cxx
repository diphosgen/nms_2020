#include "test_heat_conduction_equation_2d.hxx"

namespace Tests
{

void test_heat_conduction_equation_2D() noexcept
{
    std::cout << "|---- test_heat_conduction_equation_2D()" << std::endl;

    test_heat_conduction_equation_2D_1();
    test_heat_conduction_equation_2D_2();
}

}
