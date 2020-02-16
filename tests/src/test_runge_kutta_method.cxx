#include "test_runge_kutta_method.hxx"

namespace Tests
{

void test_runge_kutta_method() noexcept
{
    std::cout << "|---- test_runge_kutta_method()" << std::endl;

    test_forward_Euler_method();
    test_explicit_midpoint_method();
    test_explicit_Heun_second_order_method();
    test_explicit_Ralston_second_order_method();
    test_explicit_generic_second_order_method();
    test_explicit_Kutta_third_order_method();
    test_explicit_generic_third_order_method();
    test_explicit_Heun_third_order_method();
    test_explicit_Ralston_third_order_method();
    test_explicit_SSPRK3_method();
    test_classic_fourth_order_method();
    test_explicit_Ralston_fourth_order_method();
    test_classic_three_eighths_rule_method();
    test_explicit_Fehlberg_second_order_method();
    test_explicit_Bogacki_Shampine_third_order_method();
    test_explicit_Runge_Kutta_Fehlberg_fifth_order_method();
    test_explicit_Cash_Karp_fifth_order_method();
    test_Dormand_Prince_method();
    test_backward_Euler_method();
    test_implicit_midpoint_method();
    test_implicit_Crank_Nicolson_method();
    test_Gauss_Legendre_fourth_order_method();
    test_Gauss_Legendre_sixth_order_method();
    test_Qin_Zhang_symplectic_second_order_DIRK();
    test_Pareschi_Russo_second_order_DIRK();
    test_Crouzeix_third_order_DIRK();
    test_Norsett_fourth_order_DIRK();
    test_Norsett_third_order_DIRK();
    test_Lobatto_IIIA_second_order_method();
    test_Lobatto_IIIA_fourth_order_method();
    test_Lobatto_IIIA_sixth_order_method();
    test_Lobatto_IIIA_eighth_order_method();
    test_Lobatto_IIIB_second_order_method();
    test_Lobatto_IIIB_fourth_order_method();
    test_Lobatto_IIIB_sixth_order_method();
    test_Lobatto_IIIB_eighth_order_method();
    test_Lobatto_IIIC_second_order_method();
    test_Lobatto_IIIC_fourth_order_method();
    test_Lobatto_IIIC_sixth_order_method();
    test_Lobatto_IIIC_eighth_order_method();
    test_Lobatto_IIINW_second_order_method();
    test_Lobatto_IIINW_fourth_order_method();
    test_Radau_IA_third_order_method();
    test_Radau_IA_fifth_order_method();
    test_Radau_IIA_third_order_method();
    test_Radau_IIA_fifth_order_method();
    test_Ceschino_Kunzmann_third_order_DIRK();
    test_Ceschino_Kunzmann_fifth_order_DIRK();
}

}
