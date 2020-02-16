#include "test_butcher_table.hxx"

namespace Tests
{

static void test_butcher_table_1() noexcept;
static void test_butcher_table_2() noexcept;

void test_butcher_table() noexcept
{
    std::cout << "|---- test_butcher_table()" << std::endl;

    test_butcher_table_1();
    test_butcher_table_2();
}

static void test_butcher_table_1() noexcept
{
    std::cout << "      |---- test_butcher_table_1()" << std::endl;

    using namespace Math_structures;

    constexpr int stages_amount = 4;

    const Butcher_table<stages_amount> b_n = {
        {{
            {{0.0, 0.0, 0.0, 0.0}},
            {{0.5, 0.0, 0.0, 0.0}},
            {{0.0, 0.5, 0.0, 0.0}},
            {{0.0, 0.0, 1.0, 0.0}}
        }},
        {{0.0, 0.5, 0.5, 1.0}},
        {{1.0/6.0, 2.0/6.0, 2.0/6.0, 1.0/6.0}}
    };

    assert(compare_floating_numbers(b_n.get_element_a_table(0,0), 0.0));
    assert(compare_floating_numbers(b_n.get_element_a_table(0,1), 0.0));
    assert(compare_floating_numbers(b_n.get_element_a_table(0,2), 0.0));
    assert(compare_floating_numbers(b_n.get_element_a_table(0,3), 0.0));

    assert(compare_floating_numbers(b_n.get_element_a_table(1,0), 0.5));
    assert(compare_floating_numbers(b_n.get_element_a_table(1,1), 0.0));
    assert(compare_floating_numbers(b_n.get_element_a_table(1,2), 0.0));
    assert(compare_floating_numbers(b_n.get_element_a_table(1,3), 0.0));

    assert(compare_floating_numbers(b_n.get_element_a_table(2,0), 0.0));
    assert(compare_floating_numbers(b_n.get_element_a_table(2,1), 0.5));
    assert(compare_floating_numbers(b_n.get_element_a_table(2,2), 0.0));
    assert(compare_floating_numbers(b_n.get_element_a_table(2,3), 0.0));

    assert(compare_floating_numbers(b_n.get_element_a_table(3,0), 0.0));
    assert(compare_floating_numbers(b_n.get_element_a_table(3,1), 0.0));
    assert(compare_floating_numbers(b_n.get_element_a_table(3,2), 1.0));
    assert(compare_floating_numbers(b_n.get_element_a_table(3,3), 0.0));

    assert(compare_floating_numbers(b_n.get_element_c_table(0), 0.0));
    assert(compare_floating_numbers(b_n.get_element_c_table(1), 0.5));
    assert(compare_floating_numbers(b_n.get_element_c_table(2), 0.5));
    assert(compare_floating_numbers(b_n.get_element_c_table(3), 1.0));

    assert(compare_floating_numbers(b_n.get_element_b_table(0), 1.0/6.0));
    assert(compare_floating_numbers(b_n.get_element_b_table(1), 2.0/6.0));
    assert(compare_floating_numbers(b_n.get_element_b_table(2), 2.0/6.0));
    assert(compare_floating_numbers(b_n.get_element_b_table(3), 1.0/6.0));
}

static void test_butcher_table_2() noexcept
{
    std::cout << "      |---- test_butcher_table_2()" << std::endl;

    using namespace Math_structures;

    constexpr int stages_amount = 4;

    const Butcher_table<stages_amount, float> b_n = {
        {{
            {{0.0, 0.0, 0.0, 0.0}},
            {{0.5, 0.0, 0.0, 0.0}},
            {{0.0, 0.5, 0.0, 0.0}},
            {{0.0, 0.0, 1.0, 0.0}}
        }},
        {{0.0, 0.5, 0.5, 1.0}},
        {{1.0/6.0, 2.0/6.0, 2.0/6.0, 1.0/6.0}}
    };

    assert(compare_floating_numbers(b_n.get_element_a_table(0,0), 0.0));
    assert(compare_floating_numbers(b_n.get_element_a_table(0,1), 0.0));
    assert(compare_floating_numbers(b_n.get_element_a_table(0,2), 0.0));
    assert(compare_floating_numbers(b_n.get_element_a_table(0,3), 0.0));

    assert(compare_floating_numbers(b_n.get_element_a_table(1,0), 0.5));
    assert(compare_floating_numbers(b_n.get_element_a_table(1,1), 0.0));
    assert(compare_floating_numbers(b_n.get_element_a_table(1,2), 0.0));
    assert(compare_floating_numbers(b_n.get_element_a_table(1,3), 0.0));

    assert(compare_floating_numbers(b_n.get_element_a_table(2,0), 0.0));
    assert(compare_floating_numbers(b_n.get_element_a_table(2,1), 0.5));
    assert(compare_floating_numbers(b_n.get_element_a_table(2,2), 0.0));
    assert(compare_floating_numbers(b_n.get_element_a_table(2,3), 0.0));

    assert(compare_floating_numbers(b_n.get_element_a_table(3,0), 0.0));
    assert(compare_floating_numbers(b_n.get_element_a_table(3,1), 0.0));
    assert(compare_floating_numbers(b_n.get_element_a_table(3,2), 1.0));
    assert(compare_floating_numbers(b_n.get_element_a_table(3,3), 0.0));

    assert(compare_floating_numbers(b_n.get_element_c_table(0), 0.0));
    assert(compare_floating_numbers(b_n.get_element_c_table(1), 0.5));
    assert(compare_floating_numbers(b_n.get_element_c_table(2), 0.5));
    assert(compare_floating_numbers(b_n.get_element_c_table(3), 1.0));

    assert(compare_floating_numbers(b_n.get_element_b_table(0), 1.0/6.0));
    assert(compare_floating_numbers(b_n.get_element_b_table(1), 2.0/6.0));
    assert(compare_floating_numbers(b_n.get_element_b_table(2), 2.0/6.0));
    assert(compare_floating_numbers(b_n.get_element_b_table(3), 1.0/6.0));
}

}
