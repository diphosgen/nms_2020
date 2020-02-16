#include "test_matrix_fixed.hxx"

namespace Tests
{

void test_matrix_fixed() noexcept
{
    std::cout << "|---- test_matrix_fixed()" << std::endl;

    using namespace Math_structures;

    Matrix_fixed<int, 2, 3> m1{{1, 2, 3}, {5, 6, 7}};
    Matrix_fixed<int, 2, 3> m2{{9, 7, 5}, {3, 1, 9}};

    assert(m1.get_cols_amount() == 3);
    assert(m1.get_rows_amount() == 2);

    assert(m2.get_cols_amount() == 3);
    assert(m2.get_rows_amount() == 2);

    assert(!(m1 == m2));

    assert(m1[0][0] == 1);
    assert(m1[0][1] == 2);
    assert(m1[0][2] == 3);
    assert(m1[1][0] == 5);
    assert(m1[1][1] == 6);
    assert(m1[1][2] == 7);

    assert(m2[0][0] == 9);
    assert(m2[0][1] == 7);
    assert(m2[0][2] == 5);
    assert(m2[1][0] == 3);
    assert(m2[1][1] == 1);
    assert(m2[1][2] == 9);

    assert(max_abs(m1) == 7);
    assert(max_abs(m2) == 9);

    m1 = m2;

    assert(m1 == m2);

    assert(m1[0][0] == 9);
    assert(m1[0][1] == 7);
    assert(m1[0][2] == 5);
    assert(m1[1][0] == 3);
    assert(m1[1][1] == 1);
    assert(m1[1][2] == 9);

    std::array<int, 3> v10 = m1[0];
    const std::array<int, 3> v11 = m1[1];

    assert(v10[0] == 9);
    assert(v10[1] == 7);
    assert(v10[2] == 5);

    assert(v11[0] == 3);
    assert(v11[1] == 1);
    assert(v11[2] == 9);
}

}
