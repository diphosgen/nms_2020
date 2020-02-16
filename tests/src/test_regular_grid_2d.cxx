#include "test_regular_grid_2d.hxx"

namespace Tests
{

void test_regular_grid_2d() noexcept
{
    std::cout << "|---- test_regular_grid_2d()" << std::endl;

    using namespace Math_structures;

    Regular_grid_2d<int> rg_2d_1{{1, 2, 3}, {5, 6, 7}};
    Regular_grid_2d<int> rg_2d_2{{9, 7, 5}, {3, 1, 9}};

    assert(rg_2d_1.get_size_x() == 2);
    assert(rg_2d_1.get_size_y() == 3);

    assert(rg_2d_2.get_size_x() == 2);
    assert(rg_2d_2.get_size_y() == 3);

    assert(!(rg_2d_1 == rg_2d_2));

    assert(rg_2d_1[0][0] == 1);
    assert(rg_2d_1[0][1] == 2);
    assert(rg_2d_1[0][2] == 3);
    assert(rg_2d_1[1][0] == 5);
    assert(rg_2d_1[1][1] == 6);
    assert(rg_2d_1[1][2] == 7);

    assert(rg_2d_2[0][0] == 9);
    assert(rg_2d_2[0][1] == 7);
    assert(rg_2d_2[0][2] == 5);
    assert(rg_2d_2[1][0] == 3);
    assert(rg_2d_2[1][1] == 1);
    assert(rg_2d_2[1][2] == 9);

    assert(max_abs(rg_2d_1) == 7);
    assert(max_abs(rg_2d_2) == 9);

    rg_2d_1 = rg_2d_2;

    assert(rg_2d_1 == rg_2d_2);

    assert(rg_2d_1[0][0] == 9);
    assert(rg_2d_1[0][1] == 7);
    assert(rg_2d_1[0][2] == 5);
    assert(rg_2d_1[1][0] == 3);
    assert(rg_2d_1[1][1] == 1);
    assert(rg_2d_1[1][2] == 9);
}

}
