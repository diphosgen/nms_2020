#include "test_regular_grid_1d.hxx"

namespace Tests
{

void test_regular_grid_1d() noexcept
{
    std::cout << "|---- test_regular_grid_1d()" << std::endl;

    using namespace Math_structures;

    Regular_grid_1d<int> rg_1d_1{1, 2, 3, 4};
    Regular_grid_1d<int> rg_1d_2{4, 3, 3, 4};

    assert(rg_1d_1.size() == 4);
    assert(rg_1d_2.size() == 4);

    assert(rg_1d_1[0] == 1);
    assert(rg_1d_1[1] == 2);
    assert(rg_1d_1[2] == 3);
    assert(rg_1d_1[3] == 4);

    assert(rg_1d_2[0] == 4);
    assert(rg_1d_2[1] == 3);
    assert(rg_1d_2[2] == 3);
    assert(rg_1d_2[3] == 4);

    rg_1d_1 = rg_1d_2;

    assert(rg_1d_1 == rg_1d_2);

    rg_1d_1 = Regular_grid_1d<int>{0, 1, 2, 3};
    rg_1d_2 = Regular_grid_1d<int>{2, 3, 4, 5};

    assert(rg_1d_1 * rg_1d_2 == 26);

    rg_1d_1 += rg_1d_2;

    assert(rg_1d_1 == Regular_grid_1d<int>({2, 4, 6, 8}));

    rg_1d_1 -= rg_1d_2;

    assert(rg_1d_1 == Regular_grid_1d<int>({0, 1, 2, 3}));

    rg_1d_2 *= 10;

    assert(rg_1d_2 == Regular_grid_1d<int>({20, 30, 40, 50}));

    rg_1d_2 /= 10;

    assert(rg_1d_2 == Regular_grid_1d<int>({2, 3, 4, 5}));

    rg_1d_1[0] = 5;
    rg_1d_1[1] = 50;
    rg_1d_1[2] = 500;
    rg_1d_1[3] = 5000;

    assert(rg_1d_1[0] == 5);
    assert(rg_1d_1[1] == 50);
    assert(rg_1d_1[2] == 500);
    assert(rg_1d_1[3] == 5000);

    assert(max_abs(rg_1d_1) == 5000);

    Regular_grid_1d rg_1d_3 = {
        Regular_grid_1d<int>{1, 2, 3, 4, 5},
        Regular_grid_1d<int>{4, 3, 2, 1, 0},
        Regular_grid_1d<int>{9, 8, 7, 6, 5}
    };

    assert(max_abs(rg_1d_3) == 9);
}

}
