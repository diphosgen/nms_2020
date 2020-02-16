#include "test_regular_grid_3d.hxx"

namespace Tests
{

void test_regular_grid_3d() noexcept
{
    std::cout << "|---- test_regular_grid_3d()" << std::endl;

    using namespace Math_structures;

    Regular_grid_3d<int> rg_3d_1{
                                    {
                                        {1, 2, 3, 4, 5},
                                        {6, 7, 8, 9, 0},
                                        {9, 8, 7, 6, 5},
                                        {3, 3, 4, 3, 3}
                                    },
                                    {
                                        {2, 4, 6, 8, 0},
                                        {1, 3, 5, 7, 9},
                                        {1, 2, 5, 6, 9},
                                        {1, 1, 1, 2, 2}
                                    },
                                    {
                                        {9, 8, 7, 6, 5},
                                        {4, 3, 2, 1, 0},
                                        {0, 0, 2, 1, 1},
                                        {9, 7, 5, 3, 1}
                                    }
                                };

    Regular_grid_3d<int> rg_3d_2{
                                    {
                                        {0, 1, 2, 3, 4},
                                        {5, 6, 7, 8, 9},
                                        {0, 1, 0, 2, 0},
                                        {3, 0, 4, 0, 5}
                                    },
                                    {
                                        {1, 2, 1, 0, 1},
                                        {0, 9, 0, 8, 0},
                                        {7, 0, 6, 0, 5},
                                        {0, 4, 0, 3, 0}
                                    },
                                    {
                                        {0, 0, 2, 0, 0},
                                        {2, 1, 0, 1, 2},
                                        {3, 2, 9, 4, 1},
                                        {0, 1, 0, 1, 0}
                                    }
                                };

    assert(rg_3d_1.get_size_x() == 3);
    assert(rg_3d_1.get_size_y() == 4);
    assert(rg_3d_1.get_size_z() == 5);

    assert(rg_3d_2.get_size_x() == 3);
    assert(rg_3d_2.get_size_y() == 4);
    assert(rg_3d_2.get_size_z() == 5);

    assert(rg_3d_1[0][0][0] == 1);
    assert(rg_3d_1[0][0][1] == 2);
    assert(rg_3d_1[0][0][2] == 3);
    assert(rg_3d_1[0][0][3] == 4);
    assert(rg_3d_1[0][0][4] == 5);
    assert(rg_3d_1[0][1][0] == 6);
    assert(rg_3d_1[0][1][1] == 7);
    assert(rg_3d_1[0][1][2] == 8);
    assert(rg_3d_1[0][1][3] == 9);
    assert(rg_3d_1[0][1][4] == 0);
    assert(rg_3d_1[0][2][0] == 9);
    assert(rg_3d_1[0][2][1] == 8);
    assert(rg_3d_1[0][2][2] == 7);
    assert(rg_3d_1[0][2][3] == 6);
    assert(rg_3d_1[0][2][4] == 5);
    assert(rg_3d_1[0][3][0] == 3);
    assert(rg_3d_1[0][3][1] == 3);
    assert(rg_3d_1[0][3][2] == 4);
    assert(rg_3d_1[0][3][3] == 3);
    assert(rg_3d_1[0][3][4] == 3);

    assert(rg_3d_1[1][0][0] == 2);
    assert(rg_3d_1[1][0][1] == 4);
    assert(rg_3d_1[1][0][2] == 6);
    assert(rg_3d_1[1][0][3] == 8);
    assert(rg_3d_1[1][0][4] == 0);
    assert(rg_3d_1[1][1][0] == 1);
    assert(rg_3d_1[1][1][1] == 3);
    assert(rg_3d_1[1][1][2] == 5);
    assert(rg_3d_1[1][1][3] == 7);
    assert(rg_3d_1[1][1][4] == 9);
    assert(rg_3d_1[1][2][0] == 1);
    assert(rg_3d_1[1][2][1] == 2);
    assert(rg_3d_1[1][2][2] == 5);
    assert(rg_3d_1[1][2][3] == 6);
    assert(rg_3d_1[1][2][4] == 9);
    assert(rg_3d_1[1][3][0] == 1);
    assert(rg_3d_1[1][3][1] == 1);
    assert(rg_3d_1[1][3][2] == 1);
    assert(rg_3d_1[1][3][3] == 2);
    assert(rg_3d_1[1][3][4] == 2);

    assert(rg_3d_1[2][0][0] == 9);
    assert(rg_3d_1[2][0][1] == 8);
    assert(rg_3d_1[2][0][2] == 7);
    assert(rg_3d_1[2][0][3] == 6);
    assert(rg_3d_1[2][0][4] == 5);
    assert(rg_3d_1[2][1][0] == 4);
    assert(rg_3d_1[2][1][1] == 3);
    assert(rg_3d_1[2][1][2] == 2);
    assert(rg_3d_1[2][1][3] == 1);
    assert(rg_3d_1[2][1][4] == 0);
    assert(rg_3d_1[2][2][0] == 0);
    assert(rg_3d_1[2][2][1] == 0);
    assert(rg_3d_1[2][2][2] == 2);
    assert(rg_3d_1[2][2][3] == 1);
    assert(rg_3d_1[2][2][4] == 1);
    assert(rg_3d_1[2][3][0] == 9);
    assert(rg_3d_1[2][3][1] == 7);
    assert(rg_3d_1[2][3][2] == 5);
    assert(rg_3d_1[2][3][3] == 3);
    assert(rg_3d_1[2][3][4] == 1);

    assert(rg_3d_2[0][0][0] == 0);
    assert(rg_3d_2[0][0][1] == 1);
    assert(rg_3d_2[0][0][2] == 2);
    assert(rg_3d_2[0][0][3] == 3);
    assert(rg_3d_2[0][0][4] == 4);
    assert(rg_3d_2[0][1][0] == 5);
    assert(rg_3d_2[0][1][1] == 6);
    assert(rg_3d_2[0][1][2] == 7);
    assert(rg_3d_2[0][1][3] == 8);
    assert(rg_3d_2[0][1][4] == 9);
    assert(rg_3d_2[0][2][0] == 0);
    assert(rg_3d_2[0][2][1] == 1);
    assert(rg_3d_2[0][2][2] == 0);
    assert(rg_3d_2[0][2][3] == 2);
    assert(rg_3d_2[0][2][4] == 0);
    assert(rg_3d_2[0][3][0] == 3);
    assert(rg_3d_2[0][3][1] == 0);
    assert(rg_3d_2[0][3][2] == 4);
    assert(rg_3d_2[0][3][3] == 0);
    assert(rg_3d_2[0][3][4] == 5);

    assert(rg_3d_2[1][0][0] == 1);
    assert(rg_3d_2[1][0][1] == 2);
    assert(rg_3d_2[1][0][2] == 1);
    assert(rg_3d_2[1][0][3] == 0);
    assert(rg_3d_2[1][0][4] == 1);
    assert(rg_3d_2[1][1][0] == 0);
    assert(rg_3d_2[1][1][1] == 9);
    assert(rg_3d_2[1][1][2] == 0);
    assert(rg_3d_2[1][1][3] == 8);
    assert(rg_3d_2[1][1][4] == 0);
    assert(rg_3d_2[1][2][0] == 7);
    assert(rg_3d_2[1][2][1] == 0);
    assert(rg_3d_2[1][2][2] == 6);
    assert(rg_3d_2[1][2][3] == 0);
    assert(rg_3d_2[1][2][4] == 5);
    assert(rg_3d_2[1][3][0] == 0);
    assert(rg_3d_2[1][3][1] == 4);
    assert(rg_3d_2[1][3][2] == 0);
    assert(rg_3d_2[1][3][3] == 3);
    assert(rg_3d_2[1][3][4] == 0);

    assert(rg_3d_2[2][0][0] == 0);
    assert(rg_3d_2[2][0][1] == 0);
    assert(rg_3d_2[2][0][2] == 2);
    assert(rg_3d_2[2][0][3] == 0);
    assert(rg_3d_2[2][0][4] == 0);
    assert(rg_3d_2[2][1][0] == 2);
    assert(rg_3d_2[2][1][1] == 1);
    assert(rg_3d_2[2][1][2] == 0);
    assert(rg_3d_2[2][1][3] == 1);
    assert(rg_3d_2[2][1][4] == 2);
    assert(rg_3d_2[2][2][0] == 3);
    assert(rg_3d_2[2][2][1] == 2);
    assert(rg_3d_2[2][2][2] == 9);
    assert(rg_3d_2[2][2][3] == 4);
    assert(rg_3d_2[2][2][4] == 1);
    assert(rg_3d_2[2][3][0] == 0);
    assert(rg_3d_2[2][3][1] == 1);
    assert(rg_3d_2[2][3][2] == 0);
    assert(rg_3d_2[2][3][3] == 1);
    assert(rg_3d_2[2][3][4] == 0);

    assert(!(rg_3d_1 == rg_3d_2));

    assert(max_abs(rg_3d_1) == 9);
    assert(max_abs(rg_3d_2) == 9);

    rg_3d_2[0][0][0] = 1;
    rg_3d_2[0][0][1] = 2;
    rg_3d_2[0][0][2] = 3;
    rg_3d_2[0][0][3] = 4;
    rg_3d_2[0][0][4] = 5;
    rg_3d_2[0][1][0] = 6;
    rg_3d_2[0][1][1] = 7;
    rg_3d_2[0][1][2] = 8;
    rg_3d_2[0][1][3] = 9;
    rg_3d_2[0][1][4] = 0;
    rg_3d_2[0][2][0] = 9;
    rg_3d_2[0][2][1] = 8;
    rg_3d_2[0][2][2] = 7;
    rg_3d_2[0][2][3] = 6;
    rg_3d_2[0][2][4] = 5;
    rg_3d_2[0][3][0] = 3;
    rg_3d_2[0][3][1] = 3;
    rg_3d_2[0][3][2] = 4;
    rg_3d_2[0][3][3] = 3;
    rg_3d_2[0][3][4] = 3;

    rg_3d_2[1][0][0] = 2;
    rg_3d_2[1][0][1] = 4;
    rg_3d_2[1][0][2] = 6;
    rg_3d_2[1][0][3] = 8;
    rg_3d_2[1][0][4] = 0;
    rg_3d_2[1][1][0] = 1;
    rg_3d_2[1][1][1] = 3;
    rg_3d_2[1][1][2] = 5;
    rg_3d_2[1][1][3] = 7;
    rg_3d_2[1][1][4] = 9;
    rg_3d_2[1][2][0] = 1;
    rg_3d_2[1][2][1] = 2;
    rg_3d_2[1][2][2] = 5;
    rg_3d_2[1][2][3] = 6;
    rg_3d_2[1][2][4] = 9;
    rg_3d_2[1][3][0] = 1;
    rg_3d_2[1][3][1] = 1;
    rg_3d_2[1][3][2] = 1;
    rg_3d_2[1][3][3] = 2;
    rg_3d_2[1][3][4] = 2;

    rg_3d_2[2][0][0] = 9;
    rg_3d_2[2][0][1] = 8;
    rg_3d_2[2][0][2] = 7;
    rg_3d_2[2][0][3] = 6;
    rg_3d_2[2][0][4] = 5;
    rg_3d_2[2][1][0] = 4;
    rg_3d_2[2][1][1] = 3;
    rg_3d_2[2][1][2] = 2;
    rg_3d_2[2][1][3] = 1;
    rg_3d_2[2][1][4] = 0;
    rg_3d_2[2][2][0] = 0;
    rg_3d_2[2][2][1] = 0;
    rg_3d_2[2][2][2] = 2;
    rg_3d_2[2][2][3] = 1;
    rg_3d_2[2][2][4] = 1;
    rg_3d_2[2][3][0] = 9;
    rg_3d_2[2][3][1] = 7;
    rg_3d_2[2][3][2] = 5;
    rg_3d_2[2][3][3] = 3;
    rg_3d_2[2][3][4] = 1;

    assert(rg_3d_2[0][0][0] == 1);
    assert(rg_3d_2[0][0][1] == 2);
    assert(rg_3d_2[0][0][2] == 3);
    assert(rg_3d_2[0][0][3] == 4);
    assert(rg_3d_2[0][0][4] == 5);
    assert(rg_3d_2[0][1][0] == 6);
    assert(rg_3d_2[0][1][1] == 7);
    assert(rg_3d_2[0][1][2] == 8);
    assert(rg_3d_2[0][1][3] == 9);
    assert(rg_3d_2[0][1][4] == 0);
    assert(rg_3d_2[0][2][0] == 9);
    assert(rg_3d_2[0][2][1] == 8);
    assert(rg_3d_2[0][2][2] == 7);
    assert(rg_3d_2[0][2][3] == 6);
    assert(rg_3d_2[0][2][4] == 5);
    assert(rg_3d_2[0][3][0] == 3);
    assert(rg_3d_2[0][3][1] == 3);
    assert(rg_3d_2[0][3][2] == 4);
    assert(rg_3d_2[0][3][3] == 3);
    assert(rg_3d_2[0][3][4] == 3);

    assert(rg_3d_2[1][0][0] == 2);
    assert(rg_3d_2[1][0][1] == 4);
    assert(rg_3d_2[1][0][2] == 6);
    assert(rg_3d_2[1][0][3] == 8);
    assert(rg_3d_2[1][0][4] == 0);
    assert(rg_3d_2[1][1][0] == 1);
    assert(rg_3d_2[1][1][1] == 3);
    assert(rg_3d_2[1][1][2] == 5);
    assert(rg_3d_2[1][1][3] == 7);
    assert(rg_3d_2[1][1][4] == 9);
    assert(rg_3d_2[1][2][0] == 1);
    assert(rg_3d_2[1][2][1] == 2);
    assert(rg_3d_2[1][2][2] == 5);
    assert(rg_3d_2[1][2][3] == 6);
    assert(rg_3d_2[1][2][4] == 9);
    assert(rg_3d_2[1][3][0] == 1);
    assert(rg_3d_2[1][3][1] == 1);
    assert(rg_3d_2[1][3][2] == 1);
    assert(rg_3d_2[1][3][3] == 2);
    assert(rg_3d_2[1][3][4] == 2);

    assert(rg_3d_2[2][0][0] == 9);
    assert(rg_3d_2[2][0][1] == 8);
    assert(rg_3d_2[2][0][2] == 7);
    assert(rg_3d_2[2][0][3] == 6);
    assert(rg_3d_2[2][0][4] == 5);
    assert(rg_3d_2[2][1][0] == 4);
    assert(rg_3d_2[2][1][1] == 3);
    assert(rg_3d_2[2][1][2] == 2);
    assert(rg_3d_2[2][1][3] == 1);
    assert(rg_3d_2[2][1][4] == 0);
    assert(rg_3d_2[2][2][0] == 0);
    assert(rg_3d_2[2][2][1] == 0);
    assert(rg_3d_2[2][2][2] == 2);
    assert(rg_3d_2[2][2][3] == 1);
    assert(rg_3d_2[2][2][4] == 1);
    assert(rg_3d_2[2][3][0] == 9);
    assert(rg_3d_2[2][3][1] == 7);
    assert(rg_3d_2[2][3][2] == 5);
    assert(rg_3d_2[2][3][3] == 3);
    assert(rg_3d_2[2][3][4] == 1);

    assert(rg_3d_1 == rg_3d_2);
}

}
