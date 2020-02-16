#include "test_vector_unfixed.hxx"

namespace Tests
{

void test_vector_unfixed() noexcept
{
    std::cout << "|---- test_vector_unfixed()" << std::endl;

    using namespace Math_structures;

    Vector_unfixed<int> v1{1, 2, 3, 4};
    Vector_unfixed<int> v2{4, 3, 3, 4};

    assert(v1.size() == 4);
    assert(v2.size() == 4);

    assert(v1[0] == 1);
    assert(v1[1] == 2);
    assert(v1[2] == 3);
    assert(v1[3] == 4);

    assert(v2[0] == 4);
    assert(v2[1] == 3);
    assert(v2[2] == 3);
    assert(v2[3] == 4);

    v1 = v2;

    assert(v1 == v2);

    v1 = Vector_unfixed<int>{0, 1, 2, 3};
    v2 = Vector_unfixed<int>{2, 3, 4, 5};

    assert(v1 * v2 == 26);

    v1 += v2;

    assert(v1 == Vector_unfixed<int>({2, 4, 6, 8}));

    v1 -= v2;

    assert(v1 == Vector_unfixed<int>({0, 1, 2, 3}));

    v2 *= 10;

    assert(v2 == Vector_unfixed<int>({20, 30, 40, 50}));

    v2 /= 10;

    assert(v2 == Vector_unfixed<int>({2, 3, 4, 5}));

    v1[0] = 5;
    v1[1] = 50;
    v1[2] = 500;
    v1[3] = 5000;

    assert(v1[0] == 5);
    assert(v1[1] == 50);
    assert(v1[2] == 500);
    assert(v1[3] == 5000);

    assert(max_abs(v1) == 5000);

    Vector_unfixed v3 = {
        Vector_unfixed<int>{1, 2, 3, 4, 5},
        Vector_unfixed<int>{4, 3, 2, 1, 0},
        Vector_unfixed<int>{9, 8, 7, 6, 5}
    };

    assert(max_abs(v3) == 9);
}

}
