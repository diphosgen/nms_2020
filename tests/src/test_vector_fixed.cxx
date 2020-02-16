#include "test_vector_fixed.hxx"

namespace Tests
{

void test_vector_fixed() noexcept
{
    std::cout << "|---- test_vector_fixed()" << std::endl;

    using namespace Math_structures;

    constexpr Vector_fixed<int, 4> v1_{1, 2, 3, 4};
    constexpr Vector_fixed<int, 4> v2_{4, 3, 3, 4};

    static_assert(v1_.size() == 4);
    static_assert(v2_.size() == 4);

    static_assert(v1_[0] == 1);
    static_assert(v1_[1] == 2);
    static_assert(v1_[2] == 3);
    static_assert(v1_[3] == 4);

    static_assert(v2_[0] == 4);
    static_assert(v2_[1] == 3);
    static_assert(v2_[2] == 3);
    static_assert(v2_[3] == 4);

    assert(v1_.size() == 4);
    assert(v2_.size() == 4);

    assert(v1_[0] == 1);
    assert(v1_[1] == 2);
    assert(v1_[2] == 3);
    assert(v1_[3] == 4);

    assert(v2_[0] == 4);
    assert(v2_[1] == 3);
    assert(v2_[2] == 3);
    assert(v2_[3] == 4);

    Vector_fixed<int, 4> v1{v1_};
    Vector_fixed<int, 4> v2{v2_};

    v1 = v2;

    assert(v1 == v2);

    v1 = Vector_fixed<int, 4>{0, 1, 2, 3};
    v2 = Vector_fixed<int, 4>{2, 3, 4, 5};

    assert(v1 * v2 == 26);

    v1 += v2;

    assert((v1 == Vector_fixed<int, 4>({2, 4, 6, 8})));

    v1 -= v2;

    assert((v1 == Vector_fixed<int, 4>({0, 1, 2, 3})));

    v2 *= 10;

    assert((v2 == Vector_fixed<int, 4>({20, 30, 40, 50})));

    v2 /= 10;

    assert((v2 == Vector_fixed<int, 4>({2, 3, 4, 5})));

    v1[0] = 5;
    v1[1] = 50;
    v1[2] = 500;
    v1[3] = 5000;

    assert(v1[0] == 5);
    assert(v1[1] == 50);
    assert(v1[2] == 500);
    assert(v1[3] == 5000);

    assert(max_abs(v1) == 5000);

    constexpr Vector_fixed<Vector_fixed<int, 5>, 3> v3 = {
    Vector_fixed<int, 5>{{1, 2, 3, 4, 5}},
    Vector_fixed<int, 5>{{4, 3, 2, 1, 0}},
    Vector_fixed<int, 5>{{9, 8, 7, 6, 5}}
    };

    static_assert(max_abs(v3) == 9);
    assert(max_abs(v3) == 9);
}

}
