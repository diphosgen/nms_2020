#ifndef RANDOM_HPP_INCLUDED
#define RANDOM_HPP_INCLUDED

#include "basic_libs.hxx"

namespace Math_structures
{

template<
    typename Genrand = std::mt19937,
    typename = std::enable_if_t<
                    !std::is_same_v<Genrand, std::random_device>
                    >
    >
class Random
{
public:

    constexpr Random() noexcept
        :   gen_rand(static_cast<typename Genrand::result_type>(seed_default))
    { }

    constexpr Random(const Random<Genrand>& rng) noexcept
        :   gen_rand{rng.gen_rand}
    { }

    constexpr Random(Random<Genrand>&& rng) noexcept
        :   gen_rand{std::move(rng.gen_rand)}
    { }

    constexpr Random(const int& seed) noexcept
        :   gen_rand(static_cast<typename Genrand::result_type>(seed))
    { }

    ~Random() noexcept = default;

    constexpr Random<Genrand>& operator=(const Random<Genrand>& rng) noexcept
    {
        if (&rng != this) {
            gen_rand = rng.gen_rand;
        }

        return *this;
    }

    constexpr Random<Genrand>& operator=(Random<Genrand>&& rng) noexcept
    {
        if (&rng != this) {
            gen_rand = std::move(rng.gen_rand);
        }

        return *this;
    }

    constexpr void cooldown(const int& seed) noexcept
    {
        this->gen_rand.seed(static_cast<typename Genrand::result_type>(seed));
    }

    constexpr double operator()() noexcept
    {
        constexpr double temp_denominator = 1.0 / static_cast<double>(Genrand::max() - Genrand::min());
        const double res = temp_denominator * static_cast<double>(gen_rand() - Genrand::min());
        return res;
    }

    constexpr double operator()(const double begin_number, const double end_number) noexcept
    {
        assert(end_number >= begin_number);
        return begin_number + (end_number - begin_number) * (*this)();
    }

    constexpr int operator()(const int begin_number, const int end_number) noexcept
    {
        assert(end_number >= begin_number);
        return begin_number + static_cast<int>((end_number - begin_number) * (*this)());
    }

    constexpr double angle_pi() noexcept
    {
        constexpr double angle_pi = M_PI;
        constexpr double temp_denominator = 1.0 / (Genrand::max() - Genrand::min());
        return M_PI * temp_denominator * static_cast<double>(gen_rand() - Genrand::min());
    }

    constexpr double angle_2pi() noexcept
    {
        constexpr double temp_denominator = 1.0 / (Genrand::max() - Genrand::min());
        return 2 * M_PI * temp_denominator * static_cast<double>(gen_rand() - Genrand::min());
    }

private:

    static constexpr const int seed_default = 0;
    Genrand gen_rand{seed_default};

};

inline int get_seed(const int init = 0) noexcept
{
    using Nanoseconds = std::chrono::nanoseconds;
    using ClockType = std::chrono::system_clock;
    using TimeType = std::chrono::time_point<ClockType>;

    const TimeType timeStart = static_cast<TimeType>(Nanoseconds(init));
    TimeType timeEnd = ClockType::now();

    int timeInterval =
        static_cast<int>(std::chrono::duration_cast<Nanoseconds>(timeEnd - timeStart).count());

    int seed = std::abs(timeInterval);

    return seed;
}

}

#endif // RANDOM_HPP_INCLUDED
