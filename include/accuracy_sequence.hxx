#ifndef ACCURACY_SEQUENCE_HXX_INCLUDED
#define ACCURACY_SEQUENCE_HXX_INCLUDED

#include "basic_libs.hxx"
#include "basic_math_functions.hxx"

namespace Math_structures
{

class Accuracy_sequence
{
public:

    constexpr Accuracy_sequence() noexcept = default;

    constexpr Accuracy_sequence(const Accuracy_sequence&) noexcept = default;
    constexpr Accuracy_sequence(Accuracy_sequence&&) noexcept = default;

    virtual ~Accuracy_sequence() noexcept = default;

    Accuracy_sequence& operator=(const Accuracy_sequence&) noexcept = default;
    Accuracy_sequence& operator=(Accuracy_sequence&&) noexcept = default;

    virtual int get_sequence_element(const int number) const noexcept = 0;
};

class Richardson_sequence
    :   public Accuracy_sequence
{
public:

    constexpr Richardson_sequence() noexcept = default;

    constexpr Richardson_sequence(const Richardson_sequence&) noexcept = default;
    constexpr Richardson_sequence(Richardson_sequence&&) noexcept = default;

    virtual ~Richardson_sequence() noexcept = default;

    Richardson_sequence& operator=(const Richardson_sequence&) noexcept = default;
    Richardson_sequence& operator=(Richardson_sequence&&) noexcept = default;

    virtual int get_sequence_element(const int number) const noexcept override
    {
        return number + 1;
    }
};

class Romberg_sequence
    :   public Accuracy_sequence
{
public:

    constexpr Romberg_sequence() noexcept = default;

    constexpr Romberg_sequence(const Romberg_sequence&) noexcept = default;
    constexpr Romberg_sequence(Romberg_sequence&&) noexcept = default;

    virtual ~Romberg_sequence() noexcept = default;

    Romberg_sequence& operator=(const Romberg_sequence&) noexcept = default;
    Romberg_sequence& operator=(Romberg_sequence&&) noexcept = default;

    virtual int get_sequence_element(const int number) const noexcept override
    {
        const int element = (number == 0) ? 1 : std::pow(2, number);
        return element;
    }
};

}

#endif // ACCURACY_SEQUENCE_HXX_INCLUDED
