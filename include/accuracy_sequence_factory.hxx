#ifndef ACCURACY_SEQUENCE_FACTORY_HXX_INCLUDED
#define ACCURACY_SEQUENCE_FACTORY_HXX_INCLUDED

#include "accuracy_sequence.hxx"

namespace Math_structures
{

class Accuracy_sequence_factory
{
public:

    constexpr Accuracy_sequence_factory() noexcept = default;

	constexpr Accuracy_sequence_factory(const Accuracy_sequence_factory&) noexcept = default;
	constexpr Accuracy_sequence_factory(Accuracy_sequence_factory&&) noexcept = default;

    virtual ~Accuracy_sequence_factory() noexcept = default;

	Accuracy_sequence_factory&
		operator=(const Accuracy_sequence_factory&) noexcept = default;

	Accuracy_sequence_factory&
		operator=(Accuracy_sequence_factory&&) noexcept = default;

    virtual std::unique_ptr<Accuracy_sequence> create() const noexcept = 0;

};

class Accuracy_sequence_default_factory
    : public Accuracy_sequence_factory
{
public:

	constexpr Accuracy_sequence_default_factory() noexcept = default;

	constexpr Accuracy_sequence_default_factory(const Accuracy_sequence_default_factory&) noexcept = default;
	constexpr Accuracy_sequence_default_factory(Accuracy_sequence_default_factory&&) noexcept = default;

    virtual ~Accuracy_sequence_default_factory() noexcept = default;

	Accuracy_sequence_default_factory&
		operator=(const Accuracy_sequence_default_factory&) noexcept = default;

	Accuracy_sequence_default_factory&
		operator=(Accuracy_sequence_default_factory&&) noexcept = default;

    virtual std::unique_ptr<Accuracy_sequence> create() const noexcept
    {
        return std::make_unique<Romberg_sequence>();
    }

};

}

#endif // ACCURACY_SEQUENCE_FACTORY_HXX_INCLUDED
