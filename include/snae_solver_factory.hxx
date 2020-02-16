#ifndef SNAE_SOLVER_FACTORY_HXX_INCLUDED
#define SNAE_SOLVER_FACTORY_HXX_INCLUDED

#include "snae_solver.hxx"
#include "snae_solver_fpi.hxx"

namespace Math_structures
{

template<typename T_args>
class SNAE_solver_factory
{
public:

    constexpr SNAE_solver_factory() noexcept = default;

	constexpr SNAE_solver_factory(const SNAE_solver_factory<T_args>&) noexcept = default;
	constexpr SNAE_solver_factory(SNAE_solver_factory<T_args>&&) noexcept = default;

    virtual ~SNAE_solver_factory() noexcept = default;

	constexpr SNAE_solver_factory<T_args>&
		operator=(const SNAE_solver_factory<T_args>&) noexcept = default;

	constexpr SNAE_solver_factory<T_args>&
		operator=(SNAE_solver_factory<T_args>&&) noexcept = default;

    virtual std::unique_ptr<SNAE_solver<T_args>> create() const noexcept = 0;

};

template<typename T_args>
class SNAE_solver_default_factory
    : public SNAE_solver_factory<T_args>
{
public:

	constexpr SNAE_solver_default_factory() noexcept = default;

	constexpr SNAE_solver_default_factory(const SNAE_solver_default_factory<T_args>&) noexcept = default;
	constexpr SNAE_solver_default_factory(SNAE_solver_default_factory<T_args>&&) noexcept = default;

    virtual ~SNAE_solver_default_factory() noexcept = default;

	constexpr SNAE_solver_default_factory<T_args>&
		operator=(const SNAE_solver_default_factory<T_args>&) noexcept = default;

	constexpr SNAE_solver_default_factory<T_args>&
		operator=(SNAE_solver_default_factory<T_args>&&) noexcept = default;

    virtual std::unique_ptr<SNAE_solver<T_args>> create() const noexcept
    {
        return std::make_unique<SNAE_solver_FPI<T_args>>();
    }

};

}

#endif // SNAE_SOLVER_FACTORY_HXX_INCLUDED
