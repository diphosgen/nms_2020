#ifndef SLAE_SOLVER_FACTORY_HXX_INCLUDED
#define SLAE_SOLVER_FACTORY_HXX_INCLUDED

#include "slae_solver.hxx"
#include "slae_solver_gauss.hxx"

namespace Math_structures
{

template<
    typename Matrix_type = Matrix_unfixed<double>,
    typename Vector_type = Vector_unfixed<double>
    >
class SLAE_solver_factory
{
public:

    constexpr SLAE_solver_factory() noexcept = default;

	constexpr SLAE_solver_factory(const SLAE_solver_factory<Matrix_type, Vector_type>&) noexcept = default;
	constexpr SLAE_solver_factory(SLAE_solver_factory<Matrix_type, Vector_type>&&) noexcept = default;

    virtual ~SLAE_solver_factory() noexcept = default;

	constexpr SLAE_solver_factory<Matrix_type, Vector_type>&
		operator=(const SLAE_solver_factory<Matrix_type, Vector_type>&) noexcept = default;

	constexpr SLAE_solver_factory<Matrix_type, Vector_type>&
		operator=(SLAE_solver_factory<Matrix_type, Vector_type>&&) noexcept = default;

    virtual std::unique_ptr<SLAE_solver<Matrix_type, Vector_type>> create() const noexcept = 0;

};

template<
    typename Matrix_type = Matrix_unfixed<double>,
    typename Vector_type = Vector_unfixed<double>
    >
class SLAE_solver_default_factory
    : public SLAE_solver_factory<Matrix_type, Vector_type>
{
public:

	constexpr SLAE_solver_default_factory() noexcept = default;

	constexpr SLAE_solver_default_factory
        (const SLAE_solver_default_factory<Matrix_type, Vector_type>&) noexcept = default;

	constexpr SLAE_solver_default_factory
        (SLAE_solver_default_factory<Matrix_type, Vector_type>&&) noexcept = default;

    virtual ~SLAE_solver_default_factory() noexcept = default;

	constexpr SLAE_solver_default_factory<Matrix_type, Vector_type>&
		operator=(const SLAE_solver_default_factory<Matrix_type, Vector_type>&) noexcept = default;

	constexpr SLAE_solver_default_factory<Matrix_type, Vector_type>&
		operator=(SLAE_solver_default_factory<Matrix_type, Vector_type>&&) noexcept = default;

    virtual std::unique_ptr<SLAE_solver<Matrix_type, Vector_type>> create() const noexcept
    {
        return std::make_unique<SLAE_solver_Gauss<Matrix_type, Vector_type>>();
    }

};

}

#endif // SLAE_SOLVER_FACTORY_HXX_INCLUDED
