#ifndef HEAT_CONDUCTION_EQUATION_1D_HXX_INCLUDED
#define HEAT_CONDUCTION_EQUATION_1D_HXX_INCLUDED

#include "boundary_conditions.hxx"
#include "regular_grid_1d.hxx"
#include "one_step_ode_solver.hxx"

namespace Math_structures
{

template<
    Boundary_conditions lbc = Boundary_conditions::Dirichlet,
    Boundary_conditions rbc = Boundary_conditions::Dirichlet,
    typename T = double
    >
class Heat_conduction_equation_1d
{
private:

    using Function_T  = std::function<std::decay_t<T>(std::decay_t<T>)>;
    using Function_TT = std::function<std::decay_t<T>(std::decay_t<T>, std::decay_t<T>)>;

    template<Boundary_conditions condition>
    using Boundary_conditions_function = std::conditional_t<
                                                condition == Boundary_conditions::Robin,
                                                std::pair<Function_T, Function_T>,
                                                Function_T
                                                >;

public:

    using Source_function = Function_TT;
    using Initial_conditions_function = Function_T;

    using Function_lbc = Boundary_conditions_function<lbc>;
    using Function_rbc = Boundary_conditions_function<rbc>;

    using ODE_solver = One_step_ODE_solver<Regular_grid_1d<T>, T>;

    constexpr Heat_conduction_equation_1d() noexcept = delete;

    constexpr Heat_conduction_equation_1d(const Heat_conduction_equation_1d<lbc, rbc, T>&) noexcept = default;
    constexpr Heat_conduction_equation_1d(Heat_conduction_equation_1d<lbc, rbc, T>&&) noexcept = default;

    constexpr Heat_conduction_equation_1d(const Initial_conditions_function& init_conditions,
                                          const Function_lbc& function_lbc,
                                          const Function_rbc& function_rbc,
                                          const Source_function& source,
                                          const int grid_steps_amount,
                                          const T& lbc_coordinate = 0.0,
                                          const T& rbc_coordinate = 1.0,
                                          const T& thermal_diffusivity = 1.0) noexcept
        :   init_conditions{init_conditions},
            function_lbc{function_lbc},
            function_rbc{function_rbc},
            source{source},
            grid_steps_amount{grid_steps_amount},
            grid_nodes_amount{grid_steps_amount + 1},
            lbc_coordinate{lbc_coordinate},
            rbc_coordinate{rbc_coordinate},
            grid_spatial_step{(rbc_coordinate - lbc_coordinate)/grid_steps_amount},
            thermal_diffusivity{thermal_diffusivity}
    {
        assert(this->grid_steps_amount > 1);
        assert(this->rbc_coordinate > this->lbc_coordinate);
        assert(this->thermal_diffusivity > 0.0);
    }

    ~Heat_conduction_equation_1d() noexcept = default;

    constexpr Heat_conduction_equation_1d& operator=
        (const Heat_conduction_equation_1d<lbc, rbc, T>&) noexcept = default;

    constexpr Heat_conduction_equation_1d& operator=
        (Heat_conduction_equation_1d<lbc, rbc, T>&&) noexcept = default;

    Regular_grid_1d<T> get_initial_conditions_grid() const noexcept
    {
        Regular_grid_1d<T> grid_solution(this->grid_nodes_amount);

        for (int i = 0; i < this->grid_nodes_amount; ++i) {
            const T x = this->lbc_coordinate + i * this->grid_spatial_step;
            grid_solution[i] = this->init_conditions(x);
        }

        return grid_solution;
    }

    Regular_grid_1d<T> get_solution(const Regular_grid_1d<T>& prev_grid_solution,
                                    const std::unique_ptr<ODE_solver>& ode_solver_ptr,
                                    const T& t, const T& dt) const
    {
        Regular_grid_1d<T> grid_solution(this->grid_nodes_amount);

        auto f_rhs = [this](const Regular_grid_1d<T>& solution, const T& t)
            noexcept -> Regular_grid_1d<T>
        {
            return this->get_rhs_ode_system(solution, t);
        };

        grid_solution = ode_solver_ptr->get_solution(f_rhs, prev_grid_solution, t, dt);

        this->set_formation_dirichlet_boundary_conditions(grid_solution, t, dt);

        return grid_solution;
    }

private:

    const Initial_conditions_function init_conditions{};

    const Function_lbc function_lbc{};
    const Function_rbc function_rbc{};

    const Source_function source{};

    const int grid_steps_amount = 0;
    const int grid_nodes_amount{};

    const T lbc_coordinate = T(0.0);
    const T rbc_coordinate = T(1.0);

    const T grid_spatial_step{};

    const T thermal_diffusivity = T(1.0);

    Regular_grid_1d<T> get_rhs_ode_system
        (const Regular_grid_1d<T>& grid_solution, const T& t) const noexcept
    {
        Regular_grid_1d<T> rhs_ode_system(this->grid_nodes_amount);

        this->insert_boundary_conditions(rhs_ode_system, grid_solution, t);

        for (int i = 1; i < this->grid_nodes_amount-1; ++i) {
            const T h_sqr = this->grid_spatial_step * this->grid_spatial_step;
            const T sec_finite_diff = grid_solution[i+1] - 2 * grid_solution[i] + grid_solution[i-1];
            const T laplace_deriv = sec_finite_diff / h_sqr;
            const T x = this->lbc_coordinate + i * this->grid_spatial_step;
            rhs_ode_system[i] = this->thermal_diffusivity * laplace_deriv + this->source(x, t);
        }

        return rhs_ode_system;
    }

    void set_formation_dirichlet_boundary_conditions
        (Regular_grid_1d<T>& grid_solution, const T& t, const T& dt) const noexcept
    {
        if constexpr (lbc == Boundary_conditions::Dirichlet) {
            const int lbc_node = 0;
            grid_solution[lbc_node] = this->function_lbc(t + dt);
        }

        if constexpr (rbc == Boundary_conditions::Dirichlet) {
            const int rbc_node = this->grid_steps_amount;
            grid_solution[rbc_node] = this->function_rbc(t + dt);
        }
    }

    void insert_boundary_conditions(Regular_grid_1d<T>& rhs_ode_system,
                                    const Regular_grid_1d<T>& grid_solution,
                                    const T& t) const noexcept
    {
        insert_lbc(rhs_ode_system, grid_solution, t);
        insert_rbc(rhs_ode_system, grid_solution, t);
    }

    void insert_lbc(Regular_grid_1d<T>& rhs_ode_system,
                    const Regular_grid_1d<T>& grid_solution,
                    const T& t) const noexcept
    {
        if constexpr (lbc == Boundary_conditions::Dirichlet) {
            const int lbc_node = 0;
            rhs_ode_system[lbc_node] = T(0.0);
        } else if constexpr (lbc == Boundary_conditions::Neumann) {
            const int lbc_node = 0;
            const T h_sqr = this->grid_spatial_step * this->grid_spatial_step;
            const T multiplier = grid_solution[1] - grid_solution[0] -
                                 this->grid_spatial_step * this->function_lbc(t);
            rhs_ode_system[lbc_node] = 2 * this->thermal_diffusivity * multiplier / h_sqr +
                                       this->source(this->lbc_coordinate, t);
        } else if constexpr (lbc == Boundary_conditions::Robin) {
            const int lbc_node = 0;
            const T h_sqr = this->grid_spatial_step * this->grid_spatial_step;
            constexpr int c_number = 0;
            constexpr int m_number = 1;
            const T multiplier = grid_solution[1] - grid_solution[0] -
                                 this->grid_spatial_step * (std::get<m_number>(this->function_lbc)(t) -
                                                            std::get<c_number>(this->function_lbc)(t) *
                                                            grid_solution[0]);
            rhs_ode_system[lbc_node] = 2 * this->thermal_diffusivity * multiplier / h_sqr +
                                       this->source(this->lbc_coordinate, t);
        } else { }
    }

    void insert_rbc(Regular_grid_1d<T>& rhs_ode_system,
                    const Regular_grid_1d<T>& grid_solution,
                    const T& t) const noexcept
    {
        if constexpr (rbc == Boundary_conditions::Dirichlet) {
            const int rbc_node = this->grid_steps_amount;
            rhs_ode_system[rbc_node] = T(0.0);
        } else if constexpr (rbc == Boundary_conditions::Neumann) {
            const int rbc_node = this->grid_steps_amount;
            const T h_sqr = this->grid_spatial_step * this->grid_spatial_step;
            const T multiplier = grid_solution[this->grid_steps_amount-1] -
                                 grid_solution[this->grid_steps_amount] +
                                 this->grid_spatial_step * this->function_rbc(t);
            rhs_ode_system[rbc_node] = 2 * this->thermal_diffusivity * multiplier / h_sqr +
                                       this->source(this->rbc_coordinate, t);
        } else if constexpr (rbc == Boundary_conditions::Robin) {
            const int rbc_node = this->grid_steps_amount;
            const T h_sqr = this->grid_spatial_step * this->grid_spatial_step;
            constexpr int c_number = 0;
            constexpr int m_number = 1;
            const T multiplier = grid_solution[this->grid_steps_amount-1] -
                                 grid_solution[this->grid_steps_amount] +
                                 this->grid_spatial_step * (std::get<m_number>(this->function_rbc)(t) -
                                                            std::get<c_number>(this->function_rbc)(t) *
                                                            grid_solution[this->grid_steps_amount]);
            rhs_ode_system[rbc_node] = 2 * this->thermal_diffusivity * multiplier / h_sqr +
                                       this->source(this->rbc_coordinate, t);
        } else { }
    }
};

}

#endif // HEAT_CONDUCTION_EQUATION_1D_HXX_INCLUDED
