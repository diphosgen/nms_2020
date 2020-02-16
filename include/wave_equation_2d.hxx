#ifndef WAVE_EQUATION_2D_HXX_INCLUDED
#define WAVE_EQUATION_2D_HXX_INCLUDED

#include "boundary_conditions.hxx"
#include "vector_fixed.hxx"
#include "regular_grid_2d.hxx"
#include "one_step_ode_solver.hxx"

namespace Math_structures
{

template<
    Boundary_conditions lbc_x = Boundary_conditions::Dirichlet,
    Boundary_conditions rbc_x = Boundary_conditions::Dirichlet,
    Boundary_conditions lbc_y = Boundary_conditions::Dirichlet,
    Boundary_conditions rbc_y = Boundary_conditions::Dirichlet,
    typename T = double
    >
class Wave_equation_2d
{
private:

    using Function_T   = std::function<std::decay_t<T>(std::decay_t<T>)>;
    using Function_TT  = std::function<std::decay_t<T>(std::decay_t<T>, std::decay_t<T>)>;
    using Function_TTT = std::function<std::decay_t<T>(std::decay_t<T>, std::decay_t<T>, std::decay_t<T>)>;

    template<Boundary_conditions condition>
    using Boundary_conditions_function = std::conditional_t<
                                                condition == Boundary_conditions::Robin,
                                                std::pair<Function_TT, Function_TT>,
                                                Function_TT
                                                >;

public:

    static constexpr int grid_variables_amount = 2;

    static constexpr int grid_variable_displacement_number = 0;
    static constexpr int grid_variable_velocity_number = 1;

    using Grid_variable = Vector_fixed<std::decay_t<T>, grid_variables_amount>;

    using Source_function = Function_TTT;
    using Initial_conditions_function = std::function<Grid_variable(std::decay_t<T>, std::decay_t<T>)>;

    using Function_lbc_x = Boundary_conditions_function<lbc_x>;
    using Function_rbc_x = Boundary_conditions_function<rbc_x>;
    using Function_lbc_y = Boundary_conditions_function<lbc_y>;
    using Function_rbc_y = Boundary_conditions_function<rbc_y>;

    using ODE_solver = One_step_ODE_solver<Regular_grid_2d<Grid_variable>, T>;

    constexpr Wave_equation_2d() noexcept = delete;

    constexpr Wave_equation_2d
        (const Wave_equation_2d<lbc_x, rbc_x, lbc_y, rbc_y, T>&) noexcept = default;

    constexpr Wave_equation_2d
        (Wave_equation_2d<lbc_x, rbc_x, lbc_y, rbc_y, T>&&) noexcept = default;

    constexpr Wave_equation_2d(const Initial_conditions_function& init_conditions,
                               const Function_lbc_x& function_lbc_x,
                               const Function_rbc_x& function_rbc_x,
                               const Function_lbc_y& function_lbc_y,
                               const Function_rbc_y& function_rbc_y,
                               const Source_function& source,
                               const int grid_steps_amount_x,
                               const int grid_steps_amount_y,
                               const T& lbc_x_coordinate = 0.0,
                               const T& rbc_x_coordinate = 1.0,
                               const T& lbc_y_coordinate = 0.0,
                               const T& rbc_y_coordinate = 1.0,
                               const T& wave_speed = 1.0) noexcept
        :   init_conditions{init_conditions},
            function_lbc_x{function_lbc_x},
            function_rbc_x{function_rbc_x},
            function_lbc_y{function_lbc_y},
            function_rbc_y{function_rbc_y},
            source{source},
            grid_steps_amount_x{grid_steps_amount_x},
            grid_steps_amount_y{grid_steps_amount_y},
            grid_nodes_amount_x{grid_steps_amount_x + 1},
            grid_nodes_amount_y{grid_steps_amount_y + 1},
            lbc_x_coordinate{lbc_x_coordinate},
            rbc_x_coordinate{rbc_x_coordinate},
            lbc_y_coordinate{lbc_y_coordinate},
            rbc_y_coordinate{rbc_y_coordinate},
            grid_spatial_step_x{(rbc_x_coordinate - lbc_x_coordinate)/grid_steps_amount_x},
            grid_spatial_step_y{(rbc_y_coordinate - lbc_y_coordinate)/grid_steps_amount_y},
            wave_speed{wave_speed}
    {
        assert(this->grid_steps_amount_x > 1);
        assert(this->grid_steps_amount_y > 1);
        assert(this->rbc_x_coordinate > this->lbc_x_coordinate);
        assert(this->rbc_y_coordinate > this->lbc_y_coordinate);
        assert(this->wave_speed > 0.0);
    }

    ~Wave_equation_2d() noexcept = default;

    constexpr Wave_equation_2d& operator=
        (const Wave_equation_2d<lbc_x, rbc_x, lbc_y, rbc_y, T>&) noexcept = default;

    constexpr Wave_equation_2d& operator=
        (Wave_equation_2d<lbc_x, rbc_x, lbc_y, rbc_y, T>&&) noexcept = default;

    Regular_grid_2d<Grid_variable> get_initial_conditions_grid() const noexcept
    {
        Regular_grid_2d<Grid_variable> grid_solution(this->grid_nodes_amount_x,
                                                     this->grid_nodes_amount_y);

        for (int i = 0; i < this->grid_nodes_amount_x; ++i) {
            for (int j = 0; j < this->grid_nodes_amount_y; ++j) {
                const T x = this->lbc_x_coordinate + i * this->grid_spatial_step_x;
                const T y = this->lbc_y_coordinate + j * this->grid_spatial_step_y;
                grid_solution[i][j] = this->init_conditions(x, y);
            }
        }

        return grid_solution;
    }

    Regular_grid_2d<Grid_variable> get_solution(const Regular_grid_2d<Grid_variable>& prev_grid_solution,
                                                const std::unique_ptr<ODE_solver>& ode_solver_ptr,
                                                const T& t, const T& dt) const
    {
        Regular_grid_2d<Grid_variable> grid_solution(this->grid_nodes_amount_x,
                                                     this->grid_nodes_amount_y);

        auto f_rhs = [this](const Regular_grid_2d<Grid_variable>& solution, const T& t)
            noexcept -> Regular_grid_2d<Grid_variable>
        {
            return this->get_rhs_ode_system(solution, t);
        };

        grid_solution = ode_solver_ptr->get_solution(f_rhs, prev_grid_solution, t, dt);

        this->set_formation_dirichlet_boundary_conditions(grid_solution, t, dt);

        return grid_solution;
    }

private:

    const Initial_conditions_function init_conditions{};

    const Function_lbc_x function_lbc_x{};
    const Function_rbc_x function_rbc_x{};
    const Function_lbc_y function_lbc_y{};
    const Function_rbc_y function_rbc_y{};

    const Source_function source{};

    const int grid_steps_amount_x = 0;
    const int grid_steps_amount_y = 0;

    const int grid_nodes_amount_x{};
    const int grid_nodes_amount_y{};

    const T lbc_x_coordinate = T(0.0);
    const T rbc_x_coordinate = T(1.0);
    const T lbc_y_coordinate = T(0.0);
    const T rbc_y_coordinate = T(1.0);

    const T grid_spatial_step_x{};
    const T grid_spatial_step_y{};

    const T wave_speed = T(1.0);

    Regular_grid_2d<Grid_variable> get_rhs_ode_system
        (const Regular_grid_2d<Grid_variable>& grid_solution, const T& t) const noexcept
    {
        Regular_grid_2d<Grid_variable> rhs_ode_system(this->grid_nodes_amount_x,
                                                      this->grid_nodes_amount_y);

        this->insert_boundary_conditions(rhs_ode_system, grid_solution, t);

        for (int i = 1; i < this->grid_nodes_amount_x-1; ++i) {
            for (int j = 1; j < this->grid_nodes_amount_y-1; ++j) {

                const T sec_finite_diff_x = grid_solution[i+1][j][grid_variable_displacement_number] -
                                            2 * grid_solution[i][j][grid_variable_displacement_number] +
                                            grid_solution[i-1][j][grid_variable_displacement_number];

                const T sec_finite_diff_y = grid_solution[i][j+1][grid_variable_displacement_number] -
                                            2 * grid_solution[i][j][grid_variable_displacement_number] +
                                            grid_solution[i][j-1][grid_variable_displacement_number];

                const T h_x_sqr = this->grid_spatial_step_x * this->grid_spatial_step_x;
                const T h_y_sqr = this->grid_spatial_step_y * this->grid_spatial_step_y;

                const T laplace_deriv_x = sec_finite_diff_x / h_x_sqr;
                const T laplace_deriv_y = sec_finite_diff_y / h_y_sqr;

                const T laplace_deriv = laplace_deriv_x + laplace_deriv_y;

                const T x = this->lbc_x_coordinate + i * this->grid_spatial_step_x;
                const T y = this->lbc_y_coordinate + j * this->grid_spatial_step_y;

                rhs_ode_system[i][j][grid_variable_displacement_number] =
                    grid_solution[i][j][grid_variable_velocity_number];

                rhs_ode_system[i][j][grid_variable_velocity_number] =
                    (this->wave_speed * this->wave_speed) * laplace_deriv + this->source(x, y, t);

            }
        }

        return rhs_ode_system;
    }

    void set_formation_dirichlet_boundary_conditions
        (Regular_grid_2d<Grid_variable>& grid_solution, const T& t, const T& dt) const noexcept
    {
        if constexpr (lbc_x == Boundary_conditions::Dirichlet) {
            const int lbc_node_x = 0;
            for (int i = 0; i < this->grid_nodes_amount_y; ++i) {
                const T y = this->lbc_y_coordinate + i * this->grid_spatial_step_y;
                grid_solution[lbc_node_x][i][grid_variable_displacement_number] =
                    this->function_lbc_x(y, t + dt);
                grid_solution[lbc_node_x][i][grid_variable_velocity_number] = T(0.0);
            }
        }

        if constexpr (rbc_x == Boundary_conditions::Dirichlet) {
            const int rbc_node_x = this->grid_steps_amount_x;
            for (int i = 0; i < this->grid_nodes_amount_y; ++i) {
                const T y = this->lbc_y_coordinate + i * this->grid_spatial_step_y;
                grid_solution[rbc_node_x][i][grid_variable_displacement_number] =
                    this->function_rbc_x(y, t + dt);
                grid_solution[rbc_node_x][i][grid_variable_velocity_number] = T(0.0);
            }
        }

        if constexpr (lbc_y == Boundary_conditions::Dirichlet) {
            const int lbc_node_y = 0;
            for (int i = 0; i < this->grid_nodes_amount_x; ++i) {
                const T x = this->lbc_x_coordinate + i * this->grid_spatial_step_x;
                grid_solution[i][lbc_node_y][grid_variable_displacement_number] =
                    this->function_lbc_y(x, t + dt);
                grid_solution[i][lbc_node_y][grid_variable_velocity_number] = T(0.0);
            }
        }

        if constexpr (rbc_y == Boundary_conditions::Dirichlet) {
            const int rbc_node_y = this->grid_steps_amount_y;
            for (int i = 0; i < this->grid_nodes_amount_x; ++i) {
                const T x = this->lbc_x_coordinate + i * this->grid_spatial_step_x;
                grid_solution[i][rbc_node_y][grid_variable_displacement_number] =
                    this->function_lbc_y(x, t + dt);
                grid_solution[i][rbc_node_y][grid_variable_velocity_number] = T(0.0);
            }
        }
    }

    void insert_boundary_conditions(Regular_grid_2d<Grid_variable>& rhs_ode_system,
                                    const Regular_grid_2d<Grid_variable>& grid_solution,
                                    const T& t) const noexcept
    {
        insert_lbc_x(rhs_ode_system, grid_solution, t);
        insert_rbc_x(rhs_ode_system, grid_solution, t);
        insert_lbc_y(rhs_ode_system, grid_solution, t);
        insert_rbc_y(rhs_ode_system, grid_solution, t);
    }

    void insert_lbc_x(Regular_grid_2d<Grid_variable>& rhs_ode_system,
                      const Regular_grid_2d<Grid_variable>& grid_solution,
                      const T& t) const noexcept
    {
        if constexpr (lbc_x == Boundary_conditions::Dirichlet) {
            const int lbc_node_x = 0;
            for (int i = 0; i < this->grid_nodes_amount_y; ++i) {
                rhs_ode_system[lbc_node_x][i][grid_variable_displacement_number] = T(0.0);
                rhs_ode_system[lbc_node_x][i][grid_variable_velocity_number] = T(0.0);
            }
        } else if constexpr (lbc_x == Boundary_conditions::Neumann) {
            const int lbc_node_x = 0;
            for (int i = 0; i < this->grid_nodes_amount_y; ++i) {
                const T y = this->lbc_y_coordinate + i * this->grid_spatial_step_y;
                const T h_sqr_x = this->grid_spatial_step_x * this->grid_spatial_step_x;
                const T multiplier =
                    grid_solution[1][i][grid_variable_displacement_number] -
                    grid_solution[0][i][grid_variable_displacement_number] -
                    this->grid_spatial_step_x * this->function_lbc_x(y, t);
                rhs_ode_system[lbc_node_x][i][grid_variable_displacement_number] =
                    grid_solution[lbc_node_x][i][grid_variable_velocity_number];
                rhs_ode_system[lbc_node_x][i][grid_variable_velocity_number] =
                    2 * (this->wave_speed * this->wave_speed) * multiplier / h_sqr_x +
                    this->source(this->lbc_x_coordinate, y, t);
            }
        } else if constexpr (lbc_x == Boundary_conditions::Robin) {
            const int lbc_node_x = 0;
            for (int i = 0; i < this->grid_nodes_amount_y; ++i) {
                const T y = this->lbc_y_coordinate + i * this->grid_spatial_step_y;
                const T h_sqr_x = this->grid_spatial_step_x * this->grid_spatial_step_x;
                constexpr int c_number = 0;
                constexpr int m_number = 1;
                const T multiplier =
                    grid_solution[1][i][grid_variable_displacement_number] -
                    grid_solution[0][i][grid_variable_displacement_number] -
                    this->grid_spatial_step_x * (std::get<m_number>(this->function_lbc_x)(y, t) -
                                                 std::get<c_number>(this->function_lbc_x)(y, t) *
                                                 grid_solution[0][i]);
                rhs_ode_system[lbc_node_x][i][grid_variable_displacement_number] =
                    grid_solution[lbc_node_x][i][grid_variable_velocity_number];
                rhs_ode_system[lbc_node_x][i][grid_variable_velocity_number] =
                    2 * (this->wave_speed * this->wave_speed) * multiplier / h_sqr_x +
                    this->source(this->lbc_x_coordinate, y, t);
            }
        } else { }
    }

    void insert_rbc_x(Regular_grid_2d<Grid_variable>& rhs_ode_system,
                      const Regular_grid_2d<Grid_variable>& grid_solution,
                      const T& t) const noexcept
    {
        if constexpr (rbc_x == Boundary_conditions::Dirichlet) {
            const int rbc_node_x = this->grid_steps_amount_x;
            for (int i = 0; i < this->grid_nodes_amount_y; ++i) {
                rhs_ode_system[rbc_node_x][i][grid_variable_displacement_number] = T(0.0);
                rhs_ode_system[rbc_node_x][i][grid_variable_velocity_number] = T(0.0);
            }
        } else if constexpr (rbc_x == Boundary_conditions::Neumann) {
            const int rbc_node_x = this->grid_steps_amount_x;
            for (int i = 0; i < this->grid_nodes_amount_y; ++i) {
                const T y = this->lbc_y_coordinate + i * this->grid_spatial_step_y;
                const T h_sqr_x = this->grid_spatial_step_x * this->grid_spatial_step_x;
                const T multiplier =
                    grid_solution[this->grid_steps_amount_x-1][i][grid_variable_displacement_number] -
                    grid_solution[this->grid_steps_amount_x][i][grid_variable_displacement_number] +
                    this->grid_spatial_step_x * this->function_rbc_x(y, t);
                rhs_ode_system[rbc_node_x][i][grid_variable_displacement_number] =
                    grid_solution[rbc_node_x][i][grid_variable_velocity_number];
                rhs_ode_system[rbc_node_x][i][grid_variable_velocity_number] =
                    2 * (this->wave_speed * this->wave_speed) * multiplier / h_sqr_x +
                    this->source(this->rbc_x_coordinate, y, t);
            }
        } else if constexpr (rbc_x == Boundary_conditions::Robin) {
            const int rbc_node_x = this->grid_steps_amount_x;
            for (int i = 0; i < this->grid_nodes_amount_y; ++i) {
                const T y = this->lbc_y_coordinate + i * this->grid_spatial_step_y;
                const T h_sqr_x = this->grid_spatial_step_x * this->grid_spatial_step_x;
                constexpr int c_number = 0;
                constexpr int m_number = 1;
                const T multiplier =
                    grid_solution[this->grid_steps_amount_x-1][i][grid_variable_displacement_number] -
                    grid_solution[this->grid_steps_amount_x][i][grid_variable_displacement_number] +
                    this->grid_spatial_step_x * (std::get<m_number>(this->function_rbc_x)(y, t) -
                                                 std::get<c_number>(this->function_rbc_x)(y, t) *
                                                 grid_solution[this->grid_steps_amount_x][i]);
                rhs_ode_system[rbc_node_x][i][grid_variable_displacement_number] =
                    grid_solution[rbc_node_x][i][grid_variable_velocity_number];
                rhs_ode_system[rbc_node_x][i][grid_variable_velocity_number] =
                    2 * (this->wave_speed * this->wave_speed) * multiplier / h_sqr_x +
                    this->source(this->rbc_x_coordinate, y, t);
            }
        } else { }
    }

    void insert_lbc_y(Regular_grid_2d<Grid_variable>& rhs_ode_system,
                      const Regular_grid_2d<Grid_variable>& grid_solution,
                      const T& t) const noexcept
    {
        if constexpr (lbc_y == Boundary_conditions::Dirichlet) {
            const int lbc_node_y = 0;
            for (int i = 0; i < this->grid_nodes_amount_x; ++i) {
                rhs_ode_system[i][lbc_node_y][grid_variable_displacement_number] = T(0.0);
                rhs_ode_system[i][lbc_node_y][grid_variable_velocity_number] = T(0.0);
            }
        } else if constexpr (lbc_y == Boundary_conditions::Neumann) {
            const int lbc_node_y = 0;
            for (int i = 0; i < this->grid_nodes_amount_x; ++i) {
                const T x = this->lbc_x_coordinate + i * this->grid_spatial_step_x;
                const T h_sqr_y = this->grid_spatial_step_y * this->grid_spatial_step_y;
                const T multiplier =
                    grid_solution[i][1][grid_variable_displacement_number] -
                    grid_solution[i][0][grid_variable_displacement_number] -
                    this->grid_spatial_step_y * this->function_lbc_y(x, t);
                rhs_ode_system[i][lbc_node_y][grid_variable_displacement_number] =
                    grid_solution[i][lbc_node_y][grid_variable_velocity_number];
                rhs_ode_system[i][lbc_node_y][grid_variable_velocity_number] =
                    2 * (this->wave_speed * this->wave_speed) * multiplier / h_sqr_y +
                    this->source(x, this->lbc_y_coordinate, t);
            }
        } else if constexpr (lbc_y == Boundary_conditions::Robin) {
            const int lbc_node_y = 0;
            for (int i = 0; i < this->grid_nodes_amount_x; ++i) {
                const T x = this->lbc_x_coordinate + i * this->grid_spatial_step_x;
                const T h_sqr_y = this->grid_spatial_step_y * this->grid_spatial_step_y;
                constexpr int c_number = 0;
                constexpr int m_number = 1;
                const T multiplier =
                    grid_solution[i][1][grid_variable_displacement_number] -
                    grid_solution[i][0][grid_variable_displacement_number] -
                    this->grid_spatial_step_y * (std::get<m_number>(this->function_lbc_y)(x, t) -
                                                 std::get<c_number>(this->function_lbc_y)(x, t) *
                                                 grid_solution[i][0]);
                rhs_ode_system[i][lbc_node_y][grid_variable_displacement_number] =
                    grid_solution[i][lbc_node_y][grid_variable_velocity_number];
                rhs_ode_system[i][lbc_node_y][grid_variable_velocity_number] =
                    2 * (this->wave_speed * this->wave_speed) * multiplier / h_sqr_y +
                    this->source(x, this->lbc_y_coordinate, t);
            }
        } else { }
    }

    void insert_rbc_y(Regular_grid_2d<Grid_variable>& rhs_ode_system,
                      const Regular_grid_2d<Grid_variable>& grid_solution,
                      const T& t) const noexcept
    {
        if constexpr (rbc_y == Boundary_conditions::Dirichlet) {
            const int rbc_node_y = this->grid_steps_amount_y;
            for (int i = 0; i < this->grid_nodes_amount_x; ++i) {
                rhs_ode_system[i][rbc_node_y][grid_variable_displacement_number] = T(0.0);
                rhs_ode_system[i][rbc_node_y][grid_variable_velocity_number] = T(0.0);
            }
        } else if constexpr (rbc_y == Boundary_conditions::Neumann) {
            const int rbc_node_y = this->grid_steps_amount_y;
            for (int i = 0; i < this->grid_nodes_amount_x; ++i) {
                const T x = this->lbc_x_coordinate + i * this->grid_spatial_step_x;
                const T h_sqr_y = this->grid_spatial_step_y * this->grid_spatial_step_y;
                const T multiplier =
                    grid_solution[i][this->grid_steps_amount_y-1][grid_variable_displacement_number] -
                    grid_solution[i][this->grid_steps_amount_y][grid_variable_displacement_number] +
                    this->grid_spatial_step_y * this->function_rbc_y(x, t);
                rhs_ode_system[i][rbc_node_y][grid_variable_displacement_number] =
                    grid_solution[i][rbc_node_y][grid_variable_velocity_number];
                rhs_ode_system[i][rbc_node_y][grid_variable_velocity_number] =
                    2 * (this->wave_speed * this->wave_speed) * multiplier / h_sqr_y +
                    this->source(x, this->rbc_y_coordinate, t);
            }
        } else if constexpr (rbc_y == Boundary_conditions::Robin) {
            const int rbc_node_y = this->grid_steps_amount_y;
            for (int i = 0; i < this->grid_nodes_amount_x; ++i) {
                const T x = this->lbc_x_coordinate + i * this->grid_spatial_step_x;
                const T h_sqr_y = this->grid_spatial_step_y * this->grid_spatial_step_y;
                constexpr int c_number = 0;
                constexpr int m_number = 1;
                const T multiplier =
                    grid_solution[i][this->grid_steps_amount_y-1][grid_variable_displacement_number] -
                    grid_solution[i][this->grid_steps_amount_y][grid_variable_displacement_number] +
                    this->grid_spatial_step_y * (std::get<m_number>(this->function_rbc_y)(x, t) -
                                                 std::get<c_number>(this->function_rbc_y)(x, t) *
                                                 grid_solution[i][this->grid_steps_amount_y]);
                rhs_ode_system[i][rbc_node_y][grid_variable_displacement_number] =
                    grid_solution[i][rbc_node_y][grid_variable_velocity_number];
                rhs_ode_system[i][rbc_node_y][grid_variable_velocity_number] =
                    2 * (this->wave_speed * this->wave_speed) * multiplier / h_sqr_y +
                    this->source(x, this->rbc_y_coordinate, t);
            }
        } else { }
    }
};

}

#endif // WAVE_EQUATION_2D_HXX_INCLUDED
