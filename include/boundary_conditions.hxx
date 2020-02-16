#ifndef BOUNDARY_CONDITIONS_HXX_INCLUDED
#define BOUNDARY_CONDITIONS_HXX_INCLUDED

namespace Math_structures
{

enum class Boundary_conditions: int
{
    Dirichlet = 0,
    Neumann = 1,
    Robin = 2,
};

}

#endif // BOUNDARY_CONDITIONS_HXX_INCLUDED
