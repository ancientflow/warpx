#ifndef WARPX_INSERT_HALLELECTROSTATICCONFIG_H_
#define WARPX_INSERT_HALLELECTROSTATICCONFIG_H_

/* Variable-permittivity Poisson solve with a volume-anode overset mask. */
#if defined(WARPX_DIM_3D)
#define WARPX_USE_HALL_ELECTROSTATIC_MATERIALS

/* Independent far-field material solve and fine-grid Dirichlet coupling. */
#define WARPX_USE_HALL_FAR_FIELD
#endif

#if defined(WARPX_USE_HALL_FAR_FIELD) && !defined(WARPX_USE_HALL_ELECTROSTATIC_MATERIALS)
#error "WARPX_USE_HALL_FAR_FIELD requires the 3D Hall electrostatic material path"
#endif

#endif // WARPX_INSERT_HALLELECTROSTATICCONFIG_H_
