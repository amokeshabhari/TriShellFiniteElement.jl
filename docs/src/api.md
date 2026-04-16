# API Reference

## Interpolations

```@docs
IP3
IP6
```

## Constitutive Matrices

```@docs
calculate_membrane_constitutive_matrix
calculate_bending_constitutive_matrix
calculate_shear_constitutive_matrix
```

## Element Stiffness Matrices

```@docs
calculate_element_membrane_stiffness_matrix
calculate_element_bending_stiffness_matrix
calculate_element_shear_stiffness_matrix
local_elastic_stiffness_matrix!
```

## Geometric Stiffness

```@docs
calculate_element_geometric_stiffness_matrix
```

## Global Assembly

```@docs
assemble_global_Ke!
assemble_global_Kg!
```

## Coordinate Transformations

```@docs
calculation_rotation_matrix
global_nodal_coords_to_planar_coords
rotation_matrix_for_element_stiffness_drilling
```
