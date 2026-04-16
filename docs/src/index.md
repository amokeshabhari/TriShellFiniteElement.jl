# TriShellFiniteElement.jl

`TriShellFiniteElement.jl` is a Julia package that implements a **3-node triangular shell finite element** with drilling degrees of freedom (DOF). The element combines membrane, plate bending, and transverse shear stiffness into an 18-DOF shell element suitable for linear static and buckling analysis of shell structures.

---

## Key Features

- 3-node triangular shell element with **6 DOF per node** (u, v, w, θₓ, θᵧ, θᵤ)
- Composite element formulation:
  - In-plane **membrane** stiffness (plane stress)
  - Out-of-plane **plate bending** stiffness (Mindlin-Reissner)
  - **Transverse shear** stiffness with 5/6 shear correction factor
  - **Drilling DOF** stabilisation
- **Geometric stiffness** matrix for linear buckling analysis
- Local-to-global coordinate transformation for arbitrarily oriented elements in 3D space
- Integration with [Ferrite.jl](https://ferrite-fem.github.io/Ferrite.jl/) for mesh handling and DOF management

---

## Element Formulation

The element has **18 DOF total** (3 nodes × 6 DOF/node):

```
Node DOF: [u, v, w, θₓ, θᵧ, θᵤ]
```

| Symbol | Description |
|--------|-------------|
| u, v | In-plane translations (membrane) |
| w | Out-of-plane translation (bending/shear) |
| θₓ, θᵧ | Bending rotations |
| θᵤ | Drilling rotation (in-plane) |

The local stiffness is built from three contributions:

| Contribution | Size | Quadrature | Notes |
|---|---|---|---|
| Membrane | 6×6 | 1-point | Plane stress, `IP3` shape functions |
| Bending | 12×12 | 1-point | Mindlin-Reissner, `IP3` |
| Shear | 18×18 | 3-point | `IP6` shape functions, static condensation applied |

The combined bending+shear matrix uses a shear correction factor `α` that balances rotational stiffness between the two contributions. Drilling DOF are stabilised at 1/100 of the minimum bending diagonal.

---

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/runtosolve/TriShellFiniteElement.jl")
```

---

## Basic Usage

### Element-Level Stiffness Matrices

```julia
using Ferrite, Tensors, TriShellFiniteElement

E = 200_000.0   # Young's modulus (MPa)
ν = 0.30        # Poisson's ratio
t = 1.0         # Shell thickness (mm)

# Element node coordinates in local 2D system
x = [Tensors.Vec((0.0, 0.0)),
     Tensors.Vec((100.0, 0.0)),
     Tensors.Vec((0.0, 100.0))]

ip3 = IP3()
qr1 = QuadratureRule{RefTriangle}(1)
cv  = CellValues(qr1, ip3, ip3)
reinit!(cv, x)

# Constitutive matrices
Dm = calculate_membrane_constitutive_matrix(E, ν, t)
Db = calculate_bending_constitutive_matrix(E, ν, t)
Ds = calculate_shear_constitutive_matrix(E, ν, t)

# Element stiffness matrices
Ke_m = calculate_element_membrane_stiffness_matrix(Dm, cv)  # 6×6
Ke_b = calculate_element_bending_stiffness_matrix(Db, cv)   # 18×18
Ke_s = calculate_element_shear_stiffness_matrix(Ds, cv)     # 18×18
```

### Global Stiffness Assembly

```julia
using Ferrite, TriShellFiniteElement

E = 200_000.0
ν = 0.30
t = 1.0

# Build mesh and DOF handler (3D grid of triangles)
grid = generate_grid(Triangle, (10, 10), Vec((0.0, 0.0, 0.0)), Vec((1000.0, 1000.0, 0.0)))
ip   = Lagrange{RefTriangle, 1}()
dh   = DofHandler(grid)
add!(dh, :u, ip^3)   # u, v, w translations
add!(dh, :θ, ip^3)   # θₓ, θᵧ, θᵤ rotations
close!(dh)

# Quadrature rules and custom interpolations
ip3 = IP3()
ip6 = IP6()
qr1 = QuadratureRule{RefTriangle}(1)
qr3 = QuadratureRule{RefTriangle}(2)

# Assemble global elastic stiffness
K = allocate_matrix(dh)
K = assemble_global_Ke!(K, dh, qr1, qr3, ip3, ip6, E, ν, t)
```

### Buckling Analysis

After computing element stresses `σXX`, `σYY`, `τXY` from a linear static solve:

```julia
# Assemble geometric stiffness
Kg = allocate_matrix(dh)
Kg = assemble_global_Kg!(Kg, dh, σXX, σYY, τXY)

# Solve generalised eigenvalue problem  K φ = λ Kg φ
# using your preferred eigensolver (e.g. ArnoldiMethod.jl, KrylovKit.jl)
```

---

## Dependencies

| Package | Role |
|---|---|
| [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl) | Mesh, DOF handler, cell values, quadrature |
| [Tensors.jl](https://github.com/Ferrite-FEM/Tensors.jl) | Coordinate vectors and tensor operations |
| LinearAlgebra | Matrix operations (Julia standard library) |
