# TriShellFiniteElement.jl

`TriShellFiniteElement.jl` is an open-source Julia package that implements a **3-node triangular Mindlin shell finite element (MIN3)** for linear static analysis and linear buckling of thin-walled structures. The package is built on top of [Ferrite.jl](https://ferrite-fem.github.io/Ferrite.jl/) and is documented in:

> Moen, C.D., Ádány, S., and Shabhari, A. (2026). *An open-source triangular shell finite element software implementation for the study of thin-walled structures.* Annual Stability Conference, Structural Stability Research Council (SSRC), Atlanta, Georgia, April 21–24, 2026.

The element formulation follows:

> Tessler, A. and Hughes, T.J.R. (1985). *A three-node Mindlin plate element with improved transverse shear.* Computer Methods in Applied Mechanics and Engineering, 50(1), 71–101.

---

## Key Features

- 3-node triangular shell element with **6 DOF per node** (u, v, w, θₓ, θᵧ, θᵤ), totalling **18 DOF per element**
- Composite element formulation:
  - In-plane **membrane** stiffness (plane stress, linear shape functions)
  - Out-of-plane **plate bending** stiffness (Mindlin-Reissner, anisoparametric interpolation)
  - **Transverse shear** stiffness with 5/6 shear correction factor
  - **Drilling DOF** stabilisation (in-plane rotation)
- **Geometric stiffness** matrix for linear buckling analysis (Green-Lagrange nonlinear strains)
- Local-to-global coordinate transformation for arbitrarily oriented elements in 3D space
- Validated against Abaqus S3R/S4R elements and closed-form analytical solutions

---

## Theoretical Background

### Flat Shell = Membrane + Plate

A *flat shell element* superposes two independent mechanical contributions:

1. **Membrane (in-plane)**: 2D plane-stress behaviour governing in-plane translations u and v.
2. **Plate (out-of-plane bending + transverse shear)**: governing out-of-plane translation w and bending rotations ψₓ, ψᵧ.

In the local coordinate system, membrane and bending are decoupled.  They are coupled only through the local-to-global transformation when the structure contains non-coplanar elements.

### Mindlin-Reissner Plate Theory

The package adopts **Mindlin-Reissner (first-order shear deformation) plate theory**.  Unlike Kirchhoff theory, cross-sections remain plane but not necessarily perpendicular to the mid-surface—this is the plate analogue of Timoshenko beam theory.

| Theory | Cross-section assumption | Applicable range |
|---|---|---|
| Kirchhoff (Euler-Bernoulli analogue) | Plane **and** perpendicular to mid-surface | Thin plates only |
| Mindlin (Timoshenko analogue) | Plane, but **not** necessarily perpendicular | Thin and moderately thick plates |

#### Governing Strains

```
Membrane strains (z = 0):
  εₓ,ₘ = ∂u/∂x
  εᵧ,ₘ = ∂v/∂y
  γₓᵧ,ₘ = ∂u/∂y + ∂v/∂x

Bending curvatures:
  κₓₓ = ∂ψₓ/∂x
  κᵧᵧ = ∂ψᵧ/∂y
  2κₓᵧ = ∂ψₓ/∂y + ∂ψᵧ/∂x

Transverse shear strains:
  γₓᵤ = ψₓ + ∂w/∂x
  γᵧᵤ = ψᵧ − ∂w/∂y
```

#### Constitutive Matrices

**Membrane** (plane stress, thickness t):

```
Dₘ = t / (1−ν²) × [E    νE   0        ]
                   [νE    E   0        ]
                   [0     0   (1−ν)E/2 ]
```

**Bending** (flexural rigidity D = Et³/[12(1−ν²)]):

```
Dᵦ = [D    νD   0         ]
     [νD    D   0         ]
     [0     0   D(1−ν)/2  ]
```

**Transverse shear** (5/6 shear correction factor, κ = 5/6):

```
Dₛ = [κGt   0  ]
     [0    κGt  ]
```

---

### The Shear Locking Problem and MIN3

A standard isoparametric linear triangle (same shape functions for w, ψₓ, ψᵧ) **locks** in the thin limit—the transverse shear energy grows as a penalty term and drives the element to excessively stiff (orders-of-magnitude wrong) solutions.

Tessler and Hughes (1985) proposed the **anisoparametric MIN3** element to resolve this via two auxiliary criteria:

- **Criterion 1**: Constant transverse shear strains must be representable for any element size.
- **Criterion 2**: The Kirchhoff (zero-shear) limit must be attainable without spurious constraining.

These criteria require the deflection w to use a polynomial one degree higher than the rotations ψₓ, ψᵧ.

![Triangular shell element](assets/moen2026_p2_img1.png)

*Figure: The 3-node triangular shell element with displacement functions u(x,y), v(x,y), w(x,y) and cross-section rotations ψₓ(x,y), ψᵧ(x,y) defined in the local element coordinate system x–y–z.*

#### Shape Functions

**`IP3`** — Linear 3-node functions for ψₓ, ψᵧ (and membrane u, v):

```
N₁ = 1 − ξ − η
N₂ = ξ
N₃ = η
```

**`IP6`** — Quadratic 6-node functions for w:

```
N₁ = 1 − ξ − η         (corner 1)
N₂ = ξ                  (corner 2)
N₃ = η                  (corner 3)
N₄ = 4ξ(1 − ξ − η)     (mid-edge, nodes 1–2)
N₅ = 4ξη               (mid-edge, nodes 2–3)
N₆ = 4η(1 − ξ − η)     (mid-edge, nodes 3–1)
```

![Shape functions and local coordinates](assets/moen2026_p4_img1.png)

*Figure: Local reference coordinates ξ–η and the six quadratic shape functions for w.  The three mid-edge nodes are eliminated by static condensation.*

#### Static Condensation

The 3 mid-edge DOF (nodes 4, 5, 6) are eliminated numerically, reducing the 12×12 bending+shear partition to a 9×9 system.  Although static condensation can be performed analytically with symbolic tools, the resulting expressions are extremely long; numerical condensation per-element is more efficient.

#### Gaussian Quadrature

The bending stiffness (IP6) requires at minimum a **3-point Gauss rule** to integrate the quadratic shape function products exactly over the triangle.  The membrane and geometric stiffness use the simpler 1-point rule.

![3-point Gaussian quadrature weights and locations](assets/moen2026_p6_img1.png)

*Figure: Gauss quadrature point coordinates (ξᵢ, ηᵢ) and weights Wᵢ for the 3-point triangular rule.*

| Point | ξ | η | Weight |
|---|---|---|---|
| 1 | 1/6 | 1/6 | 1/6 |
| 2 | 2/3 | 1/6 | 1/6 |
| 3 | 1/6 | 2/3 | 1/6 |

---

### Drilling DOF Stabilisation

Pure shell formulations have 5 DOF per node (u, v, w, ψₓ, ψᵧ).  In 3D assemblies, the in-plane rotational DOF θᵤ (about the element normal) is needed to prevent singularity where non-coplanar elements meet.  A small stabilisation stiffness equal to **min(diagonal bending rotational stiffness) / 100** is added for each θᵤ.

---

### Geometric Stiffness and Linear Buckling

The geometric stiffness is derived from the **Green-Lagrange nonlinear membrane strains** at z = 0:

```
εₓ,ₘᴺᴸ = ½[(∂u/∂x)² + (∂v/∂x)² + (∂w/∂x)²]
εᵧ,ₘᴺᴸ = ½[(∂u/∂y)² + (∂v/∂y)² + (∂w/∂y)²]
γₓᵧ,ₘᴺᴸ = (∂u/∂x)(∂u/∂y) + (∂v/∂x)(∂v/∂y) + (∂w/∂x)(∂w/∂y)
```

Because linear shape functions are used, the membrane stress field is **uniform** within each element, described by three scalars σₓ, σᵧ, τₓᵧ.  This allows an analytical geometric stiffness expression.

**Linear buckling** is solved from:

```
K φ = λ Kᵍ φ
```

where K is the elastic stiffness, Kᵍ is the geometric stiffness from a preceding linear static analysis, λ is the critical load multiplier, and φ is the buckling mode shape.

---

### Local-to-Global Transformation

Each element stiffness is computed in a local frame aligned to the element plane, then rotated to global X–Y–Z:

```
Kᵉ,global = Rᵀ Kᵉ,local R
```

The local frame is defined with origin at node 1, y′-axis along edge 1→2, z′-axis normal to the element plane (cross product of two edges), and x′-axis completing the right-handed triad.

---

## Element Formulation Summary

| Contribution | DOF | Matrix size | Quadrature | Notes |
|---|---|---|---|---|
| Membrane **Kₑ,ₘ** | u, v at 3 corners | 6×6 | 1-point | `IP3`, plane stress |
| Bending **Kₑ,ᵦ** | w, ψₓ, ψᵧ at 3 corners (+3 mid) | 12×12 → 9×9 | 3-point | `IP6` for w; static condensation |
| Shear **Kₑ,ₛ** | same 9 as bending | 9×9 | 1-point | 5/6 correction factor |
| Combined **Kₑ** | 5 DOF/node (no drilling) | 15×15 | — | Kₑ,ₘ ⊕ (Kₑ,ᵦ + Kₑ,ₛ) |
| + Drilling | 6 DOF/node | 18×18 | — | θᵤ stabilised |
| Geometric **Kᵍ** | u, v, w at 3 corners | 9×9 | 1-point | Constant σₓ, σᵧ, τₓᵧ; analytical |

---

## Validation Examples

All examples use a **100 mm × 1000 mm plate** with E = 200,000 MPa and ν = 0.30.  Seven mesh densities are tested (45 to 4,884,365 DOF).  Results are compared to Abaqus S4R/S3R elements and analytical solutions.

![Plate boundary conditions](assets/moen2026_p7_img1.png)

*Figure: Plate boundary conditions — ends fixed in Y, midspan in Z, two corner nodes in X.*

### Out-of-Plane Bending — Midspan Distributed Load

Analytical solutions: 1.289 mm (thick, t = 10 mm) and 156.252 mm (thin, t = 1 mm).

![Convergence: out-of-plane midspan load](assets/moen2026_p8_img1.png)
![Table 1](assets/moen2026_p8_img2.png)

### Out-of-Plane Bending — Uniform Pressure

Analytical solutions: 0.80075 mm (thick) and 97.657 mm (thin).

![Convergence: uniform pressure](assets/moen2026_p9_img1.png)
![Table 2](assets/moen2026_p9_img2.png)

### In-Plane Bending — Midspan Distributed Load

Analytical solutions: 1.289 mm (thick) and 64.45 mm (thin).

![Convergence: in-plane midspan load](assets/moen2026_p10_img1.png)
![Table 3](assets/moen2026_p10_img2.png)

### In-Plane Bending — Uniform Pressure

Analytical solutions: 0.80075 mm (thick) and 40.0375 mm (thin).

![Convergence: in-plane uniform pressure](assets/moen2026_p11_img1.png)
![Table 4](assets/moen2026_p11_img2.png)

### Column Buckling

Euler column buckling under uniform end compressive stress.  Analytical critical stresses (including shear): 1603.78 N/mm² (thick) and 0.65797 N/mm² (thin).

![Convergence: column buckling](assets/moen2026_p12_img1.png)
![Table 5](assets/moen2026_p12_img2.png)

### Steel Column Base Plate with Holes

A steel column base plate with bolt holes is meshed with Gmsh.jl and analysed under tensile axial load.

![Base plate mesh (undeformed)](assets/moen2026_p13_img1.png)
![Base plate deformation](assets/moen2026_p13_img2.png)

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

E = 200_000.0;  ν = 0.30;  t = 1.0

x = [Tensors.Vec((0.0, 0.0)),
     Tensors.Vec((100.0, 0.0)),
     Tensors.Vec((0.0, 100.0))]

ip3 = IP3();  ip6 = IP6()
qr1 = QuadratureRule{RefTriangle}(1)
qr3 = QuadratureRule{RefTriangle}(2)
cv3 = CellValues(qr1, ip3, ip3)
cv6 = CellValues(qr3, ip6, ip3)
reinit!(cv3, x);  reinit!(cv6, x)

Dm = calculate_membrane_constitutive_matrix(E, ν, t)
Db = calculate_bending_constitutive_matrix(E, ν, t)
Ds = calculate_shear_constitutive_matrix(E, ν, t)

Ke_m = calculate_element_membrane_stiffness_matrix(Dm, cv3)  # 6×6
Ke_b = calculate_element_bending_stiffness_matrix(Db, cv6)   # 12×12
Ke_s = calculate_element_shear_stiffness_matrix(Ds, cv6)     # 12×12
```

### Global Stiffness Assembly

```julia
using Ferrite, TriShellFiniteElement

E = 200_000.0;  ν = 0.30;  t = 1.0

grid = generate_grid(Triangle, (10, 10), Vec((0.0, 0.0, 0.0)), Vec((1000.0, 1000.0, 0.0)))
ip   = Lagrange{RefTriangle, 1}()
dh   = DofHandler(grid)
add!(dh, :u, ip^3)
add!(dh, :θ, ip^3)
close!(dh)

ip3 = IP3();  ip6 = IP6()
qr1 = QuadratureRule{RefTriangle}(1)
qr3 = QuadratureRule{RefTriangle}(2)

K  = allocate_matrix(dh)
K  = assemble_global_Ke!(K, dh, qr1, qr3, ip3, ip6, E, ν, t)
```

### Buckling Analysis

```julia
Kg = allocate_matrix(dh)
Kg = assemble_global_Kg!(Kg, dh, σXX, σYY, τXY)
# Solve: K φ = λ Kᵍ φ
```

---

## Dependencies

| Package | Role |
|---|---|
| [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl) | Mesh, DOF handler, cell values, quadrature |
| [Tensors.jl](https://github.com/Ferrite-FEM/Tensors.jl) | Coordinate vectors and tensor operations |
| LinearAlgebra | Matrix operations (Julia standard library) |

---

## References

1. Tessler, A. and Hughes, T.J.R. (1985). *A three-node Mindlin plate element with improved transverse shear.* Computer Methods in Applied Mechanics and Engineering, **50**(1), 71–101.
2. Moen, C.D., Ádány, S., and Shabhari, A. (2026). *An open-source triangular shell finite element software implementation for the study of thin-walled structures.* SSRC Annual Stability Conference, Atlanta, Georgia.
3. Liu, Y.J. and Riggs, H.R. (2002). *Development of the MIN-N family of triangular anisoparametric Mindlin plate elements.* UHM/CE/02-02, University of Hawaii at Manoa.
4. Dunavant, D.A. (1985). *High degree efficient symmetrical Gaussian quadrature rules for the triangle.* Int. J. Numer. Methods Eng., **21**, 1129–1148.
5. Moen, C.D. and Ádány, S. (2025). *Thin shell finite element formulations implemented in open-source software.* SSRC Annual Stability Conference, Louisville, KY.
