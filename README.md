# TriShellFiniteElement.jl

`TriShellFiniteElement.jl` is a Julia package implementing the **3-node triangular Mindlin shell finite element (MIN3)** for linear static and buckling analysis of thin-walled structures. Built on [Ferrite.jl](https://ferrite-fem.github.io/Ferrite.jl/).

**Paper:** Moen, C.D., Ádány, S., and Shabhari, A. (2026). *An open-source triangular shell finite element software implementation for the study of thin-walled structures.* SSRC Annual Stability Conference, Atlanta, GA.

**Element theory:** Tessler, A. and Hughes, T.J.R. (1985). *A three-node Mindlin plate element with improved transverse shear.* Comput. Methods Appl. Mech. Eng., **50**(1), 71–101.

---

## Features

- 3-node triangular shell element — **18 DOF** (6 per node: u, v, w, θₓ, θᵧ, θᵤ)
- Anisoparametric interpolation (MIN3) to eliminate transverse shear locking
- Membrane, Mindlin-Reissner bending, and transverse shear stiffness (5/6 correction factor)
- Drilling DOF stabilisation for 3D assemblies
- Geometric stiffness matrix for linear buckling (Green-Lagrange nonlinear strains)
- Local-to-global transformation for arbitrarily oriented elements in 3D
- Validated against Abaqus S3R/S4R and analytical solutions

---

## Element Formulation

### Flat shell decomposition

The element superposes two decoupled contributions in the local frame:

| Sub-element | Active DOF | Governing physics |
|---|---|---|
| **Membrane** | u, v (in-plane translation) | 2D plane stress |
| **Plate** | w, ψₓ, ψᵧ (out-of-plane + rotations) | Mindlin-Reissner bending + transverse shear |

In-plane (membrane) and out-of-plane (bending/shear) are decoupled locally; they couple only through the local-to-global rotation when elements are non-coplanar.

### Mindlin-Reissner kinematics

Cross-sections remain plane but rotate independently of the deflection slope — the plate analogue of Timoshenko beam theory. The strain state is:

```
Membrane:       εₓ = ∂u/∂x,   εᵧ = ∂v/∂y,   γₓᵧ = ∂u/∂y + ∂v/∂x
Bending:        κₓₓ = ∂ψₓ/∂x,  κᵧᵧ = ∂ψᵧ/∂y,  2κₓᵧ = ∂ψₓ/∂y + ∂ψᵧ/∂x
Transv. shear:  γₓᵤ = ψₓ + ∂w/∂x,   γᵧᵤ = ψᵧ − ∂w/∂y
```

### Constitutive matrices

Let D = Et³/[12(1−ν²)] (flexural rigidity) and G = E/[2(1+ν)].

**Membrane** Dₘ (3×3, plane stress, thickness t):
```
    t      ⎡ 1    ν    0        ⎤
  ——————   ⎢ ν    1    0        ⎥
  (1−ν²)   ⎣ 0    0   (1−ν)/2  ⎦  × E
```

**Bending** Dᵦ (3×3):
```
  ⎡ D    νD   0        ⎤
  ⎢ νD    D   0        ⎥
  ⎣ 0     0  D(1−ν)/2  ⎦
```

**Transverse shear** Dₛ (2×2, 5/6 correction factor):
```
  ⎡ 5Gt/6   0    ⎤
  ⎣  0    5Gt/6  ⎦
```

### Shear locking and the anisoparametric (MIN3) fix

A standard isoparametric triangle (same linear shape functions for w, ψₓ, ψᵧ) locks in the thin limit — the shear energy penalty grows without bound. Tessler & Hughes (1985) resolve this with two kinematic criteria:

1. **Constant shear strains** must be representable at any element size.
2. **The Kirchhoff limit** (γ → 0) must be reachable without spurious constraining.

Both require w to be interpolated one polynomial degree higher than ψₓ, ψᵧ — hence the name *anisoparametric*.

### Shape functions

**`IP3` — linear, used for ψₓ, ψᵧ and membrane u, v:**

```
N₁ = 1 − ξ − η,   N₂ = ξ,   N₃ = η
```

**`IP6` — quadratic, used for w (6 nodes: 3 corner + 3 mid-edge):**

```
N₁ = 1 − ξ − η            N₄ = 4ξ(1 − ξ − η)
N₂ = ξ                     N₅ = 4ξη
N₃ = η                     N₆ = 4η(1 − ξ − η)
```

<<<<<<< HEAD
=======
![Local coordinates and shape functions](docs/src/assets/moen2026_p4_img1.png)

*Local reference coordinates ξ–η and the six shape functions for w. Mid-edge nodes N₄–N₆ are eliminated by static condensation, leaving 3 corner nodes × 3 DOF (w, ψₓ, ψᵧ).*

>>>>>>> 4705ff221a8ad3cf7ce4b603ecdfcfc523adf6bb
### Static condensation

The 3 mid-edge DOF (N₄–N₆) are eliminated numerically per element, reducing the 12×12 bending+shear partition to 9×9. Analytical condensation via symbolic tools yields correct but extremely long expressions; numerical condensation is more practical.

### Gaussian quadrature

| Partition | Rule | Reason |
|---|---|---|
| Membrane, geometric, shear | 1-point | Linear integrands; exact |
| Bending | 3-point | Quadratic N₄–N₆ require higher-order rule |

<<<<<<< HEAD
=======
![3-point Gauss rule for triangles](docs/src/assets/moen2026_p6_img1.png)

*3-point Gauss quadrature points and weights for the standard triangle with vertices (0,0), (1,0), (0,1): ξᵢ/ηᵢ ∈ {1/6, 2/3} and Wᵢ = 1/6.*

>>>>>>> 4705ff221a8ad3cf7ce4b603ecdfcfc523adf6bb
### Drilling DOF stabilisation

A 5-DOF-per-node element (u, v, w, ψₓ, ψᵧ) produces a singular stiffness when non-coplanar elements share a node. The in-plane rotational DOF θᵤ is added with stabilisation stiffness:

```
k_drill = min(diagonal bending-rotation entries of Kₑ) / 100
```

This prevents singularity without materially affecting deformation results.

### Geometric stiffness and buckling

The geometric stiffness **Kᵍ** is derived from second-order (Green-Lagrange) membrane strains at z = 0:

```
εₓᴺᴸ = ½[(∂u/∂x)² + (∂v/∂x)² + (∂w/∂x)²]
εᵧᴺᴸ = ½[(∂u/∂y)² + (∂v/∂y)² + (∂w/∂y)²]
γₓᵧᴺᴸ =  (∂u/∂x)(∂u/∂y) + (∂v/∂x)(∂v/∂y) + (∂w/∂x)(∂w/∂y)
```

Linear membrane shape functions → constant σₓ, σᵧ, τₓᵧ per element → analytical **Kᵍ**.

Linear buckling is posed as the generalised eigenvalue problem:

```
K φ = λ Kᵍ φ
```

where λ is the critical load multiplier and φ the buckling mode.

### Local-to-global transformation

The local frame is defined per element: origin at node 1, y′ along edge 1→2, z′ normal to the element plane (= edge 1→2 × edge 1→3 normalised), x′ completing the right-hand triad. The assembled element stiffness is transformed as:

```
Kᵉ_global = Rᵀ Kᵉ_local R
```

---

## Stiffness Matrix Summary

| Partition | DOF | Size | Quadrature |
|---|---|---|---|
| Membrane **Kₑ,ₘ** | u, v at 3 corners | 6×6 | 1-pt, `IP3` |
| Bending **Kₑ,ᵦ** | w, ψₓ, ψᵧ — 6 nodes → 3 corners | 12×12 → 9×9 | 3-pt, `IP6` |
| Shear **Kₑ,ₛ** | same 9 DOF as bending | 9×9 | 1-pt, `IP3`/`IP6` |
| Combined **Kₑ** | 5 DOF/node (no drilling) | 15×15 | — |
| + Drilling | 6 DOF/node | **18×18** | — |
| Geometric **Kᵍ** | u, v, w at 3 corners | 9×9 | 1-pt, analytical |

---

## Validation

All examples: 100 mm × 1000 mm isotropic plate, E = 200,000 MPa, ν = 0.30. Seven mesh densities (45 → 4,884,365 DOF). Benchmarked against Abaqus S4R/S3R and closed-form analytical solutions (Euler-Bernoulli / Timoshenko).

<<<<<<< HEAD
=======
![Boundary conditions](docs/src/assets/moen2026_p7_img1.png)

*Plate boundary conditions: ends fixed in Y, midspan pinned in Z, two corner nodes in X.*

>>>>>>> 4705ff221a8ad3cf7ce4b603ecdfcfc523adf6bb
### Out-of-plane bending — midspan line load

Analytical: 1.289 mm (thick, t = 10 mm) | 156.252 mm (thin, t = 1 mm)

<<<<<<< HEAD
=======
![](docs/src/assets/moen2026_p8_img1.png) ![](docs/src/assets/moen2026_p8_img2.png)

>>>>>>> 4705ff221a8ad3cf7ce4b603ecdfcfc523adf6bb
### Out-of-plane bending — uniform pressure

Analytical: 0.80075 mm (thick) | 97.657 mm (thin)

<<<<<<< HEAD
=======
![](docs/src/assets/moen2026_p9_img1.png) ![](docs/src/assets/moen2026_p9_img2.png)

>>>>>>> 4705ff221a8ad3cf7ce4b603ecdfcfc523adf6bb
### In-plane bending — midspan line load

Analytical: 1.289 mm (thick) | 64.45 mm (thin)

<<<<<<< HEAD
=======
![](docs/src/assets/moen2026_p10_img1.png) ![](docs/src/assets/moen2026_p10_img2.png)

>>>>>>> 4705ff221a8ad3cf7ce4b603ecdfcfc523adf6bb
### In-plane bending — uniform pressure

Analytical: 0.80075 mm (thick) | 40.0375 mm (thin)

<<<<<<< HEAD
=======
![](docs/src/assets/moen2026_p11_img1.png) ![](docs/src/assets/moen2026_p11_img2.png)

>>>>>>> 4705ff221a8ad3cf7ce4b603ecdfcfc523adf6bb
### Column buckling (linear)

Euler critical stress (with shear): 1603.78 N/mm² (thick) | 0.65797 N/mm² (thin)

<<<<<<< HEAD
=======
![](docs/src/assets/moen2026_p12_img1.png) ![](docs/src/assets/moen2026_p12_img2.png)

>>>>>>> 4705ff221a8ad3cf7ce4b603ecdfcfc523adf6bb
### Steel column base plate with bolt holes

Base plate meshed with [Gmsh.jl](https://github.com/JuliaFEM/Gmsh.jl); fixed BCs at bolt holes; tensile axial load applied at column stub.

<<<<<<< HEAD
=======
![](docs/src/assets/moen2026_p13_img1.png) ![](docs/src/assets/moen2026_p13_img2.png)

*Left: undeformed mesh. Right: out-of-plane deformation contour.*

>>>>>>> 4705ff221a8ad3cf7ce4b603ecdfcfc523adf6bb
---

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/runtosolve/TriShellFiniteElement.jl")
```

---

<<<<<<< HEAD
=======
## Usage

### Element-level stiffness

```julia
using Ferrite, Tensors, TriShellFiniteElement

E = 200_000.0; ν = 0.30; t = 1.0   # MPa, —, mm

# Node coordinates in the local 2D plane
x = [Tensors.Vec((0.0, 0.0)), Tensors.Vec((100.0, 0.0)), Tensors.Vec((0.0, 100.0))]

ip3 = TriShellFiniteElement.IP3()
ip6 = TriShellFiniteElement.IP6()
qr1 = QuadratureRule{RefTriangle}(1)   # 1-point (membrane / shear)
qr3 = QuadratureRule{RefTriangle}(2)   # 3-point (bending — IP6 requires this)
cv3 = CellValues(qr1, ip3, ip3);  reinit!(cv3, x)
cv6 = CellValues(qr3, ip6, ip3);  reinit!(cv6, x)

Dm = TriShellFiniteElement.calculate_membrane_constitutive_matrix(E, ν, t)  # 3×3
Db = TriShellFiniteElement.calculate_bending_constitutive_matrix(E, ν, t)   # 3×3
Ds = TriShellFiniteElement.calculate_shear_constitutive_matrix(E, ν, t)     # 2×2

Ke_m = TriShellFiniteElement.calculate_element_membrane_stiffness_matrix(Dm, cv3)  # 6×6
Ke_b = TriShellFiniteElement.calculate_element_bending_stiffness_matrix(Db, cv6)   # 12×12
Ke_s = TriShellFiniteElement.calculate_element_shear_stiffness_matrix(Ds, cv6)     # 12×12
```

### Full 18×18 local shell stiffness

```julia
x3d = [Vec((0.0,0.0,0.0)), Vec((100.0,0.0,0.0)), Vec((0.0,100.0,0.0))]
Ke  = zeros(18, 18)
TriShellFiniteElement.local_elastic_stiffness_matrix!(Ke, x3d, E, ν, t)
```

### Global assembly and buckling

```julia
using Ferrite, TriShellFiniteElement

grid = generate_grid(Triangle, (10,10), Vec((0.0,0.0,0.0)), Vec((1000.0,1000.0,0.0)))
dh   = DofHandler(grid)
ip   = Lagrange{RefTriangle,1}()
add!(dh, :u, ip^3)   # u, v, w
add!(dh, :θ, ip^3)   # θₓ, θᵧ, θᵤ
close!(dh)

ip3 = TriShellFiniteElement.IP3();  ip6 = TriShellFiniteElement.IP6()
qr1 = QuadratureRule{RefTriangle}(1);  qr3 = QuadratureRule{RefTriangle}(2)

K  = allocate_matrix(dh)
K  = TriShellFiniteElement.assemble_global_Ke!(K, dh, qr1, qr3, ip3, ip6, E, ν, t)

# After linear static solve for σXX, σYY, τXY:
Kg = allocate_matrix(dh)
Kg = TriShellFiniteElement.assemble_global_Kg!(Kg, dh, σXX, σYY, τXY)

# Solve: K φ = λ Kᵍ φ  (use ArnoldiMethod.jl or KrylovKit.jl)
```

---

>>>>>>> 4705ff221a8ad3cf7ce4b603ecdfcfc523adf6bb
## Dependencies

| Package | Role |
|---|---|
| [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl) | Mesh, DOF handler, cell values, quadrature |
| [Tensors.jl](https://github.com/Ferrite-FEM/Tensors.jl) | Coordinate vectors |
| LinearAlgebra | Standard library — matrix operations |

---

## References

1. Tessler, A. and Hughes, T.J.R. (1985). *A three-node Mindlin plate element with improved transverse shear.* Comput. Methods Appl. Mech. Eng., **50**(1), 71–101.
2. Moen, C.D., Ádány, S., and Shabhari, A. (2026). *An open-source triangular shell finite element software implementation for the study of thin-walled structures.* SSRC Annual Stability Conference, Atlanta, GA.
3. Liu, Y.J. and Riggs, H.R. (2002). *Development of the MIN-N family of triangular anisoparametric Mindlin plate elements.* Research Report UHM/CE/02-02, Univ. of Hawaii at Manoa.
4. Dunavant, D.A. (1985). *High degree efficient symmetrical Gaussian quadrature rules for the triangle.* Int. J. Numer. Methods Eng., **21**, 1129–1148.
5. Moen, C.D. and Ádány, S. (2025). *Thin shell finite element formulations implemented in open-source software.* SSRC Annual Stability Conference, Louisville, KY.
