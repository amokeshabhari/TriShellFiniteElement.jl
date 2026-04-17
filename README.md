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

### Static condensation

The 3 mid-edge DOF (N₄–N₆) are eliminated numerically per element, reducing the 12×12 bending+shear partition to 9×9. Analytical condensation via symbolic tools yields correct but extremely long expressions; numerical condensation is more practical.

### Gaussian quadrature

| Partition | Rule | Reason |
|---|---|---|
| Membrane, geometric, shear | 1-point | Linear integrands; exact |
| Bending | 3-point | Quadratic N₄–N₆ require higher-order rule |

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

### Out-of-plane bending — midspan line load

Analytical: 1.289 mm (thick, t = 10 mm) | 156.252 mm (thin, t = 1 mm)

### Out-of-plane bending — uniform pressure

Analytical: 0.80075 mm (thick) | 97.657 mm (thin)

### In-plane bending — midspan line load

Analytical: 1.289 mm (thick) | 64.45 mm (thin)

### In-plane bending — uniform pressure

Analytical: 0.80075 mm (thick) | 40.0375 mm (thin)

### Column buckling (linear)

Euler critical stress (with shear): 1603.78 N/mm² (thick) | 0.65797 N/mm² (thin)

### Steel column base plate with bolt holes

Base plate meshed with [Gmsh.jl](https://github.com/JuliaFEM/Gmsh.jl); fixed BCs at bolt holes; tensile axial load applied at column stub.

---

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/runtosolve/TriShellFiniteElement.jl")
```

---

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
