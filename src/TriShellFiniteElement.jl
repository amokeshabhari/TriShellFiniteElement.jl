module TriShellFiniteElement

using Ferrite, LinearAlgebra, Tensors

export IP3, IP6
export calculate_membrane_constitutive_matrix, calculate_bending_constitutive_matrix, calculate_shear_constitutive_matrix
export calculate_element_membrane_stiffness_matrix, calculate_element_bending_stiffness_matrix, calculate_element_shear_stiffness_matrix
export local_elastic_stiffness_matrix!
export calculate_element_geometric_stiffness_matrix
export assemble_global_Ke!, assemble_global_Kg!
export calculation_rotation_matrix, global_nodal_coords_to_planar_coords, rotation_matrix_for_element_stiffness_drilling

"""
    IP6 <: ScalarInterpolation{RefTriangle, 2}

6-node quadratic triangular interpolation on the reference triangle.

Shape functions (ξ₁, ξ₂ are the reference coordinates):
- N₁ = 1 − ξ₁ − ξ₂
- N₂ = ξ₁
- N₃ = ξ₂
- N₄ = 4ξ₁(1 − ξ₁ − ξ₂)
- N₅ = 4ξ₁ξ₂
- N₆ = 4ξ₂(1 − ξ₁ − ξ₂)

Used as the shape interpolation for the transverse shear term of the shell element.
"""
struct IP6 <: ScalarInterpolation{RefTriangle, 2}
end


function Ferrite.reference_shape_value(ip::IP6, ξ::Vec{2}, shape_number::Int)
    ξ₁ = ξ[1]
    ξ₂ = ξ[2]

    shape_number == 1 && return 1 - ξ₁ - ξ₂
    shape_number == 2 && return ξ₁
    shape_number == 3 && return ξ₂
    shape_number == 4 && return 4* ξ₁ * (1 - ξ₁ - ξ₂)
    shape_number == 5 && return 4 * ξ₁ * ξ₂
    shape_number == 6 && return 4 * ξ₂ * (1 - ξ₁ - ξ₂)

    throw(ArgumentError("no shape function $shape_number for interpolation $ip"))
end


Ferrite.getnbasefunctions(::IP6) = 6

Ferrite.adjust_dofs_during_distribution(::IP6) = false


"""
    IP3 <: ScalarInterpolation{RefTriangle, 2}

3-node linear triangular interpolation on the reference triangle.

Shape functions (ξ₁, ξ₂ are the reference coordinates):
- N₁ = 1 − ξ₁ − ξ₂
- N₂ = ξ₁
- N₃ = ξ₂

Used as both the geometry and shape interpolation for the membrane and bending terms of the shell element.
"""
struct IP3 <: ScalarInterpolation{RefTriangle, 2}
end


function Ferrite.reference_shape_value(ip::IP3, ξ::Vec{2}, shape_number::Int)
    ξ₁ = ξ[1]
    ξ₂ = ξ[2]

    shape_number == 1 && return 1 - ξ₁ - ξ₂
    shape_number == 2 && return ξ₁
    shape_number == 3 && return ξ₂

    throw(ArgumentError("no shape function $shape_number for interpolation $ip"))
end


# Ferrite.vertexdof_indices(::IP3) = ((1,2,3,4,5), (6,7,8,9,10), (11,12,13,14,15))
#organized by field u:1-9, θ:10-15
# Ferrite.vertexdof_indices(::IP3) = ((1,2,3,10,11), (4,5,6,12,13), (7,8,9,14,15))

Ferrite.getnbasefunctions(::IP3) = 3

Ferrite.adjust_dofs_during_distribution(::IP3) = false



function get_jacobian(ξ, ip, x)

    shape_number = 1
    dNdξ1 = Ferrite.reference_shape_gradient(ip, ξ, shape_number)

    shape_number = 2
    dNdξ2 = Ferrite.reference_shape_gradient(ip, ξ, shape_number)

    shape_number = 3
    dNdξ3 = Ferrite.reference_shape_gradient(ip, ξ, shape_number)


    J11 = [dNdξ1[1] dNdξ2[1] dNdξ3[1]] * [x[1][1], x[2][1], x[3][1]]
    J12 = [dNdξ1[1] dNdξ2[1] dNdξ3[1]] * [x[1][2], x[2][2], x[3][2]]
    J21 = [dNdξ1[2] dNdξ2[2] dNdξ3[2]] * [x[1][1], x[2][1], x[3][1]]
    J22 = [dNdξ1[2] dNdξ2[2] dNdξ3[2]] * [x[1][2], x[2][2], x[3][2]]


    J = [J11 J12
        J21 J22]

    return J

end


"""
    calculate_membrane_constitutive_matrix(E, ν, t)

Return the 3×3 plane-stress membrane constitutive matrix scaled by thickness `t`.

# Arguments
- `E`: Young's modulus
- `ν`: Poisson's ratio
- `t`: Shell thickness

# Returns
3×3 matrix `D` such that `{σ} = D {ε}` with `{ε} = [εₓ, εᵧ, γₓᵧ]ᵀ`.
"""
function calculate_membrane_constitutive_matrix(E, ν, t)

    G=E/(2*(1+ν))

    D=[    E/(1-ν^2) ν*E/(1-ν^2)    0
        ν*E/(1-ν^2)    E/(1-ν^2)    0
                0             0    G]

    D = D .* t

    return D

end




"""
    calculate_bending_constitutive_matrix(E, ν, t)

Return the 3×3 plate-bending constitutive matrix.

# Arguments
- `E`: Young's modulus
- `ν`: Poisson's ratio
- `t`: Shell thickness

# Returns
3×3 matrix `D` with bending stiffness constant `D₀ = Et³ / (12(1−ν²))`.
"""
function calculate_bending_constitutive_matrix(E, ν, t)

    D_const = E * t^3 / (12 * (1 - ν^2))
    D = D_const * [1.0  ν    0.0
                    ν    1.0  0.0
                    0.0  0.0  (1-ν)/2]

    return D

end



"""
    calculate_shear_constitutive_matrix(E, ν, t)

Return the 2×2 transverse shear constitutive matrix with a 5/6 shear correction factor.

# Arguments
- `E`: Young's modulus
- `ν`: Poisson's ratio
- `t`: Shell thickness

# Returns
2×2 diagonal matrix with entries `(5/6) G t` where `G = E / (2(1+ν))`.
"""
function calculate_shear_constitutive_matrix(E, ν, t)

    G=E/(2*(1+ν))

    D=[  5/6*G*t    0.0
        0.0        5/6*G*t]

    return D

end

"""
    calculate_element_membrane_stiffness_matrix(D, cv)

Return the 6×6 element membrane stiffness matrix in local coordinates.

# Arguments
- `D`: 3×3 membrane constitutive matrix (from [`calculate_membrane_constitutive_matrix`](@ref))
- `cv`: `CellValues` initialised with the element's local 2D node coordinates

# Returns
6×6 matrix `Kₑ = ∫ Bᵀ D B dA` over the element, where `B` is the 3×6 strain-displacement matrix.
"""
function calculate_element_membrane_stiffness_matrix(D, cv)

    ke = zeros(Float64, 6, 6)

    dNdx = cv.fun_values.dNdx


    for q_point in 1:getnquadpoints(cv)

        B_node_all = []

        for i in 1:3

            B_node = [dNdx[i][1]  0.0
            0.0         dNdx[i][2]
            dNdx[i][2]  dNdx[i][1]]

                            #  B_node = [dNdξ1*Jinv[1,1] + dNdξ2*Jinv[1,2]         0.0
                    #     0.0                                       dNdξ1*Jinv[2,1] + dNdξ2*Jinv[2,2]
                    #     dNdξ1*Jinv[2,1] + dNdξ2*Jinv[2,2]         dNdξ1*Jinv[1,1] + dNdξ2*Jinv[1,2]]



            push!(B_node_all, B_node)

        end

        B = hcat(B_node_all...)
        ke += B' * D * B .* getdetJdV(cv, q_point)

    end

    return ke

end



"""
    calculate_element_bending_stiffness_matrix(D, cv)

Return the 18×18 element plate-bending stiffness matrix in local coordinates.

The DOF ordering follows the Mindlin-Reissner convention: `[w, θₓ, θᵧ]` per node.
Zero rows/columns correspond to DOF not engaged by the bending strain-displacement matrix.

# Arguments
- `D`: 3×3 bending constitutive matrix (from [`calculate_bending_constitutive_matrix`](@ref))
- `cv`: `CellValues` initialised with the element's local 2D node coordinates

# Returns
18×18 matrix (rows/columns for drilling DOF are zero).
"""
function calculate_element_bending_stiffness_matrix(D, cv)

    # num_shape_functions = getnbasefunctions(ip_shape)

    ke = zeros(Float64, 18, 18)

     dNdx = cv.fun_values.dNdx

    for q_point in 1:getnquadpoints(cv)

        # ξ = qr.points[q_point]
        # J = get_jacobian(ξ, ip_geo, x)
        # Jinv = inv(J)

        B_node_all = []

        for i in 1:3

            # dNdξ1 = cv.fun_values.dNdξ[i + (q_point-1)*num_shape_functions][1]
            # dNdξ2 = cv.fun_values.dNdξ[i + (q_point-1)*num_shape_functions][2]

            B_node = [0.0       0.0             dNdx[i][1]
                        0.0       -dNdx[i][2]  0.0
                        0.0       -dNdx[i][1]  dNdx[i][2]]


            # B_node = [0.0       0.0             dNdξ1*Jinv[1,1] + dNdξ2*Jinv[1,2]
            #             0.0       -(dNdξ1*Jinv[2,1] + dNdξ2*Jinv[2,2])  0.0
            #             0.0       -(dNdξ1*Jinv[1,1] + dNdξ2*Jinv[1,2])  dNdξ1*Jinv[2,1] + dNdξ2*Jinv[2,2]]


            push!(B_node_all, B_node)

        end

        push!(B_node_all, zeros(3, 9))
        B = hcat(B_node_all...)

        println("B:", B)

        ke += B' * D * B .* getdetJdV(cv, q_point)

    end

    return ke

end




"""
    calculate_element_shear_stiffness_matrix(D, cv)

Return the 18×18 element transverse shear stiffness matrix in local coordinates.

The formulation uses a 6-node quadratic geometry (`IP6`) for the shape functions and
includes the drilling DOF rows/columns. Static condensation of the internal drilling DOF
is applied in [`local_elastic_stiffness_matrix!`](@ref).

# Arguments
- `D`: 2×2 shear constitutive matrix (from [`calculate_shear_constitutive_matrix`](@ref))
- `cv`: `CellValues` initialised with `IP6` shape functions and local 2D node coordinates

# Returns
18×18 matrix.
"""
function calculate_element_shear_stiffness_matrix(D, cv)

    num_shape_functions = getnbasefunctions(cv)

    ke = zeros(Float64, 18, 18)

    dNdx = cv.fun_values.dNdx

    println("dNdx:", dNdx)

    for q_point in 1:getnquadpoints(cv)

        # ξ = qr.points[q_point]
        # J = get_jacobian(ξ, ip_geo, x)
        # Jinv = inv(J)



        B_node_all = []

        for i=1:num_shape_functions

            # dNdξ1 = cv.fun_values.dNdξ[i + (q_point-1)*num_shape_functions][1]
            # dNdξ2 = cv.fun_values.dNdξ[i + (q_point-1)*num_shape_functions][2]

            B_node = [dNdx[i, q_point][1]        0.0     0.0
                      dNdx[i, q_point][2]         0.0     0.0]


            #   B_node = [dNdξ1*Jinv[1,1] + dNdξ2*Jinv[1,2]         0.0     0.0
            #                     dNdξ1*Jinv[2,1] + dNdξ2*Jinv[2,2]         0.0     0.0]

            if i <= 3

                N = Ferrite.shape_value(cv, q_point, i)


                    println("N " * string(i) * ":", N)


                # N = Ferrite.reference_shape_value(ip_shape, ξ, i)

                B_node += [0.0     0.0     N
                           0.0     -N     0.0]

            end

            push!(B_node_all, B_node)

        end

        B = hcat(B_node_all...)

        println("B_shear:", B)

        # println("B' * D * B ", B' * D * B)

        # ind=[1:10; 13; 16]
        # BDB = B' * D * B
        # BDB = BDB[ind, ind]

        ke += B' * D * B .* getdetJdV(cv, q_point)

        #  ke += BDB .* getdetJdV(cv, q_point)

        # println("B' * D * B .* getdetJdV(cv, q_point):", B' * D * B .* getdetJdV(cv, q_point))

        # println("getdetJdV(cv, q_point):", getdetJdV(cv, q_point))

    end

    return ke

end



"""
    local_elastic_stiffness_matrix!(qr1, qr3, ip3, ip6, E, ν, t, x)

Assemble the 18×18 local elastic stiffness matrix for a single shell element.

Combines membrane, plate-bending, and transverse shear contributions:
1. Membrane stiffness (1-point quadrature, `IP3`)
2. Bending stiffness (1-point quadrature, `IP3`)
3. Shear stiffness (3-point quadrature, `IP6`), with static condensation of drilling DOF
4. Shear correction factor applied to balance bending vs. shear rotational stiffness
5. Drilling DOF stabilisation (1/100 of minimum bending diagonal)

# Arguments
- `qr1`: 1-point `QuadratureRule{RefTriangle}` for membrane/bending
- `qr3`: 3-point `QuadratureRule{RefTriangle}` for shear
- `ip3`: [`IP3`](@ref) interpolation
- `ip6`: [`IP6`](@ref) interpolation
- `E`: Young's modulus
- `ν`: Poisson's ratio
- `t`: Shell thickness
- `x`: Vector of three 2D local node coordinates (`Tensors.Vec{2}`)

# Returns
18×18 local element stiffness matrix.
"""
function local_elastic_stiffness_matrix!(qr1, qr3, ip3, ip6, E, ν, t, x)

    #####membrane
    println("x:", x)

    cv = CellValues(qr1, ip3, ip3)
    reinit!(cv, x)

    Dm = TriShellFiniteElement.calculate_membrane_constitutive_matrix(E, ν, t)

    println("Dm:", Dm)

    # D = Dm
    # ip_geo = ip3
    # ip_shape = ip3
    # qr = qr1
    ke_m = TriShellFiniteElement.calculate_element_membrane_stiffness_matrix(Dm, cv)

    println("ke_m:", ke_m)

    ######bending
    cv = CellValues(qr1, ip3, ip3)
    reinit!(cv, x)

    Db = TriShellFiniteElement.calculate_bending_constitutive_matrix(E, ν, t)

    # D = Db
    # ip_geo = ip3
    # ip_shape = ip3
    # qr = qr1
    ke_b = TriShellFiniteElement.calculate_element_bending_stiffness_matrix(Db, cv)

    #remove zeros
    indices = [1:10; 13; 16]
    ke_b = ke_b[indices, indices]

     println("ke_b:", ke_b)

    ######shear
    cv = CellValues(qr3, ip6, ip3)
    reinit!(cv, x)

    Ds = TriShellFiniteElement.calculate_shear_constitutive_matrix(E, ν, t)

    # D = Ds
    # ip_geo = ip3
    # ip_shape = ip6
    # qr = qr3
    ke_s = TriShellFiniteElement.calculate_element_shear_stiffness_matrix(Ds, cv)

    #remove zeros
    indices = [1:10; 13; 16]
    ke_s = ke_s[indices, indices]


       #static condensation
    inda=1:9
    indi=10:12
    ke_s = ke_s[inda,inda]-ke_s[inda,indi]*inv(ke_s[indi,indi])*ke_s[indi,inda]
    ke_b = ke_b[inda, inda]


    #shear correction
    alpha=sum(diag(ke_s[4:9,4:9]))/sum(diag(ke_b[4:9,4:9]))
    Cs=0.0
    ke_bs=(1/(1+Cs*alpha))*ke_s + ke_b

    # ke_bs = ke_b + ke_s

    #remove zeros
    # indices = [1:10; 13; 16]
    # ke_bs = ke_bs[indices, indices]



    ke = zeros(Float64, 15, 15)

    induv=[1 2 6 7 11 12]
    indwt=[3 4 5 8 9 10 13 14 15]
    ke[induv,induv]=ke_m
    ke[indwt,indwt]=ke_bs

    #add drilling dof

    k15 = ke

    ind=[1:5; 7:11; 13:17]
    k18 = zeros(Float64, 18, 18)
    k18[ind,ind] = k15

    kd = diag(k15)
    ind=[4 5 9 10 14 15];
    stif=minimum(kd[ind])/100;
    k18[6,6] = stif;
    k18[12,12] = stif;
    k18[18,18] = stif;

    # # transform stiffness matrix back to global X-Y-Z

    # P1=nodes(i1,1:3); P2=nodes(i2,1:3); P3=nodes(i3,1:3);
    # norm_vec=cross(P2-P1,P3-P1);
    # norm_vec=norm_vec/norm(norm_vec);
    # j3=norm_vec;
    # j1=(P2-P1)/norm(P2-P1);
    # j2=cross(j3,j1);
    # T=[j1' j2' j3'];


    # #reorder from component to fields, Ferrite default
    # ind_field = [1, 2, 3, 6, 7, 8, 11, 12, 13, 4, 5, 9, 10, 14, 15]
    # ke = ke[ind_field, ind_field]

    ke = k18

    return ke

end


"""
    assemble_global_Ke!(Ke, dh, qr1, qr3, ip3, ip6, E, ν, t)

Assemble the global elastic stiffness matrix for a shell structure.

For each element, the function:
1. Computes the local coordinate system from 3D node positions via [`calculation_rotation_matrix`](@ref)
2. Projects nodes to the local 2D plane via [`global_nodal_coords_to_planar_coords`](@ref)
3. Builds the 18×18 local stiffness via [`local_elastic_stiffness_matrix!`](@ref)
4. Rotates the local matrix to global coordinates using the 18×18 block-diagonal rotation matrix
5. Reorders DOF from component-based to Ferrite's field-based convention and assembles

# Arguments
- `Ke`: Pre-allocated global sparse matrix (e.g., from `Ferrite.allocate_matrix(dh)`)
- `dh`: `DofHandler` with fields `:u` (translations) and `:θ` (rotations)
- `qr1`: 1-point `QuadratureRule{RefTriangle}`
- `qr3`: 3-point `QuadratureRule{RefTriangle}`
- `ip3`: [`IP3`](@ref) interpolation
- `ip6`: [`IP6`](@ref) interpolation
- `E`: Young's modulus
- `ν`: Poisson's ratio
- `t`: Shell thickness

# Returns
Assembled global stiffness matrix `Ke`.
"""
function assemble_global_Ke!(Ke, dh, qr1, qr3, ip3, ip6, E, ν, t)

    assembler = start_assemble(Ke)
    for cell in CellIterator(dh)

        x_global = getcoordinates(cell)
        println("x_global:", x_global)

        T = calculation_rotation_matrix(x_global)
        println("T:", T)

        x_local = global_nodal_coords_to_planar_coords(x_global, T)

        println("x_local:", x_local)

        ke_local = TriShellFiniteElement.local_elastic_stiffness_matrix!(qr1, qr3, ip3, ip6, E, ν, t, x_local)

        #rotate element stiffness matrix back to global coordinates!
        Te = rotation_matrix_for_element_stiffness_drilling(T)
        ke_global = Te * ke_local * Te'

        # #reorder from component to fields, Ferrite default
        ind_field = [1, 2, 3, 7, 8, 9, 13, 14, 15, 4, 5, 6, 10, 11, 12, 16, 17, 18]
        ke_global = ke_global[ind_field, ind_field]

        assemble!(assembler, celldofs(cell), ke_global)
    end
    return Ke
end


####

function generate_Nuvw_derivative(dNdξ_d)

    Nuvw_d = zeros(Float64, 3, 15)
    Nuvw_d[1, 1:2:5] .= dNdξ_d
    Nuvw_d[2, 2:2:6] .= dNdξ_d
    Nuvw_d[3, 7:3:13] .= dNdξ_d

    return Nuvw_d

end


"""
    calculate_element_geometric_stiffness_matrix(cv, σxx, σyy, τxy, T)

Return the 15×15 element geometric stiffness matrix for linear buckling analysis.

Stress resultants are supplied at each quadrature point and transformed from the
global coordinate system to the local element system using the rotation matrix `T`.
Membrane (u, v) and bending/shear (w, θ) contributions are assembled separately
into the 15×15 matrix (no drilling DOF row/column).

# Arguments
- `cv`: `CellValues` initialised with local 2D node coordinates
- `σxx`: Vector of global X-normal stresses at each quadrature point
- `σyy`: Vector of global Y-normal stresses at each quadrature point
- `τxy`: Vector of global XY-shear stresses at each quadrature point
- `T`: 3×3 rotation matrix from [`calculation_rotation_matrix`](@ref)

# Returns
15×15 geometric stiffness matrix in local coordinates.
"""
function calculate_element_geometric_stiffness_matrix(cv, σxx, σyy, τxy, T)

    # reinit!(cv, x)

    # num_shape_functions = getnbasefunctions(ip_shape)

    # kgx = zeros(Float64, 15, 15)
    # kgy = zeros(Float64, 15, 15)
    # kgxy = zeros(Float64, 15, 15)

    kuv = zeros(Float64, 6, 6)
    kwt = zeros(Float64, 9, 9)
    kg = zeros(Float64, 15, 15)

    dNdx = cv.fun_values.dNdx

    println("dNdx:", dNdx)

    for q_point in 1:getnquadpoints(cv)

        # ξ = qr.points[q_point]
        # J = TriShellFiniteElement.get_jacobian(ξ, ip_geo, x)
        # Jinv = inv(J)

        # dNdξ_x = [cv.fun_values.dNdξ[i + (q_point-1)*num_shape_functions][1] for i=1:num_shape_functions]
        # dNdξ_y = [cv.fun_values.dNdξ[i + (q_point-1)*num_shape_functions][2] for i=1:num_shape_functions]

        # Nx = generate_Nuvw_derivative(dNdξ_x)
        # Ny = generate_Nuvw_derivative(dNdξ_y)


        dNx = [dNdx[1, q_point][1]  0.0 dNdx[2, q_point][1] 0.0 dNdx[3, q_point][1] 0.0
                0.0  dNdx[1, q_point][1]    0.0 dNdx[2, q_point][1] 0.0 dNdx[3, q_point][1] ]

        println("dNx_uv:", dNx)

        dNy = [dNdx[1, q_point][2]  0.0 dNdx[2, q_point][2] 0.0 dNdx[3, q_point][2] 0.0
            0.0  dNdx[1, q_point][2]    0.0 dNdx[2, q_point][2] 0.0 dNdx[3, q_point][2] ]

        println("dNy_uv:", dNy)

           GGuvx=dNx'*dNx;
            GGuvy=dNy'*dNy;
            GGuvxy=dNx'*dNy+dNy'*dNx;


              println("GGuvx:", GGuvx)
              println("GGuvy:", GGuvy)
              println("GGuvxy:", GGuvxy)

        # dNx[1, 1] = dNdx[1, q_point][1]
        # dNx[1, 4] = dNdx[2, q_point][1]
        # dNx[1, 7] = dNdx[3, q_point][1]

        # dNy = zeros(Float64, 3, 9)

        # dNy[1, 1] = dNdx[1, q_point][2]
        # dNy[1, 4] = dNdx[2, q_point][2]
        # dNy[1, 7] = dNdx[3, q_point][2]




        dNx = zeros(Float64, 3, 9)

        dNx[1, 1] = dNdx[1, q_point][1]
        dNx[1, 4] = dNdx[2, q_point][1]
        dNx[1, 7] = dNdx[3, q_point][1]

        dNy = zeros(Float64, 3, 9)

        dNy[1, 1] = dNdx[1, q_point][2]
        dNy[1, 4] = dNdx[2, q_point][2]
        dNy[1, 7] = dNdx[3, q_point][2]

        GGwtx=dNx'*dNx;
        GGwty=dNy'*dNy;
        GGwtxy=dNx'*dNy+dNy'*dNx;

                println("GGwtx:", GGwtx)
              println("GGwty:", GGwty)
              println("GGwtxy:", GGwtxy)


        # str_mat = [σxx[q_point]     τxy[q_point]
        #            τxy[q_point]     σyy[q_point]]


        # str_mat = T[1:2, 1:2]' * str_mat * T[1:2, 1:2]

        str_vec=[σxx[q_point], σyy[q_point], τxy[q_point]]

        println("str_vec", str_vec)

        c=T[1,1]
        s=T[2,1]
        T2=[c^2 s^2 2*s*c; s^2 c^2 -2*s*c; -s*c s*c c^2-s^2]
        str_vec=T2*str_vec;

        println("T2", T2)
        println("str_vec_trans", str_vec)

        kuv += (GGuvx*str_vec[1] + GGuvy*str_vec[2] + GGuvxy*str_vec[3] ) * getdetJdV(cv, q_point)
        kwt += (GGwtx*str_vec[1] + GGwty*str_vec[2] + GGwtxy*str_vec[3] ) * getdetJdV(cv, q_point)

    end

    induv = [1 2 6 7 11 12]
    indwt = [3 4 5 8 9 10 13 14 15]
    kg[induv,induv] = kuv
    kg[indwt,indwt] = kwt

    return kg

end


# function geometric_stiffness_matrix!(cv, qr1, ip3, x, σ_element)

#     reinit!(cv, x)
#     kg = calculate_element_geometric_stiffness_matrix(cv, ip3, ip3, qr1, x, σ_element)

#     return kg

# end



"""
    assemble_global_Kg!(Kg, dh, σXX, σYY, τXY)

Assemble the global geometric stiffness matrix for linear buckling analysis.

Stress resultants (from a prior linear static solve) are passed as element-level
vectors. For each element the function transforms stresses to local coordinates,
computes the element geometric stiffness via
[`calculate_element_geometric_stiffness_matrix`](@ref), rotates it to global
coordinates, and assembles it into `Kg`.

# Arguments
- `Kg`: Pre-allocated global sparse matrix
- `dh`: `DofHandler`
- `σXX`: Vector of per-element global X-normal stresses
- `σYY`: Vector of per-element global Y-normal stresses
- `τXY`: Vector of per-element global XY-shear stresses

# Returns
Assembled global geometric stiffness matrix `Kg`.
"""
function assemble_global_Kg!(Kg, dh, σXX, σYY, τXY)

    assembler = start_assemble(Kg)
    i = 1
    for cell in CellIterator(dh)

        x_global = getcoordinates(cell)
        println("x_global:", x_global)

        T = calculation_rotation_matrix(x_global)
        println("T:", T)

        x_local = global_nodal_coords_to_planar_coords(x_global, T)

        println("x_local:", x_local)


        #convert stresses from global coordinate system to local coordinate system

        σXX_element = σXX[i]
        σYY_element = σYY[i]
        τXY_element = τXY[i]

        σxx_element
        σyy_element
        τxy_element


        kg_local = calculate_element_geometric_stiffness_matrix(cv, σxx_element, σyy_element, τxy_element, T)

        #rotate element stiffness matrix back to global coordinates!
        Te = rotation_matrix_for_element_stiffness_drilling(T)
        kg_global = Te * kg_local * Te'

        # #reorder from component to fields, Ferrite default
        ind_field = [1, 2, 3, 7, 8, 9, 13, 14, 15, 4, 5, 6, 10, 11, 12, 16, 17, 18]
        kg_global = kg_global[ind_field, ind_field]

        assemble!(assembler, celldofs(cell), kg_global)
    end
    return Kg
end




# function assemble_global_Kg!(Kg, dh, qr1, ip3, σ_global)

#     cv = CellValues(qr1, ip3, ip3)
#     #need to convert global stress σ to local stress at some point
#     assembler = start_assemble(Kg)
#     i = 1
#     for cell in CellIterator(dh)

#         #3D
#         x_global = getcoordinates(cell)

#         T = calculation_rotation_matrix(x_global)

#         str_mat_global = [σ_global[i][1] σ_global[i][3]
#                           σ_global[i][3] σ_global[i][2]]

#         str_mat_local =  T[1:2,1:2]' * str_mat_global * T[1:2, 1:2]

#         σ_local = [str_mat_local[1, 1], str_mat_local[2, 2], str_mat_local[1, 2]]

#         x_local = global_nodal_coords_to_planar_coords(x_global, T)

#         kg = geometric_stiffness_matrix!(cv, qr1, ip3, x_local, σ_local)

#         #rotate element stiffness matrix back to global coordinates!
#         Te = rotation_matrix_for_element_stiffness(T)
#         kg_global = Te * kg * Te'

#         assemble!(assembler, celldofs(cell), kg_global)
#         i += 1
#     end
#     return Kg
# end


"""
    calculation_rotation_matrix(node)

Compute the 3×3 rotation matrix that defines the local coordinate system of a
triangular shell element from its three 3D node positions.

The local axes are:
- **j₁**: unit vector along edge P1 → P2
- **j₃**: unit normal to the element plane (`(P2−P1) × (P3−P1)`)
- **j₂**: `j₃ × j₁` (completes the right-handed system)

The returned matrix `T = [j₁ j₂ j₃]` transforms from local to global coordinates:
`x_global = T * x_local`.

# Arguments
- `node`: Vector of three `Tensors.Vec{3}` global node coordinates

# Returns
3×3 rotation matrix `T`.
"""
function calculation_rotation_matrix(node)

    P1 = node[1]
    P2 = node[2]
    P3 = node[3]

    norm_vec=cross(P2-P1,P3-P1)
    norm_vec=norm_vec/norm(norm_vec)
    j3=norm_vec
    j1=(P2-P1)/norm(P2-P1)
    j2=cross(j3,j1)
    T=[j1 j2 j3]

    return T

end


"""
    global_nodal_coords_to_planar_coords(cell_nodes_global, T)

Project three 3D global node coordinates into the element's local 2D plane.

The first node is taken as the origin. Each node is transformed by `T'` and
the third component (normal direction) is dropped, returning 2D local coordinates
suitable for use with `CellValues`.

# Arguments
- `cell_nodes_global`: Vector of three `Tensors.Vec{3}` global node coordinates
- `T`: 3×3 rotation matrix from [`calculation_rotation_matrix`](@ref)

# Returns
Vector of three `Tensors.Vec{2}` local planar coordinates.
"""
function global_nodal_coords_to_planar_coords(cell_nodes_global, T)

    P1 = cell_nodes_global[1]
    P2 = cell_nodes_global[2]
    P3 = cell_nodes_global[3]

    P1a=T'*(P1-P1)
    P2a=T'*(P2-P1)
    P3a=T'*(P3-P1)

    cell_local = [P1a, P2a, P3a]

    cell_nodes_local = [Tensors.Vec((cell_local[i][1], cell_local[i][2])) for i in eachindex(cell_local)]

    return cell_nodes_local

end


function rotation_matrix_for_element_stiffness_no_drilling(T3)

    #no drilling dof yet

    T2=T3[1:2,1:2]

    T = Matrix(1.0I, 15, 15)

    ind=[1, 2, 3]
    T[ind,ind] = T3
    ind=[4 5 6]
    T[ind,ind] = T3
    ind=[7 8 9]
    T[ind,ind] = T3
    ind=[10, 11]
    T[ind,ind] = T2
    ind=[12, 13]
    T[ind,ind] = T2
    ind=[14, 15]
    T[ind,ind] = T2

    return T

end


"""
    rotation_matrix_for_element_stiffness_drilling(T3)

Build the 18×18 block-diagonal rotation matrix for transforming an 18-DOF shell
element stiffness matrix between local and global coordinates.

The 3×3 rotation matrix `T3` is placed on the diagonal in six 3×3 blocks —
one block per set of translational DOF `[u, v, w]` and one per set of rotational
DOF `[θₓ, θᵧ, θᵤ]` at each of the three nodes.

# Arguments
- `T3`: 3×3 element rotation matrix from [`calculation_rotation_matrix`](@ref)

# Returns
18×18 block-diagonal rotation matrix `Tₑ` such that
`Kₑ_global = Tₑ * Kₑ_local * Tₑᵀ`.
"""
function rotation_matrix_for_element_stiffness_drilling(T3)

    #no drilling dof yet

    # T2=T3[1:2,1:2]

    T = Matrix(1.0I, 18, 18)

    ind=[1, 2, 3]
    T[ind,ind] = T3
    ind=[7, 8, 9]
    T[ind,ind] = T3
    ind=[13, 14, 15]
    T[ind,ind] = T3
    ind=[4, 5, 6]
    T[ind,ind] = T3
    ind=[10, 11, 12]
    T[ind,ind] = T3
    ind=[16, 17, 18]
    T[ind,ind] = T3

    return T

end


end # module TriShellFiniteElement
