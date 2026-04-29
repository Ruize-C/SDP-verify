module VerifySDP

# Core numerical routines for the SDP verification experiments.
#
# Organization:
#   1. model and moment-set definitions;
#   2. local one- and two-dimer operators;
#   3. fermionic basis utilities and Hamiltonian construction;
#   4. exact-diagonalization observables;
#   5. reduced two-site-marginal SDP construction;
#   6. finite-difference response and experiment helpers.

using LinearAlgebra
using SparseArrays
using JuMP
using SCS
using MosekTools

export MomentSpec,
    DimerModel,
    SpinfulHubbardDimerModel,
    dimer_model,
    spinful_hubbard_dimer_model,
    physical_observables,
    solve_sdp_observables,
    moment_ranges_on_optimal_face,
    run_probe,
    run_spinful_hubbard_probe,
    hopping_op,
    density_density_op,
    imbalance_op

const Mat = Matrix{Float64}

"""
    MomentSpec(name, sites, op)

One selected two-site moment. `sites` gives the physical lattice sites and
`op` is the local matrix in the corresponding two-site occupation basis.
"""
struct MomentSpec
    name::String
    sites::Tuple{Int, Int}
    op::Mat
end

"""
Spinless dimer-chain model used for small smoke tests of the SDP machinery.
"""
struct DimerModel
    L::Int
    Np::Int
    td::Float64
    perturbation::Symbol
    selected::Vector{MomentSpec}
end

"""
Spinful Hubbard dimer-chain model used in the reported experiments.

The reference point is a product of isolated strong dimers. `perturbation`
selects the weak inter-dimer term used in `H(theta)`.
"""
struct SpinfulHubbardDimerModel
    L::Int
    Np::Int
    td::Float64
    U::Float64
    stagger::Float64
    perturbation::Symbol
    selected::Vector{MomentSpec}
end

"""
    dimer_model(; ...)

Build a small spinless dimer-chain test instance.
"""
function dimer_model(; L::Int = 4, Np::Int = div(L, 2), td::Real = 1.0,
    perturbation::Symbol = :density, selected::Symbol = :dimer_hopping)
    iseven(L) || error("L must be even for the dimer reference model.")
    0 <= Np <= L || error("Np must satisfy 0 <= Np <= L.")

    moments = MomentSpec[]
    if selected == :dimer_hopping
        for i in 1:2:L
            push!(moments, MomentSpec("hop_$(i)_$(i + 1)", (i, i + 1), hopping_op(1.0)))
        end
    elseif selected == :dimer_imbalance
        for i in 1:2:L
            push!(moments, MomentSpec("imbalance_$(i)_$(i + 1)", (i, i + 1), imbalance_op(1.0)))
        end
    elseif selected == :dimer_imbalance_and_hopping
        for i in 1:2:L
            push!(moments, MomentSpec("imbalance_$(i)_$(i + 1)", (i, i + 1), imbalance_op(1.0)))
            push!(moments, MomentSpec("hop_$(i)_$(i + 1)", (i, i + 1), hopping_op(1.0)))
        end
    elseif selected == :inter_density
        for i in 2:2:(L - 2)
            push!(moments, MomentSpec("dens_$(i)_$(i + 1)", (i, i + 1), density_density_op(1.0)))
        end
    elseif selected == :dimer_imbalance_and_inter_density
        for i in 1:2:L
            push!(moments, MomentSpec("imbalance_$(i)_$(i + 1)", (i, i + 1), imbalance_op(1.0)))
        end
        for i in 2:2:(L - 2)
            push!(moments, MomentSpec("dens_$(i)_$(i + 1)", (i, i + 1), density_density_op(1.0)))
        end
    else
        error("unknown selected moment set: $selected")
    end
    return DimerModel(L, Np, Float64(td), perturbation, moments)
end

"""
    spinful_hubbard_dimer_model(; ...)

Build a spinful Hubbard dimer-chain instance. `selected` chooses the moment
family used in the SDP response and mismatch checks.
"""
function spinful_hubbard_dimer_model(; L::Int = 4, Np::Int = L, td::Real = 1.0,
    U::Real = 4.0, stagger::Real = 0.0,
    perturbation::Symbol = :hopping, selected::Symbol = :charge_imbalance)
    iseven(L) || error("L must be even for the dimer reference model.")
    0 <= Np <= 2L || error("Np must satisfy 0 <= Np <= 2L.")

    moments = MomentSpec[]
    if selected == :charge_imbalance
        for i in 1:2:L
            push!(moments, MomentSpec("charge_imbalance_$(i)_$(i + 1)",
                (i, i + 1), spinful_charge_imbalance_op(1.0)))
        end
    elseif selected == :bond_kinetic
        for i in 1:2:L
            push!(moments, MomentSpec("bond_kinetic_$(i)_$(i + 1)",
                (i, i + 1), spinful_pair_hopping_op(1.0)))
        end
    elseif selected == :charge_and_bond
        for i in 1:2:L
            push!(moments, MomentSpec("charge_imbalance_$(i)_$(i + 1)",
                (i, i + 1), spinful_charge_imbalance_op(1.0)))
            push!(moments, MomentSpec("bond_kinetic_$(i)_$(i + 1)",
                (i, i + 1), spinful_pair_hopping_op(1.0)))
        end
    elseif selected == :double_occupancy
        for i in 1:2:L
            push!(moments, MomentSpec("double_occupancy_$(i)_$(i + 1)",
                (i, i + 1), spinful_double_occupancy_sum_op(1.0)))
        end
    elseif selected == :double_occupancy_sites
        for i in 1:2:L
            push!(moments, MomentSpec("double_occupancy_left_$(i)",
                (i, i + 1), spinful_double_occupancy_left_op(1.0)))
            push!(moments, MomentSpec("double_occupancy_right_$(i + 1)",
                (i, i + 1), spinful_double_occupancy_right_op(1.0)))
        end
    elseif selected == :double_occupancy_imbalance
        for i in 1:2:L
            push!(moments, MomentSpec("double_occupancy_imbalance_$(i)_$(i + 1)",
                (i, i + 1), spinful_double_occupancy_imbalance_op(1.0)))
        end
    elseif selected == :inter_bond_kinetic
        for i in 2:2:(L - 1)
            push!(moments, MomentSpec("inter_bond_kinetic_$(i)_$(i + 1)",
                (i, i + 1), spinful_pair_hopping_op(1.0)))
        end
    elseif selected == :inter_density
        for i in 2:2:(L - 1)
            push!(moments, MomentSpec("inter_density_$(i)_$(i + 1)",
                (i, i + 1), spinful_inter_density_op(1.0)))
        end
    elseif selected == :inter_charge_imbalance
        for i in 2:2:(L - 1)
            push!(moments, MomentSpec("inter_charge_imbalance_$(i)_$(i + 1)",
                (i, i + 1), spinful_charge_imbalance_op(1.0)))
        end
    elseif selected == :dimer_spin_z_corr
        for i in 1:2:L
            push!(moments, MomentSpec("spin_z_corr_$(i)_$(i + 1)",
                (i, i + 1), spinful_spin_z_corr_op(1.0)))
        end
    elseif selected == :dimer_spin_exchange
        for i in 1:2:L
            push!(moments, MomentSpec("spin_exchange_$(i)_$(i + 1)",
                (i, i + 1), spinful_spin_exchange_op(1.0)))
        end
    elseif selected == :dimer_spin_dot
        for i in 1:2:L
            push!(moments, MomentSpec("spin_dot_$(i)_$(i + 1)",
                (i, i + 1), spinful_spin_dot_op(1.0)))
        end
    elseif selected == :dimer_pair_hopping
        for i in 1:2:L
            push!(moments, MomentSpec("pair_hopping_$(i)_$(i + 1)",
                (i, i + 1), spinful_pair_hopping_pair_op(1.0)))
        end
    elseif selected == :dimer_singlet_projector
        for i in 1:2:L
            push!(moments, MomentSpec("singlet_projector_$(i)_$(i + 1)",
                (i, i + 1), spinful_singlet_projector_op(1.0)))
        end
    elseif selected == :dimer_triplet0_projector
        for i in 1:2:L
            push!(moments, MomentSpec("triplet0_projector_$(i)_$(i + 1)",
                (i, i + 1), spinful_triplet0_projector_op(1.0)))
        end
    else
        error("unknown selected moment set: $selected")
    end
    return SpinfulHubbardDimerModel(L, Np, Float64(td), Float64(U), Float64(stagger), perturbation, moments)
end

# ---------------------------------------------------------------------------
# Local operators
# ---------------------------------------------------------------------------

function hopping_op(coeff::Real)
    H = zeros(4, 4)
    # Local two-site basis: |00>, |01>, |10>, |11>.
    H[2, 3] = coeff
    H[3, 2] = coeff
    return H
end

function density_density_op(coeff::Real)
    H = zeros(4, 4)
    H[4, 4] = coeff
    return H
end

function imbalance_op(coeff::Real)
    H = zeros(4, 4)
    # n_left - n_right in the basis |00>, |01>, |10>, |11>.
    H[2, 2] = -coeff
    H[3, 3] = coeff
    return H
end

function number_op()
    return [0.0 0.0; 0.0 1.0]
end

function spinful_number_op()
    # Single-site basis: |0>, |up>, |down>, |up down>.
    return Diagonal([0.0, 1.0, 1.0, 2.0]) |> Matrix
end

function fock_to_ordered_index(state::UInt, k::Int)
    idx = 0
    for pos in 1:k
        bit = Int((state >> (pos - 1)) & UInt(1))
        idx |= bit << (k - pos)
    end
    return idx + 1
end

function ordered_index_to_fock(idx::Int, k::Int)
    idx0 = idx - 1
    state = UInt(0)
    for pos in 1:k
        bit = (idx0 >> (k - pos)) & 1
        if bit == 1
            state |= UInt(1) << (pos - 1)
        end
    end
    return state
end

function add_local_hopping!(A::Mat, p::Int, q::Int, coeff::Float64)
    k = round(Int, log2(size(A, 1)))
    for col in axes(A, 2)
        state = ordered_index_to_fock(col, k)
        r = apply_cdag_c(state, p, q)
        r === nothing && continue
        sign, new_state = r
        row = fock_to_ordered_index(new_state, k)
        A[row, col] += coeff * sign
    end
    return A
end

function add_local_density!(A::Mat, p::Int, coeff::Float64)
    k = round(Int, log2(size(A, 1)))
    for col in axes(A, 2)
        state = ordered_index_to_fock(col, k)
        if occupied(state, p)
            A[col, col] += coeff
        end
    end
    return A
end

function add_local_density_density!(A::Mat, p::Int, q::Int, coeff::Float64)
    k = round(Int, log2(size(A, 1)))
    for col in axes(A, 2)
        state = ordered_index_to_fock(col, k)
        if occupied(state, p) && occupied(state, q)
            A[col, col] += coeff
        end
    end
    return A
end

function apply_local_sequence(state::UInt, sequence::Vector{Tuple{Symbol, Int}})
    sign = 1.0
    out = state
    for (kind, mode) in sequence
        r = kind == :annihilate ? apply_annihilate(out, mode) :
            kind == :create ? apply_create(out, mode) :
            error("unknown local fermion operation: $kind")
        r === nothing && return nothing
        s, out = r
        sign *= s
    end
    return sign, out
end

function add_local_sequence!(A::Mat, sequence::Vector{Tuple{Symbol, Int}}, coeff::Float64)
    k = round(Int, log2(size(A, 1)))
    for col in axes(A, 2)
        state = ordered_index_to_fock(col, k)
        r = apply_local_sequence(state, sequence)
        r === nothing && continue
        sign, new_state = r
        row = fock_to_ordered_index(new_state, k)
        A[row, col] += coeff * sign
    end
    return A
end

function spinful_pair_hopping_op(coeff::Real)
    A = zeros(16, 16)
    c = Float64(coeff)
    # Local two-site mode order: left up, left down, right up, right down.
    add_local_hopping!(A, 1, 3, c)
    add_local_hopping!(A, 3, 1, c)
    add_local_hopping!(A, 2, 4, c)
    add_local_hopping!(A, 4, 2, c)
    return A
end

function spinful_inter_density_op(coeff::Real)
    A = zeros(16, 16)
    c = Float64(coeff)
    for p in (1, 2), q in (3, 4)
        add_local_density_density!(A, p, q, c)
    end
    return A
end

function spinful_charge_imbalance_op(coeff::Real)
    A = zeros(16, 16)
    c = Float64(coeff)
    add_local_density!(A, 1, c)
    add_local_density!(A, 2, c)
    add_local_density!(A, 3, -c)
    add_local_density!(A, 4, -c)
    return A
end

function spinful_double_occupancy_sum_op(coeff::Real)
    A = zeros(16, 16)
    c = Float64(coeff)
    add_local_density_density!(A, 1, 2, c)
    add_local_density_density!(A, 3, 4, c)
    return A
end

function spinful_double_occupancy_left_op(coeff::Real)
    A = zeros(16, 16)
    add_local_density_density!(A, 1, 2, Float64(coeff))
    return A
end

function spinful_double_occupancy_right_op(coeff::Real)
    A = zeros(16, 16)
    add_local_density_density!(A, 3, 4, Float64(coeff))
    return A
end

function spinful_double_occupancy_imbalance_op(coeff::Real)
    return spinful_double_occupancy_left_op(coeff) .- spinful_double_occupancy_right_op(coeff)
end

function spinful_spin_z_corr_op(coeff::Real)
    A = zeros(16, 16)
    c = Float64(coeff) / 4.0
    add_local_density_density!(A, 1, 3, c)
    add_local_density_density!(A, 1, 4, -c)
    add_local_density_density!(A, 2, 3, -c)
    add_local_density_density!(A, 2, 4, c)
    return A
end

function spinful_spin_exchange_op(coeff::Real)
    A = zeros(16, 16)
    c = Float64(coeff)
    # S_L^+ S_R^- + S_L^- S_R^+
    add_local_sequence!(A, [(:annihilate, 3), (:create, 4), (:annihilate, 2), (:create, 1)], c)
    add_local_sequence!(A, [(:annihilate, 4), (:create, 3), (:annihilate, 1), (:create, 2)], c)
    return A
end

function spinful_spin_dot_op(coeff::Real)
    return spinful_spin_z_corr_op(coeff) .+ 0.5 .* spinful_spin_exchange_op(coeff)
end

function spinful_pair_hopping_pair_op(coeff::Real)
    A = zeros(16, 16)
    c = Float64(coeff)
    # Moves an onsite pair between the two sites of a dimer.
    add_local_sequence!(A, [(:annihilate, 3), (:annihilate, 4), (:create, 2), (:create, 1)], c)
    add_local_sequence!(A, [(:annihilate, 1), (:annihilate, 2), (:create, 4), (:create, 3)], c)
    return A
end

function local_projector_op(indices::Vector{Int}, coeff::Real)
    A = zeros(16, 16)
    v = zeros(16)
    for (idx, val) in zip(indices[1:2:end], indices[2:2:end])
        v[idx] = val
    end
    return Float64(coeff) .* (v * v')
end

function spinful_singlet_projector_op(coeff::Real)
    # Local ordered basis modes are (left up, left down, right up, right down).
    # |up,down> has local state 1001, index 10; |down,up> has state 0110, index 7.
    return local_projector_op([10, 1, 7, -1] , Float64(coeff) / 2.0)
end

function spinful_triplet0_projector_op(coeff::Real)
    return local_projector_op([10, 1, 7, 1], Float64(coeff) / 2.0)
end

function spinful_hubbard_dimer_op(td::Real, U::Real, stagger::Real = 0.0)
    return -Float64(td) .* spinful_pair_hopping_op(1.0) .+
        Float64(U) .* spinful_double_occupancy_sum_op(1.0) .+
        Float64(stagger) .* spinful_charge_imbalance_op(1.0)
end

# ---------------------------------------------------------------------------
# Fermionic basis utilities
# ---------------------------------------------------------------------------

function basis_states(L::Int, Np::Int)
    return [UInt(s) for s in 0:(UInt(1) << L)-UInt(1) if count_ones(s) == Np]
end

@inline function occupied(state::UInt, i::Int)
    return ((state >> (i - 1)) & UInt(1)) == UInt(1)
end

@inline function fermion_sign(state::UInt, i::Int)
    mask = (UInt(1) << (i - 1)) - UInt(1)
    return isodd(count_ones(state & mask)) ? -1.0 : 1.0
end

function apply_annihilate(state::UInt, i::Int)
    occupied(state, i) || return nothing
    return (fermion_sign(state, i), state & ~(UInt(1) << (i - 1)))
end

function apply_create(state::UInt, i::Int)
    occupied(state, i) && return nothing
    return (fermion_sign(state, i), state | (UInt(1) << (i - 1)))
end

function apply_cdag_c(state::UInt, i::Int, j::Int)
    r1 = apply_annihilate(state, j)
    r1 === nothing && return nothing
    s1, state1 = r1
    r2 = apply_create(state1, i)
    r2 === nothing && return nothing
    s2, state2 = r2
    return (s1 * s2, state2)
end

# Add coeff * (c_i^dag c_j + h.c.) to a fixed-particle-number Hamiltonian.
function add_hopping!(H::Mat, states::Vector{UInt}, index::Dict{UInt, Int}, i::Int, j::Int, coeff::Float64)
    for (col, state) in pairs(states)
        for (a, b) in ((i, j), (j, i))
            r = apply_cdag_c(state, a, b)
            r === nothing && continue
            sign, new_state = r
            row = index[new_state]
            H[row, col] += coeff * sign
        end
    end
    return H
end

# Add a diagonal density-density term coeff * n_i n_j.
function add_density_density!(H::Mat, states::Vector{UInt}, i::Int, j::Int, coeff::Float64)
    for (col, state) in pairs(states)
        if occupied(state, i) && occupied(state, j)
            H[col, col] += coeff
        end
    end
    return H
end

function dimer_pairs(L::Int)
    return [(i, i + 1) for i in 1:2:L]
end

function inter_pairs(L::Int)
    return [(i, i + 1) for i in 2:2:(L - 1)]
end

mode_index(site::Int, spin::Int) = 2 * (site - 1) + spin

# A spinful site is represented by two fermionic modes: up and down.
# For a two-site operator, this returns the four mode indices in local order.
function spinful_pair_modes(sites::Tuple{Int, Int})
    i, j = sites
    return [mode_index(i, 1), mode_index(i, 2), mode_index(j, 1), mode_index(j, 2)]
end

# ---------------------------------------------------------------------------
# Hamiltonian construction
# ---------------------------------------------------------------------------

function hamiltonian_matrix(model::DimerModel; theta::Real = 0.0, fields::AbstractVector{<:Real} = Float64[])
    states = basis_states(model.L, model.Np)
    index = Dict(s => k for (k, s) in pairs(states))
    H = zeros(length(states), length(states))

    for (i, j) in dimer_pairs(model.L)
        add_hopping!(H, states, index, i, j, -model.td)
    end

    if theta != 0
        for (i, j) in inter_pairs(model.L)
            if model.perturbation == :density
                add_density_density!(H, states, i, j, Float64(theta))
            elseif model.perturbation == :hopping
                add_hopping!(H, states, index, i, j, -Float64(theta))
            else
                error("unknown perturbation: $(model.perturbation)")
            end
        end
    end

    if !isempty(fields)
        length(fields) == length(model.selected) || error("field length must match selected moments.")
        for (h, moment) in zip(fields, model.selected)
            add_pair_operator!(H, states, index, moment.sites, -Float64(h) .* moment.op)
        end
    end
    return H, states
end

function hamiltonian_matrix(model::SpinfulHubbardDimerModel; theta::Real = 0.0,
    fields::AbstractVector{<:Real} = Float64[])
    nmode = 2 * model.L
    states = basis_states(nmode, model.Np)
    index = Dict(s => k for (k, s) in pairs(states))
    H = zeros(length(states), length(states))

    for (i, j) in dimer_pairs(model.L)
        for spin in 1:2
            add_hopping!(H, states, index, mode_index(i, spin), mode_index(j, spin), -model.td)
        end
        for site in (i, j)
            add_density_density!(H, states, mode_index(site, 1), mode_index(site, 2), model.U)
        end
        if model.stagger != 0
            for spin in 1:2, (site, coeff) in ((i, model.stagger), (j, -model.stagger))
                mode = mode_index(site, spin)
                add_density_density!(H, states, mode, mode, coeff)
            end
        end
    end

    if theta != 0
        for (i, j) in inter_pairs(model.L)
            if model.perturbation == :density
                for spin_i in 1:2, spin_j in 1:2
                    add_density_density!(H, states, mode_index(i, spin_i), mode_index(j, spin_j), Float64(theta))
                end
            elseif model.perturbation == :hopping
                for spin in 1:2
                    add_hopping!(H, states, index, mode_index(i, spin), mode_index(j, spin), -Float64(theta))
                end
            else
                error("unknown perturbation: $(model.perturbation)")
            end
        end
    end

    if !isempty(fields)
        length(fields) == length(model.selected) || error("field length must match selected moments.")
        for (h, moment) in zip(fields, model.selected)
            add_spinful_pair_operator!(H, states, index, moment.sites, -Float64(h) .* moment.op)
        end
    end
    return H, states
end

# Add a two-site spinless local operator to a many-body Hamiltonian.
function add_pair_operator!(H::Mat, states::Vector{UInt}, index::Dict{UInt, Int},
    sites::Tuple{Int, Int}, op::Mat)
    i, j = sites
    for (col, state) in pairs(states)
        local_col = local_index(state, sites)
        for local_row in 1:4
            amp = op[local_row, local_col]
            abs(amp) <= 1e-14 && continue
            new_state = replace_local_bits(state, sites, local_row)
            count_ones(new_state) == count_ones(state) || continue
            row = get(index, new_state, 0)
            row == 0 && continue
            H[row, col] += amp
        end
    end
    return H
end

# Add a two-site spinful local operator. This is a four-mode operator because
# each physical site has up and down fermionic modes.
function add_spinful_pair_operator!(H::Mat, states::Vector{UInt}, index::Dict{UInt, Int},
    sites::Tuple{Int, Int}, op::Mat)
    modes = spinful_pair_modes(sites)
    for (col, state) in pairs(states)
        local_col = subsystem_index(Int(state), modes) + 1
        for local_row in axes(op, 1)
            amp = op[local_row, local_col]
            abs(amp) <= 1e-14 && continue
            new_state = replace_modes(state, modes, local_row)
            count_ones(new_state) == count_ones(state) || continue
            row = get(index, new_state, 0)
            row == 0 && continue
            H[row, col] += amp
        end
    end
    return H
end

function replace_modes(state::UInt, modes::Vector{Int}, local_row::Int)
    out = state
    idx0 = local_row - 1
    k = length(modes)
    for (pos, mode) in pairs(modes)
        bit = (idx0 >> (k - pos)) & 1
        out = bit == 1 ? (out | (UInt(1) << (mode - 1))) : (out & ~(UInt(1) << (mode - 1)))
    end
    return out
end

function local_index(state::UInt, sites::Tuple{Int, Int})
    i, j = sites
    a = occupied(state, i) ? 1 : 0
    b = occupied(state, j) ? 1 : 0
    return 2a + b + 1
end

function replace_local_bits(state::UInt, sites::Tuple{Int, Int}, local_row::Int)
    i, j = sites
    a = (local_row - 1) ÷ 2
    b = (local_row - 1) % 2
    out = state
    out = a == 1 ? (out | (UInt(1) << (i - 1))) : (out & ~(UInt(1) << (i - 1)))
    out = b == 1 ? (out | (UInt(1) << (j - 1))) : (out & ~(UInt(1) << (j - 1)))
    return out
end

# ---------------------------------------------------------------------------
# Exact diagonalization and reduced density matrices
# ---------------------------------------------------------------------------

function ground_state(model::DimerModel; theta::Real = 0.0, fields::AbstractVector{<:Real} = Float64[])
    H, states = hamiltonian_matrix(model; theta, fields)
    F = eigen(Symmetric(H))
    return F.values[1], F.vectors[:, 1], states, F.values
end

function ground_state(model::SpinfulHubbardDimerModel; theta::Real = 0.0,
    fields::AbstractVector{<:Real} = Float64[])
    H, states = hamiltonian_matrix(model; theta, fields)
    F = eigen(Symmetric(H))
    return F.values[1], F.vectors[:, 1], states, F.values
end

function full_wavefunction(vec::Vector{Float64}, states::Vector{UInt}, L::Int)
    psi = zeros(Float64, 1 << L)
    for (amp, state) in zip(vec, states)
        psi[Int(state) + 1] = amp
    end
    return psi
end

function rdm(psi::Vector{Float64}, L::Int, sites::Vector{Int})
    k = length(sites)
    dim = 1 << k
    rho = zeros(dim, dim)
    site_set = Set(sites)
    for s in 0:(1 << L)-1
        outside = outside_bits(s, L, site_set)
        a = subsystem_index(s, sites)
        amps = psi[s + 1]
        abs(amps) <= 1e-14 && continue
        for t in 0:(1 << L)-1
            outside_bits(t, L, site_set) == outside || continue
            b = subsystem_index(t, sites)
            rho[a + 1, b + 1] += amps * psi[t + 1]
        end
    end
    return rho
end

function outside_bits(state::Int, L::Int, sites::Set{Int})
    out = 0
    shift = 0
    for i in 1:L
        i in sites && continue
        bit = (state >> (i - 1)) & 1
        out |= bit << shift
        shift += 1
    end
    return out
end

function subsystem_index(state::Int, sites::Vector{Int})
    out = 0
    for (pos, site) in pairs(sites)
        bit = (state >> (site - 1)) & 1
        out |= bit << (length(sites) - pos)
    end
    return out
end

"""
    physical_observables(model; theta=0.0, fields=[])

Compute exact ground-state energy, selected moments, and spectral gap by exact
diagonalization.
"""
function physical_observables(model::DimerModel; theta::Real = 0.0, fields::AbstractVector{<:Real} = Float64[])
    E0, gs, states, evals = ground_state(model; theta, fields)
    psi = full_wavefunction(gs, states, model.L)
    moments = Float64[]
    for m in model.selected
        rho = rdm(psi, model.L, collect(m.sites))
        push!(moments, tr(m.op * rho))
    end
    gap = length(evals) >= 2 ? evals[2] - evals[1] : NaN
    return (; energy = E0, moments, gap)
end

function physical_observables(model::SpinfulHubbardDimerModel; theta::Real = 0.0,
    fields::AbstractVector{<:Real} = Float64[])
    E0, gs, states, evals = ground_state(model; theta, fields)
    psi = full_wavefunction(gs, states, 2 * model.L)
    moments = Float64[]
    for m in model.selected
        rho = rdm(psi, 2 * model.L, spinful_pair_modes(m.sites))
        push!(moments, tr(m.op * rho))
    end
    gap = length(evals) >= 2 ? evals[2] - evals[1] : NaN
    return (; energy = E0, moments, gap)
end

# ---------------------------------------------------------------------------
# Reduced two-site-marginal SDP
# ---------------------------------------------------------------------------

function all_pairs(L::Int)
    return [(i, j) for i in 1:(L - 1) for j in (i + 1):L]
end

function matrix_unit_basis(C::Int)
    basis = Mat[]
    for a in 1:C, b in 1:C
        E = zeros(C, C)
        E[a, b] = 1.0
        push!(basis, E)
    end
    return basis
end

# Linear JuMP expression for tr(A * X), where X is a PSD matrix variable.
function trace_expr(A::AbstractMatrix{<:Real}, X)
    expr = zero(A[1, 1]) * X[1, 1]
    for i in axes(A, 1), j in axes(A, 2)
        abs(A[i, j]) <= 1e-14 && continue
        expr += A[i, j] * X[j, i]
    end
    return expr
end

# Pair-local objective matrices for the spinless model.
function pair_objectives(model::DimerModel; theta::Real = 0.0, fields::AbstractVector{<:Real} = Float64[])
    objs = Dict{Tuple{Int, Int}, Mat}()
    for p in all_pairs(model.L)
        objs[p] = zeros(4, 4)
    end
    for (i, j) in dimer_pairs(model.L)
        objs[(i, j)] .+= -model.td .* hopping_op(1.0)
    end
    if theta != 0
        for (i, j) in inter_pairs(model.L)
            if model.perturbation == :density
                objs[(i, j)] .+= Float64(theta) .* density_density_op(1.0)
            elseif model.perturbation == :hopping
                objs[(i, j)] .+= -Float64(theta) .* hopping_op(1.0)
            else
                error("unknown perturbation: $(model.perturbation)")
            end
        end
    end
    if !isempty(fields)
        for (h, moment) in zip(fields, model.selected)
            objs[moment.sites] .+= -Float64(h) .* moment.op
        end
    end
    return objs
end

# Pair-local objective matrices for the spinful Hubbard model.
function pair_objectives(model::SpinfulHubbardDimerModel; theta::Real = 0.0,
    fields::AbstractVector{<:Real} = Float64[])
    objs = Dict{Tuple{Int, Int}, Mat}()
    for p in all_pairs(model.L)
        objs[p] = zeros(16, 16)
    end
    for p in dimer_pairs(model.L)
        objs[p] .+= spinful_hubbard_dimer_op(model.td, model.U, model.stagger)
    end
    if theta != 0
        for p in inter_pairs(model.L)
            if model.perturbation == :density
                objs[p] .+= Float64(theta) .* spinful_inter_density_op(1.0)
            elseif model.perturbation == :hopping
                objs[p] .+= -Float64(theta) .* spinful_pair_hopping_op(1.0)
            else
                error("unknown perturbation: $(model.perturbation)")
            end
        end
    end
    if !isempty(fields)
        for (h, moment) in zip(fields, model.selected)
            objs[moment.sites] .+= -Float64(h) .* moment.op
        end
    end
    return objs
end

local_dim(::DimerModel) = 2
local_dim(::SpinfulHubbardDimerModel) = 4

site_number_op(::DimerModel) = number_op()
site_number_op(::SpinfulHubbardDimerModel) = spinful_number_op()

# Construct a JuMP model with either SCS or MOSEK. MOSEK is used for the
# reported high-accuracy runs; SCS is useful for license-free smoke tests.
function make_jump_model(; solver::Symbol = :scs, verbose::Bool = false,
    eps_abs::Real = 1e-7, eps_rel::Real = 1e-7, max_iters::Integer = 50_000)
    if solver == :scs
        jump_model = Model(SCS.Optimizer)
        set_optimizer_attribute(jump_model, "eps_abs", Float64(eps_abs))
        set_optimizer_attribute(jump_model, "eps_rel", Float64(eps_rel))
        set_optimizer_attribute(jump_model, "max_iters", max_iters)
    elseif solver == :mosek
        jump_model = Model(MosekTools.Optimizer)
        set_optimizer_attribute(jump_model, "MSK_DPAR_INTPNT_CO_TOL_PFEAS", Float64(eps_abs))
        set_optimizer_attribute(jump_model, "MSK_DPAR_INTPNT_CO_TOL_DFEAS", Float64(eps_abs))
        set_optimizer_attribute(jump_model, "MSK_DPAR_INTPNT_CO_TOL_REL_GAP", Float64(eps_rel))
        set_optimizer_attribute(jump_model, "MSK_IPAR_INTPNT_MAX_ITERATIONS", max_iters)
    else
        error("unknown solver: $solver")
    end

    if verbose
        unset_silent(jump_model)
    else
        set_silent(jump_model)
    end
    return jump_model
end

# Build the reduced SDP over all two-site marginals.
#
# The constraints enforce PSD marginals, trace normalization, fixed particle
# number, and consistency between overlapping one-site reductions.
function build_sdp_model(model; theta::Real = 0.0,
    fields::AbstractVector{<:Real} = Float64[], verbose::Bool = false,
    eps_abs::Real = 1e-7, eps_rel::Real = 1e-7, max_iters::Integer = 50_000,
    solver::Symbol = :scs)
    C = local_dim(model)
    basis = matrix_unit_basis(C)
    d = length(basis)
    L = model.L
    allp = all_pairs(L)

    jump_model = make_jump_model(; solver, verbose, eps_abs, eps_rel, max_iters)

    rhoi = Any[]
    for i in 1:L
        R = @variable(jump_model, [1:C, 1:C], PSD)
        push!(rhoi, R)
    end

    rhoij = Dict{Tuple{Int, Int}, Any}()
    for p in allp
        R2 = @variable(jump_model, [1:(C * C), 1:(C * C)], PSD)
        rhoij[p] = R2
    end

    M = @variable(jump_model, [1:(L * d + 1), 1:(L * d + 1)], PSD)

    for i in 1:L
        @constraint(jump_model, tr_expr_identity(rhoi[i]) == 1.0)
    end
    @constraint(jump_model, sum(trace_expr(site_number_op(model), rhoi[i]) for i in 1:L) == model.Np)

    for p in allp
        i, j = p
        R2 = rhoij[p]
        for a in 1:C, ap in 1:C
            @constraint(jump_model,
                sum(R2[two_index(a, b, C), two_index(ap, b, C)] for b in 1:C) == rhoi[i][a, ap])
        end
        for b in 1:C, bp in 1:C
            @constraint(jump_model,
                sum(R2[two_index(a, b, C), two_index(a, bp, C)] for a in 1:C) == rhoi[j][b, bp])
        end
    end

    last = L * d + 1
    @constraint(jump_model, M[last, last] == 1.0)
    for i in 1:L
        for α in 1:d
            row = block_index(i, α, d)
            @constraint(jump_model, M[row, last] == trace_expr(transpose(basis[α]), rhoi[i]))
            @constraint(jump_model, M[last, row] == M[row, last])
        end
        for α in 1:d, β in 1:d
            row = block_index(i, α, d)
            col = block_index(i, β, d)
            A = transpose(basis[α]) * basis[β]
            @constraint(jump_model, M[row, col] == trace_expr(A, rhoi[i]))
        end
    end

    for p in allp
        i, j = p
        R2 = rhoij[p]
        for α in 1:d, β in 1:d
            row = block_index(i, α, d)
            col = block_index(j, β, d)
            A = kron(transpose(basis[α]), basis[β])
            @constraint(jump_model, M[row, col] == trace_expr(A, R2))
            @constraint(jump_model, M[col, row] == M[row, col])
        end
    end

    pair_objs = pair_objectives(model; theta, fields)
    energy_expr = sum(trace_expr(pair_objs[p], rhoij[p]) for p in allp)

    return (; jump_model, rhoij, energy_expr)
end

"""
    solve_sdp_observables(model; theta=0.0, fields=[], solver=:scs, ...)

Solve the reduced SDP and return the objective value, selected SDP moments, and
solver status.
"""
function solve_sdp_observables(model; theta::Real = 0.0,
    fields::AbstractVector{<:Real} = Float64[], verbose::Bool = false,
    eps_abs::Real = 1e-7, eps_rel::Real = 1e-7, max_iters::Integer = 50_000,
    solver::Symbol = :scs)
    problem = build_sdp_model(model; theta, fields, verbose, eps_abs, eps_rel, max_iters, solver)
    jump_model = problem.jump_model
    rhoij = problem.rhoij

    @objective(jump_model, Min, problem.energy_expr)
    optimize!(jump_model)

    term = termination_status(jump_model)
    primal = primal_status(jump_model)
    optval = objective_value(jump_model)
    moments = Float64[]
    for spec in model.selected
        push!(moments, value(trace_expr(spec.op, rhoij[spec.sites])))
    end

    return (; energy = optval, moments, termination = string(term), primal_status = string(primal))
end

# Optimize one selected moment over a near-optimal SDP face. This is a
# diagnostic helper for checking moment non-uniqueness at fixed SDP energy.
function optimize_moment_on_face(model, spec::MomentSpec, sense::Symbol;
    theta::Real = 0.0, energy_upper::Real, verbose::Bool = false,
    eps_abs::Real = 1e-7, eps_rel::Real = 1e-7, max_iters::Integer = 50_000,
    solver::Symbol = :scs)
    problem = build_sdp_model(model; theta, verbose, eps_abs, eps_rel, max_iters, solver)
    jump_model = problem.jump_model
    rhoij = problem.rhoij
    @constraint(jump_model, problem.energy_expr <= Float64(energy_upper))

    moment_expr = trace_expr(spec.op, rhoij[spec.sites])
    if sense == :min
        @objective(jump_model, Min, moment_expr)
    elseif sense == :max
        @objective(jump_model, Max, moment_expr)
    else
        error("sense must be :min or :max")
    end
    optimize!(jump_model)
    return (;
        value = value(moment_expr),
        termination = string(termination_status(jump_model)),
        primal_status = string(primal_status(jump_model)),
    )
end

function moment_ranges_on_optimal_face(model; theta::Real = 0.0,
    energy_tol::Real = 1e-6, verbose::Bool = false,
    eps_abs::Real = 1e-7, eps_rel::Real = 1e-7, max_iters::Integer = 50_000,
    solver::Symbol = :scs)
    base = solve_sdp_observables(model; theta, verbose, eps_abs, eps_rel, max_iters, solver)
    energy_upper = base.energy + Float64(energy_tol)
    rows = NamedTuple[]
    for spec in model.selected
        lo = optimize_moment_on_face(model, spec, :min; theta, energy_upper, verbose,
            eps_abs, eps_rel, max_iters, solver)
        hi = optimize_moment_on_face(model, spec, :max; theta, energy_upper, verbose,
            eps_abs, eps_rel, max_iters, solver)
        push!(rows, (;
            moment = spec.name,
            sites = spec.sites,
            base_value = base.moments[length(rows) + 1],
            min_value = lo.value,
            max_value = hi.value,
            width = hi.value - lo.value,
            min_status = lo.termination,
            max_status = hi.termination,
        ))
    end
    return (; base_energy = base.energy, energy_upper, base_status = base.termination, rows)
end

function tr_expr_identity(X)
    return sum(X[i, i] for i in axes(X, 1))
end

function two_index(a::Int, b::Int, C::Int)
    return (a - 1) * C + b
end

function block_index(site::Int, α::Int, d::Int)
    return (site - 1) * d + α
end

# ---------------------------------------------------------------------------
# Finite-difference response and experiment helpers
# ---------------------------------------------------------------------------

function finite_response(f, s::Int; eps::Float64 = 1e-4)
    out0 = f(zeros(s))
    χ = zeros(length(out0), s)
    for a in 1:s
        hp = zeros(s)
        hm = zeros(s)
        hp[a] = eps
        hm[a] = -eps
        χ[:, a] = (f(hp) - f(hm)) ./ (2eps)
    end
    return χ
end

"""
    run_probe(...)

Run the small spinless smoke experiment. This is not the main reported model,
but it exercises the same physical-vs-SDP comparison pipeline.
"""
function run_probe(; L::Int = 4, theta_grid = [-1e-2, -5e-3, -2e-3, 2e-3, 5e-3, 1e-2],
    eps_h::Float64 = 1e-4, perturbation::Symbol = :density, selected::Symbol = :dimer_imbalance,
    eps_abs::Real = 1e-7, eps_rel::Real = 1e-7, max_iters::Integer = 50_000,
    solver::Symbol = :scs)
    model = dimer_model(; L, perturbation, selected)
    s = length(model.selected)

    phys0 = physical_observables(model; theta = 0.0)
    sdp0 = solve_sdp_observables(model; theta = 0.0, eps_abs, eps_rel, max_iters, solver)

    phys_response = finite_response(h -> physical_observables(model; theta = 0.0, fields = h).moments, s; eps = eps_h)
    sdp_response = finite_response(h -> solve_sdp_observables(model; theta = 0.0, fields = h, eps_abs, eps_rel, max_iters, solver).moments, s; eps = eps_h)

    mismatch_rows = NamedTuple[]
    for θ in theta_grid
        phys = physical_observables(model; theta = θ)
        sdp = solve_sdp_observables(model; theta = θ, eps_abs, eps_rel, max_iters, solver)
        δ = sdp.moments .- phys.moments
        push!(mismatch_rows, (theta = θ, delta_norm = norm(δ), ratio = norm(δ) / abs(θ),
            sdp_energy = sdp.energy, phys_energy = phys.energy,
            sdp_status = sdp.termination))
    end

    return (;
        model,
        reference = (;
            phys_energy = phys0.energy,
            sdp_energy = sdp0.energy,
            energy_gap = sdp0.energy - phys0.energy,
            physical_gap = phys0.gap,
            phys_moments = phys0.moments,
            sdp_moments = sdp0.moments,
            moment_error = norm(sdp0.moments .- phys0.moments),
            sdp_status = sdp0.termination,
        ),
        response = (;
            phys = phys_response,
            sdp = sdp_response,
            diff_norm = norm(sdp_response - phys_response),
        ),
        mismatch_rows,
    )
end

"""
    run_spinful_hubbard_probe(...)

Common driver for the spinful Hubbard scripts. It returns reference exactness,
response matching, and small-theta mismatch diagnostics in one named tuple.
"""
function run_spinful_hubbard_probe(; L::Int = 4, U::Real = 4.0,
    theta_grid = [-1e-2, -5e-3, -2e-3, 2e-3, 5e-3, 1e-2],
    eps_h::Float64 = 1e-4, perturbation::Symbol = :hopping,
    selected::Symbol = :charge_imbalance,
    stagger::Real = 0.0,
    eps_abs::Real = 1e-7, eps_rel::Real = 1e-7, max_iters::Integer = 50_000,
    solver::Symbol = :scs)
    model = spinful_hubbard_dimer_model(; L, U, stagger, perturbation, selected)
    s = length(model.selected)

    phys0 = physical_observables(model; theta = 0.0)
    sdp0 = solve_sdp_observables(model; theta = 0.0, eps_abs, eps_rel, max_iters, solver)

    phys_response = finite_response(h -> physical_observables(model; theta = 0.0, fields = h).moments, s; eps = eps_h)
    sdp_response = finite_response(h -> solve_sdp_observables(model; theta = 0.0, fields = h, eps_abs, eps_rel, max_iters, solver).moments, s; eps = eps_h)

    mismatch_rows = NamedTuple[]
    for θ in theta_grid
        phys = physical_observables(model; theta = θ)
        sdp = solve_sdp_observables(model; theta = θ, eps_abs, eps_rel, max_iters, solver)
        δ = sdp.moments .- phys.moments
        push!(mismatch_rows, (theta = θ, delta_norm = norm(δ), ratio = norm(δ) / abs(θ),
            sdp_energy = sdp.energy, phys_energy = phys.energy,
            sdp_status = sdp.termination))
    end

    return (;
        model,
        reference = (;
            phys_energy = phys0.energy,
            sdp_energy = sdp0.energy,
            energy_gap = sdp0.energy - phys0.energy,
            physical_gap = phys0.gap,
            phys_moments = phys0.moments,
            sdp_moments = sdp0.moments,
            moment_error = norm(sdp0.moments .- phys0.moments),
            sdp_status = sdp0.termination,
        ),
        response = (;
            phys = phys_response,
            sdp = sdp_response,
            diff_norm = norm(sdp_response - phys_response),
        ),
        mismatch_rows,
    )
end

end
