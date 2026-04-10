"""
Helpers for marker-annotation setup.

Annotation handling in JWAS intentionally has two phases:

1. `get_genotypes` validates and stores the raw annotation design information and
   the user-supplied starting `Pi` without assuming how many traits the genotype
   term will affect.
2. `build_model` finalizes any annotation state that depends on the model, such
   as the number of traits and the joint inclusion states used by multi-trait
   BayesC.

This keeps the user-facing API consistent across single-trait and multi-trait
annotated methods: the same annotation matrix and starting `Pi` always go into
`get_genotypes`, while only model-dependent internal state is deferred.
"""

const ANNOTATED_BAYESC_MT_STATES = ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0))

annotated_bayesc_mt_default_row() = Float64[0.0, 0.0, 0.0, 1.0]

function annotated_bayesc_mt_state_keys()
    return [Float64[state...] for state in ANNOTATED_BAYESC_MT_STATES]
end

function annotated_bayesc_mt_state_index(state::AbstractVector{<:Real})
    length(state) == 2 || error("Annotated multi-trait BayesC v1 expects 2-trait state labels.")
    d1 = Int(round(Float64(state[1])))
    d2 = Int(round(Float64(state[2])))
    if d1 == 0 && d2 == 0
        return 1
    elseif d1 == 1 && d2 == 0
        return 2
    elseif d1 == 0 && d2 == 1
        return 3
    elseif d1 == 1 && d2 == 1
        return 4
    end
    error("Annotated multi-trait BayesC v1 expects binary 2-trait state labels.")
end

function annotated_bayesc_mt_row_from_dict(pi::AbstractDict)
    row = zeros(Float64, 4)
    for (state, prob) in pi
        row[annotated_bayesc_mt_state_index(state)] = Float64(prob)
    end
    isapprox(sum(row), 1.0; atol=1e-8) || error("Summation of probabilities of Pi is not equal to one.")
    return row
end

function annotated_bayesc_mt_pi_dict(row::AbstractVector{<:Real})
    length(row) == 4 || error("Annotated multi-trait BayesC v1 expects four joint prior probabilities.")
    probs = Float64.(row)
    keys = annotated_bayesc_mt_state_keys()
    return Dict(keys[i] => probs[i] for i in eachindex(keys))
end

function validate_bayesc_mt_start_row(start_row::AbstractVector{<:Real})
    length(start_row) == 4 || error("Annotated multi-trait BayesC v1 expects four joint prior probabilities.")
    trait1_mass = Float64(start_row[2]) + Float64(start_row[4])
    trait2_mass = Float64(start_row[3]) + Float64(start_row[4])
    shared_mass = Float64(start_row[4])

    trait1_mass > 0.0 || error("Annotated multi-trait BayesC requires positive startup prior mass in states {10,11} for trait 1.")
    trait2_mass > 0.0 || error("Annotated multi-trait BayesC requires positive startup prior mass in states {01,11} for trait 2.")
    shared_mass > 0.0 || error("Annotated multi-trait BayesC requires positive startup prior mass in shared state 11.")
    return nothing
end

function initialize_bayesc_single_trait_annotations!(genotypei::Genotypes)
    if genotypei.annotations === false || genotypei.method != "BayesC" || genotypei.ntraits != 1
        return nothing
    end

    raw_pi = genotypei.annotation_start_pi
    raw_pi isa AbstractDict && error("Annotated BayesC genotypes initialized with a joint Pi dictionary cannot be rebuilt for a single-trait model. Use a fresh get_genotypes call with scalar/vector Pi for single-trait analysis.")

    start_pi = if raw_pi isa AbstractVector
        length(raw_pi) == genotypei.nMarkers || error("Annotated BayesC starting Pi vector length $(length(raw_pi)) must match the number of markers ($(genotypei.nMarkers)).")
        collect(Float64, raw_pi)
    else
        fill(Float64(raw_pi), genotypei.nMarkers)
    end

    design_matrix = genotypei.annotations.design_matrix
    genotypei.annotations = MarkerAnnotations(design_matrix)
    genotypei.π = start_pi
    return nothing
end

function initialize_bayesc_mt_annotations!(genotypei::Genotypes)
    if genotypei.annotations === false || genotypei.method != "BayesC" || genotypei.ntraits == 1
        return nothing
    end
    genotypei.ntraits == 2 || error("Annotated multi-trait BayesC currently supports exactly 2 traits.")

    raw_pi = genotypei.annotation_start_pi
    start_row = if raw_pi isa AbstractDict
        annotated_bayesc_mt_row_from_dict(raw_pi)
    elseif raw_pi == 0.0
        # Preserve the legacy multi-trait BayesC default: before seeing data,
        # all markers start in the all-traits-active state when Pi is omitted.
        annotated_bayesc_mt_default_row()
    elseif raw_pi isa AbstractVector{<:Real} &&
           length(raw_pi) == genotypei.nMarkers &&
           all(iszero, Float64.(raw_pi))
        annotated_bayesc_mt_default_row()
    else
        error("Annotated multi-trait BayesC requires Pi=0.0 or a joint Pi dictionary.")
    end
    validate_bayesc_mt_start_row(start_row)

    design_matrix = genotypei.annotations.design_matrix
    coeffs = zeros(Float64, size(design_matrix, 2), 3)
    genotypei.annotations = MarkerAnnotations(
        design_matrix;
        nsteps=3,
        nclasses=4,
        coefficients=coeffs,
        snp_pi=repeat(reshape(start_row, 1, :), genotypei.nMarkers, 1),
    )
    genotypei.π = annotated_bayesc_mt_pi_dict(start_row)
    return nothing
end

"""
    finalize_marker_annotation_setup!(genotypei)

Finalize any annotation state that depends on the model rather than the raw
genotype input.

`get_genotypes` handles annotation validation and stores the raw design
information. This finalizer is called from `build_model` after `genotypei.ntraits`
is known, which is necessary for methods such as annotated multi-trait BayesC
whose nested probit structure and joint-state labels depend on the modeled trait
count.
"""
function finalize_marker_annotation_setup!(genotypei::Genotypes)
    if genotypei.annotations === false
        return nothing
    end

    if genotypei.method == "BayesC"
        if genotypei.ntraits == 1
            initialize_bayesc_single_trait_annotations!(genotypei)
        elseif genotypei.ntraits > 1
            initialize_bayesc_mt_annotations!(genotypei)
        end
    end
    return nothing
end
