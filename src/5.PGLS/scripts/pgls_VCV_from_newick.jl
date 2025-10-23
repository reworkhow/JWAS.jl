

#   Optional: ENV["SIGMA2"]="1.0" to rescale V by sigma^2

using Phylo  # for parsing Newick and tree operations
using LinearAlgebra
using CSV
using DataFrames
using Plots

function read_newick(path::AbstractString)
    open(path, "r") do io
        return read(io, String)
    end
end

newick_str = read_newick("/Users/dannox/Desktop/JWAS-PGLS-branch/data/HAL.newick")
ree = parsenewick(newick_str)
tips = getleafnames(ree)
h = heightstoroot(ree)
D = distances(ree)
test = h[1]
tt = h[2]
deest= D[1,2]


function vcv_from_tree(tree; sigma2::Float64=1.0)
    tips = getleafnames(tree)      # tip order
    h    = heightstoroot(tree)     # root-to-tip heights (aligned with tips)
    D    = distances(tree)         # patristic distances between tips
    n    = length(tips)
    V    = Matrix{Float64}(undef, n, n)
    @inbounds for i in 1:n, j in 1:n
        V[i, j] = 0.5 * (h[i] + h[j] - D[i, j])   # Cov(i,j) = height(MRCA(i,j))
    end
    V .*= sigma2
    V = 0.5 .* (V .+ V')           # enforce symmetry
    return tips, V, h
end

function save_v_and_tips(outprefix::AbstractString, tips::Vector{String}, V::AbstractMatrix{<:Real})
    df = DataFrame(V, :auto)
    # rename columns to tips
    for i in 1:length(tips)
        rename!(df, names(df)[i] => tips[i])
    end
    insertcols!(df, 1, :tip => tips)
    CSV.write(outprefix * "_VCV.csv", df)

    open(outprefix * "_tips.txt", "w") do io
        for t in tips
            println(io, t)
        end
    end
end

function plot_heatmap(outprefix::AbstractString, tips::Vector{String}, V::AbstractMatrix{<:Real})
    n = length(tips)
    showlabels = n ≤ 60
    xt = showlabels ? (1:n, tips) : :auto
    yt = showlabels ? (1:n, tips) : :auto
    plt = heatmap(V;
        xlabel = "Tips", ylabel = "Tips",
        title = "Phylogenetic VCV heatmap (n=$(n))",
        xticks = xt, yticks = yt,
        colorbar = true, aspect_ratio = :equal
    )
    savefig(plt, outprefix * "_VCV_heatmap.png")
end

function main()
    # if length(ARGS) != 2
    #     error("Usage: julia make_vcv_save_and_plot.jl /path/to/tree.newick /path/to/out/prefix")
    # end
    # newick_path, outprefix = ARGS
    # newick_str = read_newick(newick_path)
    newick_str = read_newick("/Users/dannox/Desktop/QTLrocks/EGP/scripts/PGLS/HAL.newick")
    outprefix = "/Users/dannox/Desktop/JWAS-PGLS-branch/output"
    tree = parsenewick(newick_str)

    sigma2 = tryparse(Float64, get(ENV, "SIGMA2", "1.0"))
    sigma2 === nothing && error("ENV[\"SIGMA2\"] must be a Float64 if set.")

    tips, V, h = vcv_from_tree(tree; sigma2 = sigma2)

    # Diagnostics
    min_eig = minimum(eigvals(Symmetric(V)))
    println("Tips: ", length(tips), " | VCV size: ", size(V), " | sigma^2 = ", sigma2)
    println("issymmetric(V): ", issymmetric(V), " | min eigenvalue: ", min_eig)  # Should be symmetric
    println("Diag vs heights (first 3): ", [(V[i,i], h[i]) for i in 1:min(3, length(tips))]) # Diagonal should match root-to-tip heights

    # Save outputs
    save_v_and_tips(outprefix, tips, V)
    plot_heatmap(outprefix, tips, V)
    println("Saved: $(outprefix)_VCV.csv")
    println("Saved: $(outprefix)_tips.txt")
    println("Saved: $(outprefix)_VCV_heatmap.png")
end

main()
