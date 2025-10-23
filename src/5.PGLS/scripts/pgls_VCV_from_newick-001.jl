using LinearAlgebra



function read_newick(path::AbstractString)
    open(path, "r") do io
        return read(io, String)
    end
end


# newick_str = read_newick("/Users/dannox/Desktop/JWAS-PGLS-branch/data/HAL.newick")
newick_str = "(
  (Human:0.8, Chimpanzee:0.7):0.25,
  (
    (Mouse:0.6, Rat:0.5):0.35,
    (Cow:0.7, Pig:0.6):0.25
  ):0.1,
  (
    (Dog:0.55, Cat:0.45):0.4,
    (
      (Horse:0.6, Donkey:0.55):0.3,
      (Elephant:0.75, Manatee:0.65):0.25
    ):0.15
  ):0.1,
  (
    (Chicken:0.9, Duck:0.85):0.1,
    (Frog:0.7, Lizard:0.6):0.2
  ):0.05
):0.0;
  "





# ---------- Minimal tree structs ----------
struct Node
    name::String        # name for leaves and for internal nodes
    parent::Int         # index of parent node (0 for root)
    children::Vector{Int} # indices of node's children
    bl::Float64            # branch length from parent to this node (root has 0)
end



# ---------- Newick parser (labels + branch lengths) ----------
function parse_newick(s::AbstractString)
    s = strip(s)
    if !isempty(s) && s[end] == ';'
        s = s[1:end-1]
    end
    nodes = Node[]               # 1-based indices
    function push_node(name, parent, bl)
        push!(nodes, Node(name, parent, Int[], bl))
        return length(nodes)
    end

    stack = Int[]                # stack of current internal nodes
    current = 0
    i, n = 1, lastindex(s)
    last_closed = 0

    read_label(i) = begin
        j = i
        while j <= n && !(s[j] in (',','(',')',':'))
            j += 1
        end
        strip(s[i:j-1]), j
    end
    read_float(i) = begin
        j = i
        # allow digits plus + - . e E (scientific notation)
        while j <= n && (isdigit(s[j]) || s[j] in ('+','-','e','E','.'))
            j += 1
        end
        val = parse(Float64, s[i:j-1])
        return val, j
    end

    # Create a dummy root to attach the first internal
    root = push_node("", 0, 0.0)
    current = root

    while i <= n
        while i <= n && isspace(s[i])
            i += 1
        end
        i > n && break
        c = s[i]
        if c == '('
            # start a new internal node under current
            nd = push_node("", current, 0.0)
            push!(nodes[current].children, nd)
            push!(stack, nd)
            current = nd
            i += 1
        elseif c == ')'
            # close current internal
            last_closed = current
            pop!(stack)
            current = isempty(stack) ? root : stack[end]
            i += 1
            # Skip whitespace before optional label/length
            while i <= n && isspace(s[i])
                i += 1
            end
            # optional internal label
            if i <= n && !(s[i] in (',',')',':'))
                label, i = read_label(i)
                nodes[last_closed] = Node(label, nodes[last_closed].parent, nodes[last_closed].children, nodes[last_closed].bl)
            end
            # optional branch length for this internal node
            if i <= n && s[i] == ':'
                i += 1
                bl, i = read_float(i)
                nodes[last_closed] = Node(nodes[last_closed].name, nodes[last_closed].parent, nodes[last_closed].children, bl)
            end
        elseif c == ','
            i += 1
        else
            # leaf: Name[:bl]
            label, i = read_label(i)
            bl = 0.0
            if i <= n && s[i] == ':'
                i += 1
                bl, i = read_float(i)
            end
            nd = push_node(label, current, bl)
            push!(nodes[current].children, nd)
        end
    end

    # find true root (the unique node with parent==0) and set its bl=0
    roots = [k for k in eachindex(nodes) if nodes[k].parent == 0]
    length(roots) == 1 || error("Newick parse error: expected 1 root, found $(length(roots)).")
    r = roots[1]
    nodes[r] = Node(nodes[r].name, nodes[r].parent, nodes[r].children, 0.0)
    return nodes, r
end

# ---------- Utilities ----------
tip_indices(nodes) = [i for i in eachindex(nodes) if isempty(nodes[i].children)]
tip_names(nodes, tips) = [nodes[i].name for i in tips]

function depths_from_root(nodes, root)
    depth = zeros(Float64, length(nodes))
    stack = [root]
    while !isempty(stack)
        v = pop!(stack)
        for c in nodes[v].children
            depth[c] = depth[v] + nodes[c].bl
            push!(stack, c)
        end
    end
    depth
end

function ancestors(nodes, v)
    a = Int[]
    while v != 0
        push!(a, v)
        v = nodes[v].parent
    end
    a
end

function mrca(nodes, a, b)
    A = Set(ancestors(nodes, a))
    for v in ancestors(nodes, b)
        if v in A
            return v
        end
    end
    error("No MRCA found (tree must be connected).")
end

# ---------- Build BM VCV ----------
"""
    vcv_bm(newick_str; sigma2=1.0)

Return (tips, Σ, h):
- tips: Vector{String} — tip names in the order used
- Σ:    Matrix{Float64} — BM VCV (rate sigma2)
- h:    Vector{Float64} — root→tip distances (variances when sigma2=1)
"""
function vcv_bm(newick_str::AbstractString; sigma2::Float64=1.0)
    nodes, root = parse_newick(newick_str)
    depth = depths_from_root(nodes, root)
    tips = tip_indices(nodes)
    names = tip_names(nodes, tips)
    n = length(tips)

    # Pairwise MRCA depths => covariance
    Σ = zeros(Float64, n, n)
    for i in 1:n
        ti = tips[i]
        Σ[i,i] = depth[ti]
        for j in i+1:n
            tj = tips[j]
            a = mrca(nodes, ti, tj)
            Σ[i,j] = depth[a]
            Σ[j,i] = depth[a]
        end
    end
    Σ .*= sigma2
    return names, Σ, depth[tips]
end

nodes, root = parse_newick(newick_str)
nodes
root
# # ---------- Example usage ----------
# if abspath(PROGRAM_FILE) == @__FILE__
#     # Replace this with your file:
#     # newick_str = read("/path/to/your/tree.newick", String)
#     newick_str = "(
#         (Human:0.8, Chimpanzee:0.7):0.25,
#         (
#           (Mouse:0.6, Rat:0.5):0.35,
#           (Cow:0.7, Pig:0.6):0.25
#         ):0.1,
#         (
#           (Dog:0.55, Cat:0.45):0.4,
#           (
#             (Horse:0.6, Donkey:0.55):0.3,
#             (Elephant:0.75, Manatee:0.65):0.25
#           ):0.15
#         ):0.1,
#         (
#           (Chicken:0.9, Duck:0.85):0.1,
#           (Frog:0.7, Lizard:0.6):0.2
#         ):0.05
#       ):0.0;
#         "

#     tips, Σ, h = vcv_bm(newick_str; sigma2=1.0)
#     println("Tip order: ", join(tips, ", "))
#     println("Σ (BM VCV):")
#     show(stdout, "text/plain", Σ); println()
#     println("Diag vs heights (should match):")
#     for i in 1:length(tips)
#         println(rpad(tips[i], 8), "  Vii=$(Σ[i,i])  hi=$(h[i])")
#     end

#     # Optional correlation matrix
#     sd = sqrt.(diag(Σ))
#     Corr = Σ ./ (sd * sd')
#     println("\nCorrelation:")
#     show(stdout, "text/plain", Corr); println()
# end

newick_str = "(
        (Human:0.8, Chimpanzee:0.7):0.25,
        (
          (Mouse:0.6, Rat:0.5):0.35,
          (Cow:0.7, Pig:0.6):0.25
        ):0.1,
        (
          (Dog:0.55, Cat:0.45):0.4,
          (
            (Horse:0.6, Donkey:0.55):0.3,
            (Elephant:0.75, Manatee:0.65):0.25
          ):0.15
        ):0.1,
        (
          (Chicken:0.9, Duck:0.85):0.1,
          (Frog:0.7, Lizard:0.6):0.2
        ):0.05
      ):0.0;
        "
tips, V, h = vcv_bm(newick_str; sigma2=1.0)

