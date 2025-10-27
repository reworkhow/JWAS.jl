using CSV, DataFrames



function read_newick(path::AbstractString)
    open(path, "r") do io
        return read(io, String)
    end
end


newick_str = read_newick("/Users/dannox/Desktop/JWAS-PGLS-branch/data/HAL.newick")
# newick_str = "((Human:0.8,Chimpanzee:0.7)Primates:0.25,((Mouse:0.6,Rat:0.5)Pests:0.35,(Cow:0.7,Pig:0.6):0.25)Farm:0.1,((Dog:0.55,Cat:0.45)Pets:0.4,((Horse:0.6,Donkey:0.55)Farming:0.3,(Elephant:0.75,Manatee:0.65)Wild:0.25):0.15):0.1,((Chicken:0.9,Duck:0.85)Food:0.1,(Frog:0.7,Lizard:0.6)Rawr:0.2):0.05)root;"
# newick_str = "((Dog:0.55,Cat:0.45)Pets:0.4,((Horse:0.6,Donkey:0.55)Farming:0.3,(Elephant:0.75,Manatee:0.65)Wild:0.25):0.15):0.1,((Chicken:0.9,Duck:0.85)Food:0.1)root;"
# species_in_pheno = ["Human", "Chimpanzee", "Mouse", "Rat", "Cow", "Pig", "Dog", "Cat", "Horse", "Donkey", "Elephant", "Manatee", "Chicken", "Duck"]

# Load species in phenotype data to intersect species
pheno = "/Users/dannox/Desktop/JWAS/segway_orthologs/brain_size_phenotypes.csv"
phenotypes = CSV.read(pheno, DataFrame, delim = ',', header=true,missingstrings=["NA"])
species_in_pheno = phenotypes[:,1]


############################
# parsing state / storage  #
############################

# We'll parse the Newick manually and build parallel arrays that
# describe each node in the tree (internal nodes and leaves).
# Each node gets:
#   - a unique node_id
#   - a label (taxon name or "(internal)" or "root")
#   - its parent node_id (or missing if it's the actual root)
#   - its branch length to that parent (Float64 or missing)

node_id = 0
isalnum(c::Char) = isletter(c) || isdigit(c)  # helper to read names

node_ids        = Int[]                                      # node_id
node_labels     = String[]                                   # node's label
node_startpos   = Int[]                                      # index in Newick string where node's label/group starts
node_parent     = Vector{Union{Int,Missing}}()               # parent node_id
node_branchlen  = Vector{Union{Float64,Missing}}()           # branch length to parent

internal_stack = Int[]  # stack of "currently open" internal nodes as we descend with '('
last_closed_internal = 0  # node_id of most recently closed internal node after ')'
expecting_label_for_closed = false      # after ")", we might see a label, e.g. ")root"
expecting_length_for_closed = false     # after ")", we might see ":0.25" giving that node's branch length

#############
# helpers   #
#############

# push_node!: create a new node in our arrays.
# label         = "(internal)" for unnamed clades, "Dog" for leaves, "root" for final root
# start_pos     = character index in the Newick string
# parent_id     = node_id of this node's parent, or missing
# brlen         = branch length from parent to this node, or missing if not known yet
function push_node!(label::String, start_pos::Int;
                    parent_id::Union{Int,Missing},
                    brlen::Union{Float64,Missing})
    global node_id
    node_id += 1
    nid = node_id
    push!(node_ids, nid)
    push!(node_labels, label)
    push!(node_startpos, start_pos)
    push!(node_parent, parent_id)
    push!(node_branchlen, brlen)
    return nid
end

# set_branchlen!: once we know the branch length for an internal node
# (which arrives after the ')'), retroactively fill it in.
function set_branchlen!(nid::Int, brlen::Float64)
    for idx in eachindex(node_ids)
        if node_ids[idx] == nid
            node_branchlen[idx] = brlen
            return
        end
    end
end

# label_internal!: after we close an internal node with ")"
# the next token might be its label, e.g. ")root"
function label_internal!(nid::Int, newlabel::String)
    for idx in eachindex(node_ids)
        if node_ids[idx] == nid
            node_labels[idx] = newlabel
            return
        end
    end
end

# new_internal!: called when we see "(".
# That means "start a new internal clade".
# We create a node for that clade, parented to the node that's currently
# on top of the stack (if any), then push it onto the stack.
function new_internal!(where)
    parent_id = isempty(internal_stack) ? missing : internal_stack[end]
    nid = push_node!("(internal)", where; parent_id=parent_id, brlen=missing)
    push!(internal_stack, nid)
    return nid
end

# new_leaf!: called when we see a taxon label (like "Dog") in the string.
# We also grab its branch length (like ":0.55") if present.
function new_leaf!(label, where, brlen_to_parent::Union{Float64,Missing})
    parent_id = isempty(internal_stack) ? missing : internal_stack[end]
    nid = push_node!(label, where; parent_id=parent_id, brlen=brlen_to_parent)
    return nid
end

##################
# main parse     #
##################

# We now walk left-to-right through the Newick string and
# build the tree incrementally.

i = 1
while i <= lastindex(newick_str)
    ch = newick_str[i]

    if ch == '('
        # "(" means: start a new internal node (a subtree/clade).
        new_internal!(i)
        expecting_label_for_closed  = false
        expecting_length_for_closed = false
        last_closed_internal = 0
        i += 1

    elseif ch == ')'
        # ")" means: we just finished the current internal node.
        # Pop it off the stack so new children won't attach to it.
        last_closed_internal = pop!(internal_stack)
        # After ")", Newick can attach:
        #   - a branch length for that entire clade ":0.25"
        #   - a label for that clade "root"
        # so we set these flags and handle them in the next tokens.
        expecting_label_for_closed  = true
        expecting_length_for_closed = true
        i += 1

    elseif ch == ':'
        # ":" always starts a branch length number, e.g. ":0.55"
        # We parse the number manually.
        i += 1
        num_start = i
        while i <= lastindex(newick_str) &&
              (isdigit(newick_str[i]) || newick_str[i] in ".eE-+")
            i += 1
        end
        brlen_val = parse(Float64, newick_str[num_start:i-1])

        # This branch length belongs to either:
        #   - a just-closed internal node (after ")")
        #   - or a leaf (handled in the leaf branch below)
        if expecting_length_for_closed && last_closed_internal != 0
            set_branchlen!(last_closed_internal, brlen_val)
            expecting_length_for_closed = false
        end

    elseif ch == ',' || ch == ';'
        # "," separates siblings.
        # ";" ends the tree.
        # After a comma/semicolon, we can't still be naming/lengthening
        # a node that was just closed by ")"
        expecting_label_for_closed  = false
        expecting_length_for_closed = false
        last_closed_internal = 0
        i += 1

    elseif isletter(ch)
        # We hit the start of a label (taxon name or internal node name).
        # Consume alphanumeric/underscore characters to get the whole label.
        start_i = i
        while i <= lastindex(newick_str) &&
              (isalnum(newick_str[i]) || newick_str[i] == '_')
            i += 1
        end
        label = newick_str[start_i:i-1]

        if expecting_label_for_closed && last_closed_internal != 0
            # This label (e.g. "root") names the internal node we just closed.
            label_internal!(last_closed_internal, label)
            expecting_label_for_closed  = false
        else
            # Otherwise this is a leaf/tip label (e.g. "Dog").
            # A tip label may be immediately followed by a ":<number>" branch len,
            # so grab it if present.
            leaf_brlen = missing
            if i <= lastindex(newick_str) && newick_str[i] == ':'
                i += 1
                num_start = i
                while i <= lastindex(newick_str) &&
                      (isdigit(newick_str[i]) || newick_str[i] in ".eE-+")
                    i += 1
                end
                leaf_brlen = parse(Float64, newick_str[num_start:i-1])
            end
            new_leaf!(label, start_i, leaf_brlen)

            expecting_label_for_closed  = false
            expecting_length_for_closed = false
            last_closed_internal = 0
        end

    else
        # Anything else (shouldn't really happen in clean Newick),
        # just advance.
        i += 1
    end
end

##########################################
# fix parenting for early-created clades #
##########################################

# At parse time, we assign parent-child relationships as we go.
# BUT we may create an internal node before we've actually seen its
# eventual true root label "root". After parsing, some nodes may still
# have parent_id = missing even though they're not really the root.
#
# Strategy:
#   1. Find the real root (the node labeled "root", or the first missing parent).
#   2. Any other node that still has missing parent gets re-parented to that root.

root_candidates = findall(i -> node_labels[i] == "root", eachindex(node_labels))
true_root_idx = !isempty(root_candidates) ? root_candidates[1] :
                findfirst(i -> ismissing(node_parent[i]), eachindex(node_parent))
true_root_id = node_ids[true_root_idx]

for j in eachindex(node_parent)
    if j != true_root_idx && ismissing(node_parent[j])
        node_parent[j] = true_root_id
    end
end

####################
# build df_nodes   #
####################

# df_nodes is our final per-node table.
# Columns:
#   node_id
#   label
#   start_index         (for debugging / trace-back into Newick string)
#   parent_id
#   is_child            (0 for root, 1 otherwise)
#   branchlen_to_parent
df_nodes = begin
    is_child_vec = [ismissing(p) ? 0 : 1 for p in node_parent]
    DataFrame(
        node_id              = node_ids,
        label                = node_labels,
        start_index          = node_startpos,
        parent_id            = node_parent,
        is_child             = is_child_vec,
        branchlen_to_parent  = node_branchlen,
    )
end

#########################################
# compute dist_from_root for each node  #
#########################################

# dist_from_root for node n = sum of branchlen_to_parent along the path
# from the true root down to n.
#
# We'll compute it recursively with memoization.

parent_of = Dict(df_nodes.node_id .=> df_nodes.parent_id)
brlen_of  = Dict(df_nodes.node_id .=> df_nodes.branchlen_to_parent)

dist_from_root_cache = Dict{Int,Float64}()

function dist_from_root(nid::Int)::Float64
    # If we've already computed this node's distance, reuse it.
    if haskey(dist_from_root_cache, nid)
        return dist_from_root_cache[nid]
    end

    # Base case: if this node has no parent (it's the root),
    # its distance from the root is 0.0.
    p = parent_of[nid]
    d = if ismissing(p)
        0.0
    else
        # Otherwise: parent's distance + this edge length.
        edge_len = brlen_of[nid]
        dist_from_root(p) + (ismissing(edge_len) ? 0.0 : edge_len)
    end

    dist_from_root_cache[nid] = d
    return d
end

# Add dist_from_root as a column in df_nodes for convenience.
df_nodes.dist_from_root = [dist_from_root(nid) for nid in df_nodes.node_id]

#########################################
# build species x species distance mat  #
#########################################

# Goal: make an N x N matrix of pairwise patristic distances
# between all *tips* ("leaves", i.e. real species), not internal nodes.
#
# Patristic distance between A and B is:
#     dist(A,B) = dist_root(A) + dist_root(B) - 2*dist_root(LCA(A,B))
#
# We'll:
#   1. define which nodes are "tips"
#   2. build parent/child maps
#   3. write helpers to find ancestors, LCA, and patristic distance
#   4. fill the matrix
#
# Note: we skip "(internal)" and "root" when building tip list.

# define tips (actual leaves, not labeled internal clades)
children = Dict(n => Int[] for n in df_nodes.node_id)
for (child, parent) in zip(df_nodes.node_id, df_nodes.parent_id)
    if !ismissing(parent)
        push!(children[parent], child)
    end
end

is_leaf = [isempty(children[n]) for n in df_nodes.node_id]
species_mask   = is_leaf .& (df_nodes.label .!= "root")
species_labels = df_nodes.label[species_mask]
species_ids    = df_nodes.node_id[species_mask]   # their node_ids

# build fresh dicts for THIS tree
parent_map = Dict(df_nodes.node_id .=> df_nodes.parent_id)
dist_map   = Dict(df_nodes.node_id .=> df_nodes.dist_from_root)

# ancestors(n, parent_map): return a vector of all ancestors of n,
# starting with n itself, then its parent, up to the root.
function ancestors(nid::Int, parent_map)
    anc = Int[]
    cur = nid
    while true
        push!(anc, cur)
        p = parent_map[cur]
        if ismissing(p)
            break
        end
        cur = p
    end
    return anc
end

# lca(a,b,...): find the "lowest common ancestor" of a and b.
# We take the intersection of ancestor sets, and choose the
# one with the largest dist_from_root (i.e. deepest in the tree).
function lca(n1::Int, n2::Int, parent_map, dist_map)
    a1 = ancestors(n1, parent_map)
    a2 = ancestors(n2, parent_map)
    a2set = Set(a2)

    best = nothing
    bestd = -Inf

    for a in a1
        if (a in a2set)
            d = dist_map[a]
            if d > bestd
                best = a
                bestd = d
            end
        end
    end

    return best
end

# patristic distance between two tip node_ids, using:
#   dist(A)+dist(B)-2*dist(LCA(A,B))
# we take abs() just in case of tiny floating negs like -1e-16.
function patristic(n1::Int, n2::Int, parent_map, dist_map)
    n1 == n2 && return 0.0
    a = lca(n1, n2, parent_map, dist_map)
    d = dist_map[n1] + dist_map[n2] - 2 * dist_map[a]
    return abs(d)
end

# Build N x N matrix of pairwise distances for all species.
# but set the diagonal to each tip's root-to-tip distance.
n = length(species_ids)
tip_root = [dist_map[species_ids[i]] for i in 1:n]

distmat = Matrix{Float64}(undef, n, n)
for ii in 1:n, jj in 1:n
    if ii == jj
        distmat[ii, jj] = tip_root[ii]            # diagonal: root→tip
    else
        distmat[ii, jj] = patristic(species_ids[ii], species_ids[jj], parent_map, dist_map)
    end
end

# Turn that matrix into a labeled DataFrame (rows+cols = species).
df_dist = DataFrame(distmat, :auto)
rename!(df_dist,
    Dict(Symbol("x$i") => Symbol(species_labels[i]) for i in 1:n)
)
df_dist[!, :Species] = species_labels
select!(df_dist, [:Species; Symbol.(species_labels)...])


# Pick species in the order they appear in the phenotype file
keep = [s for s in species_in_pheno if s in df_dist.Species]

# Build a lookup: species name -> row index in df_dist
row_idx = Dict(df_dist.Species .=> eachindex(df_dist.Species))

# Reorder rows to match `keep`, and reorder columns the same way
row_order = [row_idx[s] for s in keep]
df_dist_pheno = df_dist[row_order, [:Species; Symbol.(keep)...]]


############
# outputs  #
############

println("df_nodes:")
show(df_nodes, allrows=true, allcols=true)

println("\n\ndf_dist (pairwise patristic distances among tips):")
show(df_dist, allrows=true, allcols=true)

# Save outputs as TSV files
CSV.write("/Users/dannox/Desktop/JWAS-PGLS-branch/output/parsed_nodes.tsv", df_nodes, delim='\t')
CSV.write("/Users/dannox/Desktop/JWAS-PGLS-branch/output/pairwise_distances.tsv", df_dist_pheno, delim='\t')




# Testing distance function
# function dist_by_label(a::String, b::String)
#     ida = species_ids[findfirst(==(a), species_labels)]
#     idb = species_ids[findfirst(==(b), species_labels)]
#     return patristic(ida, idb, parent_map, dist_map)
# end

# println("Human–Mouse distance = ", dist_by_label("Human", "Mouse"))
# should print 2.1
