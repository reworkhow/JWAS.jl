using JWAS, DataFrames, CSV, Statistics, Random, ROC, Plots, HypothesisTests, Printf

# ----------------------------- Config -----------------------------------------
# Folder where your splitter saved files (from your save script)
splits_dir   = "/Users/dannox/Desktop/JWAS/irene_data/data/CV_20fold_80-20splits"
# Choose which trait column (from the saved *phenos* CSVs) to model:
trait        = "Brain_resid"
tissue       = "Cortex"
# trait        = "vlStatus"
# tissue       = "Liver"
# trait_type   = "B"  # binary trait
# trait_type   = "C"  # continuous trait
# Where to write JWAS fold outputs:
outdir_root  = "/Users/dannox/Desktop/JWAS/irene_data/$(tissue)_$(trait)_20fold_80-20Split_500-permutations-001"

# MCMC settings
chain_length = 10_000
save_every   = 50
seed         = 55
burn_in      = 1_000


permute_n = 500

# ----------------------------- Helpers ----------------------------------------
# Read a fold’s files produced by your splitter
function read_fold_tables(splits_dir::AbstractString, fold::Int, trait::AbstractString)
    tag = string(fold)
    ph_train = CSV.read(joinpath(splits_dir, "fold_$(tag)_train_phenos.csv"), DataFrame)
    ph_train = ph_train[!, [:ID, Symbol(trait)]]
    ph_test  = CSV.read(joinpath(splits_dir, "fold_$(tag)_test_phenos.csv"),  DataFrame)
    ph_test  = ph_test[!, [:ID, Symbol(trait)]]
    g_train  = CSV.read(joinpath(splits_dir, "fold_$(tag)_$(tissue)_train_genotypes.csv"), DataFrame)
    g_test   = CSV.read(joinpath(splits_dir, "fold_$(tag)_$(tissue)_test_genotypes.csv"),  DataFrame)

    # # Combine genotypes so JWAS can (often) emit EBVs for both train+test IDs
    # geno_all = vcat(g_train, g_test)
    # # Deduplicate if any overlap
    # geno_all = unique(geno_all, :ID)

    return ph_train, ph_test, g_train, g_test
end

# ----------------------------- Discover #folds --------------------------------
function discover_nfold(splits_dir::AbstractString)
    files = readdir(splits_dir)
    # accept padded or unpadded: fold_1_*, fold_01_*, fold_20_*, etc.
    pat = r"^fold_(\d+)_train_phenos\.csv$"
    matches = filter(f -> occursin(pat, f), files)
    @assert !isempty(matches) "No fold_*_train_phenos.csv files found in $(splits_dir)"
    nums = [parse(Int, match(pat, f).captures[1]) for f in matches]
    return maximum(nums)
end




# ----------------------------- Run all folds, continuous trait ----------------------------------
nfold = discover_nfold(splits_dir)
# println("Detected nfold = ", nfold)
all_acc_np = Float64[]
all_h2_np = Float64[]
all_acc_p = Float64[]
all_h2_p = Float64[]
mkdir(outdir_root)
for fold in 1:nfold
    tag = Int(fold)
    println("=== Fold $fold / $nfold ===")

    ph_train, ph_test, g_train, g_test = read_fold_tables(splits_dir, tag, trait)

    global genotypes = get_genotypes(g_train; method="BayesC", center=false, quality_control=false)

    model = build_model("$(trait) = intercept + genotypes")

    fold_outdir = joinpath(outdir_root, "fold_$(tag)")

    out = runMCMC(
        model,
        ph_train;
        output_folder = fold_outdir,
        chain_length = chain_length,
        burnin = burn_in,
        output_samples_frequency = save_every,
        seed = seed
    )

    # Marker effects and design matrices
    meff = out["marker effects genotypes"]
    markers = Symbol.(meff.Marker_ID)
    Z_test  = Matrix(g_test[!, markers])

    Zt = Matrix(g_test[!, Symbol.(meff.Marker_ID)])

    yhat = (Zt * meff.Estimate) .+ out["location parameters"].Estimate[1]

    acc = cor(yhat, ph_test[!, Symbol(trait)])

    preds_test = insertcols!(select(ph_test, [:ID, Symbol.("$(trait)")]), :EBV => yhat)


    push!(all_acc_np, acc)
    push!(all_h2_np, out["heritability"].Estimate[1])


    CSV.write(joinpath(fold_outdir, "predictions.csv"), preds_test)
    open(joinpath(fold_outdir, "accuracy.txt"), "w") do io
        write(io, string(acc, "\n"))
    end
    

    ID_order = ph_test.ID
    

    y_perm = copy(ph_train[!, Symbol(trait)])
    for i in 1:permute_n
        y_perm = shuffle!(y_perm)
        ph_train_perm = DataFrame("ID" => ph_train.ID, "$(trait)" => y_perm)
        global genotypes = get_genotypes(g_train; method="BayesC", center=false, quality_control=false)
        model_p = build_model("$(trait) = intercept + genotypes")
        out_p = runMCMC(
            model_p,
            ph_train_perm;
            output_folder = joinpath(fold_outdir, "000-permutation_$(i)"),
            chain_length = chain_length,
            burnin = burn_in,
            output_samples_frequency = save_every,
            seed = seed
        )

        # Marker effects
        meff = out_p["marker effects genotypes"]
        markers = Symbol.(meff.Marker_ID)
        Z_test  = Matrix(g_test[!, markers])

        Zt = Matrix(g_test[!, Symbol.(meff.Marker_ID)])

        yhat_p = (Zt * meff.Estimate) .+ out_p["location parameters"].Estimate[1]

        acc = cor(yhat_p, ph_test[!, Symbol(trait)])

        preds_test = insertcols!(select(ph_test, [:ID, Symbol.("$(trait)")]), :EBV => yhat_p)



        acc_perm = cor(yhat_p, ph_test[!, Symbol(trait)])
        push!(all_acc_p, acc_perm)
        push!(all_h2_p, out_p["heritability"].Estimate[1])

    end

    if fold == nfold
        println("CV summary across all folds:")

        all_acc_np = filter(!isnan, all_acc_np)  # remove NaNs for mean/SD calculation
        cv_mean_np = mean(skipmissing(all_acc_np))
        cv_sd_np   = std(skipmissing(all_acc_np))
        acc_n_np   = length(all_acc_np)
        println("JWAS BayesC CV accuracy (mean ± SD): ",
                round(cv_mean_np, digits=4), " ± ", round(cv_sd_np, digits=4))
        open(joinpath(outdir_root, "cv_NP_accuracy_summary.txt"), "w") do io
            write(io, "N=$(acc_n_np)\nmean=$(cv_mean_np)\nsd=$(cv_sd_np)\nvalues=$(join(all_acc_np, ','))\n")
        end

        all_h2_np  = filter(!isnan, all_h2_np)  # remove NaNs for mean/SD calculation
        h2_mean_np = mean(skipmissing(all_h2_np))
        h2_sd_np   = std(skipmissing(all_h2_np))
        h2_n_np    = length(all_h2_np)
        println("JWAS BayesC CV heritability (mean ± SD): ",
                round(h2_mean_np, digits=4), " ± ", round(h2_sd_np, digits=4))
        open(joinpath(outdir_root, "cv_NP_heritability_summary.txt"), "w") do io
            write(io, "N=$(h2_n_np)\nmean=$(h2_mean_np)\nsd=$(h2_sd_np)\nvalues=$(join(all_h2_np, ','))\n")
        end

        all_acc_p = filter(!isnan, all_acc_p)  # remove NaNs for mean/SD calculation
        cv_mean_p = mean(skipmissing(all_acc_p))
        cv_sd_p   = std(skipmissing(all_acc_p))
        acc_n_p   = length(all_acc_p)
        println("JWAS BayesC CV-permuted accuracy (mean ± SD): ",
                round(cv_mean_p, digits=4), " ± ", round(cv_sd_p, digits=4))
        open(joinpath(outdir_root, "cv_permute_accuracy_summary.txt"), "w") do io
            write(io, "N=$(acc_n_p)\nmean=$(cv_mean_p)\nsd=$(cv_sd_p)\nvalues=$(join(all_acc_p, ','))\n")
        end

        all_h2_p  = filter(!isnan, all_h2_p)  # remove NaNs for mean/SD calculation
        h2_mean_p = mean(skipmissing(all_h2_p))
        h2_sd_p   = std(skipmissing(all_h2_p))
        h2_n_p    = length(all_h2_p)
        println("JWAS BayesC CV-permuted heritability (mean ± SD): ",
                round(h2_mean_p, digits=4), " ± ", round(h2_sd_p, digits=4))
        open(joinpath(outdir_root, "cv_permute_heritability_summary.txt"), "w") do io
            write(io, "N=$(h2_n_p)\nmean=$(h2_mean_p)\nsd=$(h2_sd_p)\nvalues=$(join(all_h2_p, ','))\n")
        end


    end
end



test_acc = UnequalVarianceTTest(all_acc_p, all_acc_np)
p_acc = pvalue(test_acc)

histogram(
    all_acc_p, 
    bins=250, 
    normalize = true,
    alpha     = 0.5,
    label = "Null distribution",
    xlabel="Accuracy", 
    ylabel="Frequency", 
    title="$(permute_n) permutations CV (n=$(nfold)) for $(tissue) $(trait)", 
    legend=true,
    color = :lightblue
)
vline!([mean(all_acc_p)], label="Mean of Null", linestyle=:dash, color=:blue)


histogram!(
    all_acc_np;
    bins      = 10,
    normalize = true,
    alpha     = 0.5,
    label     = "CV distribution",
    color = :red
)

vline!([mean(all_acc_np)], label="Mean of CV", linestyle=:dash, color=:red)


annotate!(
    -0.8, 2.9,
    text("p-value = " * @sprintf("%.2e", p_acc), :left, 10)
)

savefig(joinpath(outdir_root, "acc_$(tissue)_$(trait)_$(permute_n)permutations.png"))



test_h2 = UnequalVarianceTTest(all_h2_p, all_h2_np)
p_h2 = pvalue(test_h2)

histogram(
    all_h2_p, 
    bins=250, 
    normalize = true,
    alpha     = 0.5,
    label = "Null distribution",
    xlabel="Heritability", 
    ylabel="Frequency", 
    title="$(permute_n) permutations CV (n=$(nfold)) for $(tissue) $(trait)", 
    legend=true,
    color = :lightblue
)

vline!([mean(all_h2_p)], label="Mean of Null", linestyle=:dash, color=:blue)

histogram!(
    all_h2_np;
    bins      = 10,
    normalize = true,
    alpha     = 0.5,
    label     = "CV distribution",
    color = :red
)

vline!([mean(all_h2_np)], label="Mean of CV", linestyle=:dash, color=:red)


annotate!(
    0.2, 120,
    text("p-value = " * @sprintf("%.2e", p_h2), :left, 10)
)

savefig(joinpath(outdir_root, "h2_$(tissue)_$(trait)_$(permute_n)permutations.png"))



closeall()







# ----------------------------- Run all folds, binary trait ----------------------------------
nfold = discover_nfold(splits_dir)
println("Detected nfold = ", nfold)
all_auc = Float64[]
all_h2 = Float64[]
mkdir(outdir_root)
for fold in 1:nfold
    tag = Int(fold)
    println("=== Fold $fold / $nfold ===")

    ph_train, ph_test, g_train, g_test = read_fold_tables(splits_dir, tag, trait)

    global genotypes = get_genotypes(g_train; method="BayesC", center=false, quality_control=false)

    model_equation = "$(trait) = intercept + genotypes"
    model = build_model(model_equation, categorical_trait=true)

    fold_outdir = joinpath(outdir_root, "fold_$(tag)")

    out = runMCMC(
        model,
        ph_train;
        output_folder = fold_outdir,
        chain_length = chain_length,
        output_samples_frequency = save_every,
        seed = seed
    )

    # Marker effects and design matrices
    meff = out["marker effects genotypes"]
    markers = Symbol.(meff.Marker_ID)
    Z_test  = Matrix(g_test[!, markers])
    Z_train = Matrix(g_train[!, markers])

    # Scores (can be any monotone transform for ROC)
    β̂ = meff.Estimate
    intercrept = out["location parameters"].Estimate[1]
    score_test  = Vector{Float64}((Z_test  * β̂) .+ intercrept)
    score_train = Vector{Float64}((Z_train * β̂) .+ intercrept)

    y_test_raw  = ph_test[!,  Symbol(trait)]
    y_train_raw = ph_train[!, Symbol(trait)]

    # If your column is already 0/1 Ints:
    y_test_bool  = y_test_raw .== 1
    y_train_bool = y_train_raw .== 1


    # --- drop missings if any exist in either vector (defensive) ---
    mask_test  = .!ismissing.(y_test_bool) .& .!ismissing.(score_test)
    mask_train = .!ismissing.(y_train_bool) .& .!ismissing.(score_train)


    y_test_b   = collect(skipmissing(y_test_bool[mask_test]))
    s_test     = score_test[mask_test]

    y_train_b  = collect(skipmissing(y_train_bool[mask_train]))
    s_train    = score_train[mask_train]

    # --- ROC/AUC on TEST (labels must be Bool) ---
    roc_test = roc(y_test_b, s_test)
    auc_test = AUC(roc_test)
    auc_test_3sf = round(auc_test; sigdigits=3)
    push!(all_auc, auc_test)
    push!(all_h2, out["heritability"].Estimate[1])


    p = plot(roc_test, label="$(trait)", title="ROC with AUC = $(auc_test_3sf)")
    savefig(p, joinpath(fold_outdir, "$(tissue)_$(trait)_roc_curve.png"))

    preds_test = insertcols!(select(ph_test, [:ID, Symbol(trait)]), :Prediction => s_test)

    CSV.write(joinpath(fold_outdir, "predictions.csv"), preds_test)
    open(joinpath(fold_outdir, "AUC.txt"), "w") do io
        write(io, string(auc_test, "\n"))
    end
    if fold == nfold
        println("CV summary across all folds:")
        all_acc = filter(!isnan, all_auc)  # remove NaNs for mean/SD calculation
        cv_mean = mean(skipmissing(all_acc))
        cv_sd   = std(skipmissing(all_acc))
        acc_n   = length(all_acc)
        println("JWAS BayesC CV AUC (mean ± SD): ",
                round(cv_mean, digits=4), " ± ", round(cv_sd, digits=4))
        open(joinpath(outdir_root, "cv_AUC_summary.txt"), "w") do io
            write(io, "N=$(acc_n)\nmean=$(cv_mean)\nsd=$(cv_sd)\nvalues=$(join(all_acc, ','))\n")
        end
        all_h2  = filter(!isnan, all_h2)  # remove NaNs for mean/SD calculation
        h2_mean = mean(skipmissing(all_h2))
        h2_sd   = std(skipmissing(all_h2))
        h2_n    = length(all_h2)
        println("JWAS BayesC CV heritability (mean ± SD): ",
                round(h2_mean, digits=4), " ± ", round(h2_sd, digits=4))
        open(joinpath(outdir_root, "cv_heritability_summary.txt"), "w") do io
            write(io, "N=$(h2_n)\nmean=$(h2_mean)\nsd=$(h2_sd)\nvalues=$(join(all_h2, ','))\n")
        end
    end
end






# ----------------------------- Line by line code -------------------------------------
# tag = string(1)
# ph_train, ph_test, g_train, g_test = read_fold_tables(splits_dir, tag, trait)
# genotypes = get_genotypes(g_train; method="BayesC", center=false, quality_control=false)
# model = build_model("$(trait) = intercept + genotypes")
# mkdir(outdir_root)
# fold_outdir = joinpath(outdir_root, "fold_$(tag)")
# out = runMCMC(
#     model,
#     ph_train;
#     output_folder = fold_outdir,
#     chain_length = chain_length,
#     output_samples_frequency = save_every,
#     seed = seed
# )
# meff = out["marker effects genotypes"]
# rename!(meff, Dict(:Estimate=>:u, :Marker_ID=>:marker)) 
# Zt = Matrix(g_test[!, Symbol.(meff.marker)])
# u = meff.u
# yhat_EBV = Zt * u
# acc = cor(yhat_EBV, ph_test[!, Symbol(trait)])
# preds_test = insertcols!(select(ph_test, [:ID, Symbol.("$(trait)")]), :EBV => yhat_EBV)
# h2 = out["heritability"].Estimate[1]
# nfold = 5


# nfold = 1
# all_auc = Float64[]
# all_h2  = Float64[]
# outdir = "/Users/dannox/Desktop/JWAS/irene_data/000-temp_cv_results"
# mkdir(outdir)


# ph_train, ph_test, g_train, g_test = read_fold_tables(splits_dir, nfold, trait)
# genotypes = get_genotypes(g_train; method="BayesC", center=false, quality_control=false)
# fold_outdir = joinpath(outdir, "fold_$(nfold)")
# model = build_model("$(trait) = intercept + genotypes")
# out = runMCMC(
#     model,
#     ph_train;
#     output_folder = fold_outdir,
#     chain_length = chain_length,
#     output_samples_frequency = save_every,
#     seed = seed
# )

# # Marker effects and design matrices
# meff = out["marker effects genotypes"]
# markers = Symbol.(meff.Marker_ID)
# Z_test  = Matrix(g_test[!, markers])
# Z_train = Matrix(g_train[!, markers])

# # Scores (can be any monotone transform for ROC)
# β̂ = meff.Estimate
# score_test  = Vector{Float64}(Z_test  * β̂)
# score_train = Vector{Float64}(Z_train * β̂)

# y_test_raw  = ph_test[!,  Symbol(trait)]
# y_train_raw = ph_train[!, Symbol(trait)]

# # If your column is already 0/1 Ints:
# y_test_bool  = y_test_raw .== 1
# y_train_bool = y_train_raw .== 1


# # --- drop missings if any exist in either vector (defensive) ---
# mask_test  = .!ismissing.(y_test_bool) .& .!ismissing.(score_test)
# mask_train = .!ismissing.(y_train_bool) .& .!ismissing.(score_train)


# y_test_b   = collect(skipmissing(y_test_bool[mask_test]))
# s_test     = score_test[mask_test]

# y_train_b  = collect(skipmissing(y_train_bool[mask_train]))
# s_train    = score_train[mask_train]

# predictions = insertcols!(select(ph_test, [:ID, Symbol(trait)]), :Score => s_test)


# # --- ROC/AUC on TEST (labels must be Bool) ---
# r_test = roc(y_test_b, s_test)
# auc_test = AUC(r_test)
# roc_train = roc(y_train_b, s_train)
# auc_train = AUC(roc_train)
# push!(all_auc, auc_val)
# push!(all_h2, out["heritability"].Estimate[1])



# function noisy(label; λ=0.0)
#     if label
#         return 1 - λ*rand()
#     else
#         return λ*rand()
#     end
# end

# labels = rand(Bool, 200);

# scores(λ) = map(labels) do label
#         noisy(label, λ=λ)
#     end

# test = scores(0.6)

# roc_good = roc(test, labels);

# area_good = AUC(roc_good)
# # 3 sigfigs for area_good
# area_good_3sf = round(area_good; sigdigits=3)
# lbl = "AUC = $(area_good_3sf)"

# # Show plot
# display(plot(r_test, label="$(trait)", title="ROC AUC = $(auc_test)"))

