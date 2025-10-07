using HTTP;
using CSV, DataFrames, JWAS, JWAS.Datasets, Statistics;
using Plots;

function dataset_big(file_name::AbstractString; dataset_url::AbstractString="https://raw.githubusercontent.com/zhaotianjing/bio_protocol/main/data")
    rdaname = joinpath(dataset_url, string(file_name))
    http_obj = HTTP.get(rdaname)
    return http_obj
end

phenofile = dataset_big("phenotypes.txt")
pedfile = dataset_big("pedigree.txt")
genofile = dataset_big("genotypes_5k.txt")
genofile_ssbr = dataset_big("genotypes_2kind.txt")

phenotypes = CSV.read(phenofile.body, DataFrame, delim=',', header=true, missingstrings=["."])
phenotypes[!, :dam] = phenotypes[!, :ID]
pedfile = CSV.read(pedfile.body, DataFrame, delim=',', header=true, missingstrings=["."])
pedigree = get_pedigree(pedfile, separator=",", header=true);
genofile = CSV.read(genofile.body, DataFrame, delim=',', header=true, missingstrings=["."])
genofile_ssbr = CSV.read(genofile_ssbr.body, DataFrame, delim=',', header=true, missingstrings=["."])

if ispath("mytest") == true
    rm("mytest", recursive=true)
end
mkdir("mytest/")
cd("mytest/")

PA_dict = Dict{String,Float64}()
for single_step in [false, true]
    for test_method in ["BayesA", "BayesB", "BayesC", "RR-BLUP", "BayesL", "GBLUP", "PBLUP"]
        # ST analysis
        newdir = "ST_" * (single_step ? "SS" : "") * test_method * "/"
        mkdir(newdir)
        cd(newdir)
        if test_method in ["BayesC", "BayesB"]
            test_estimatePi = true
        else
            test_estimatePi = false
        end

        printstyled("\n\n\n\n\n\n\n\nTest single-trait $test_method analysis using $(single_step ? "in" : "")complete genomic data\n\n\n", bold=true, color=:green)
        if test_method != "PBLUP"
            G3 = 1.0
            if single_step
                global geno = get_genotypes(genofile_ssbr, G3, header=true, separator=',', method=test_method, estimatePi=test_estimatePi)
            else
                global geno = get_genotypes(genofile, G3, header=true, separator=',', method=test_method, estimatePi=test_estimatePi)
            end
        end

        if test_method != "PBLUP"
            model_equation1 = "t1 = intercept + ID + dam + geno"
        else
            model_equation1 = "t1 = intercept + ID + dam"
        end

        R = 1.0
        model1 = build_model(model_equation1, R)

        G2 = [1.0 0.5
            0.5 1.0]
        set_random(model1, "ID dam", pedigree, G2)


        if single_step == false # genomic models or PBLUP
            out1 = runMCMC(model1, phenotypes, heterogeneous_residuals=false, #estimate all variances==true by default
                double_precision=true,
                chain_length=1000, output_samples_frequency=10,
                printout_frequency=100, seed=314)
            results = innerjoin(out1["EBV_t1"], phenotypes, on=:ID)
            ind_id = findall(x -> !ismissing(x), results[!, :t1])
            accuracy = cor(results[ind_id, :EBV], results[ind_id, :t1])
            PA_dict[newdir] = accuracy
        elseif single_step == true && test_method != "PBLUP" && test_method != "GBLUP"
            # genomic model, no PBLUP/GBLUP
            out1 = runMCMC(model1, phenotypes, heterogeneous_residuals=false, #estimate all variances==true by default
                chain_length=1000, output_samples_frequency=10, printout_frequency=100,
                single_step_analysis=true, pedigree=pedigree, seed=314)
            results = innerjoin(out1["EBV_t1"], phenotypes, on=:ID)
            ind_id = findall(x -> !ismissing(x), results[!, :t1])
            accuracy = cor(results[ind_id, :EBV], results[ind_id, :t1])
            PA_dict[newdir] = accuracy
        end
        # GWAS
        if test_method != "PBLUP" && test_method != "GBLUP"
            gwas1 = GWAS("results/MCMC_samples_marker_effects_geno_t1.txt")
            show(gwas1)
            println()
            #gwas2=GWAS("results/MCMC_samples_marker_effects_geno_y1.txt",mapfile,model1)
            #show(gwas2)
        end
        cd("..")

        printstyled("\n\n\n\n\n\n\n\nTest multi-trait $test_method analysis using $(single_step ? "in" : "")complete genomic data\n\n\n", bold=true, color=:green)
        # MT analysis
        newdir = "MT_" * (single_step ? "SS" : "") * test_method * "/"
        mkdir(newdir)
        cd(newdir)

        if test_method != "PBLUP"
            G3 = [1.0 0.5 0.5
                0.5 1.0 0.5
                0.5 0.5 1.0]
            if single_step
                global geno = get_genotypes(genofile_ssbr, G3, header=true, separator=',', method=test_method, estimatePi=test_estimatePi)
            else
                global geno = get_genotypes(genofile, G3, header=true, separator=',', method=test_method, estimatePi=test_estimatePi)
            end
        end

        if test_method != "PBLUP"
            model_equation2 = "t1 = intercept + ID + dam + geno
                               t2 = intercept + ID + geno
                               t3 = intercept + ID + geno"
        else
            model_equation2 = "t1 = intercept + ID + dam
                               t2 = intercept + ID
                               t3 = intercept + ID"
        end


        R = [1.0 0.5 0.5
            0.5 1.0 0.5
            0.5 0.5 1.0]
        model2 = build_model(model_equation2, R)

        G2 = [1.0 0.5 0.5 0.5
            0.5 1.0 0.5 0.5
            0.5 0.5 1.0 0.5
            0.5 0.5 0.5 1.0]
        set_random(model2, "ID dam", pedigree, G2)


        if single_step == false
            out2 = runMCMC(model2, phenotypes, heterogeneous_residuals=false, double_precision=true, #estimate all variances==true by default
                chain_length=1000, output_samples_frequency=10, printout_frequency=100, seed=314)
            results = innerjoin(out2["EBV_t1"], phenotypes, on=:ID)
            ind_id = findall(x -> !ismissing(x), results[!, :t1])
            accuracy = cor(results[ind_id, :EBV], results[ind_id, :t1])
            PA_dict[newdir] = accuracy
        elseif single_step == true && test_method != "PBLUP" && test_method != "GBLUP"
            out2 = runMCMC(model2, phenotypes, heterogeneous_residuals=false, double_precision=true, #estimate all variances==true by default
                chain_length=1000, output_samples_frequency=10, printout_frequency=100,
                single_step_analysis=true, pedigree=pedigree, seed=314)
            results = innerjoin(out2["EBV_t1"], phenotypes, on=:ID)
            ind_id = findall(x -> !ismissing(x), results[!, :t1])
            accuracy = cor(results[ind_id, :EBV], results[ind_id, :t1])
            PA_dict[newdir] = accuracy
        end
        if test_method != "PBLUP" && test_method != "GBLUP"
            gwas1 = GWAS("results/MCMC_samples_marker_effects_geno_t1.txt")
            show(gwas1)
            println()
            #gwas2=GWAS(model2,mapfile,"results/MCMC_samples_marker_effects_geno_y1.txt","results/MCMC_samples_marker_effects_geno_y2.txt",genetic_correlation=true)
            #show(gwas2)
        end
        cd("..")
    end
end

PA_dict

#PA_df_RRMnew = DataFrame(PredictionAccuracy=collect(values(PA_dict)), Methods=collect(keys(PA_dict)))
#PA_df_RRMnew[!, :Methods] .= replace.(PA_df_RRMnew[!, :Methods], "/" => "")
#PA_df_RRMnew
#PA_df_normal = CSV.read("PA_bigdata.csv", DataFrame)
#df_joined = innerjoin(PA_df_normal, PA_df_RRMnew, on=:Methods, makeunique=true)
#histogram(abs.(df_joined[!, :PredictionAccuracy] .- df_joined[!, :PredictionAccuracy_1]))
#maximum(abs.(df_joined[!, :PredictionAccuracy] .- df_joined[!, :PredictionAccuracy_1]))