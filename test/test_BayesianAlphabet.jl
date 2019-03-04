
using CSV,DataFrames,JWAS,JWAS.Datasets
phenofile       = Datasets.dataset("example","phenotypes.txt")
phenofile_ssbr  = Datasets.dataset("example","phenotypes_ssbr.txt")
pedfile    = Datasets.dataset("example","pedigree.txt")
genofile   = Datasets.dataset("example","genotypes.txt")
phenotypes = CSV.read(phenofile,delim = ',',header=true)
phenotypes_ssbr = CSV.read(phenofile_ssbr,delim = ',',header=true)
pedigree   = get_pedigree(pedfile,separator=",",header=true);

for single_step in [false,true]
      for test_method in ["BayesC","BayesB","RR-BLUP","GBLUP","non_genomic"]
            if test_method in ["BayesC","BayesB"]
                  test_estimatePi = true
            else
                  test_estimatePi = false
            end

            printstyled("\n\n\n\n\n\n\n\nTest single-trait Bayesian Alphabet analysis using complete genomic data\n\n\n",bold=true,color=:red)

            model_equation1  ="y1 = intercept + x1*x3 + x2 + x3 + ID + dam";
            R      = 1.0
            model1 = build_model(model_equation1,R);

            set_covariate(model1,"x1");
            G1 = 1.0
            G2 = [1.0 0.5
                  0.5 1.0]
            set_random(model1,"x2",G1);
            set_random(model1,"ID dam",pedigree,G2);

            if test_method != "non_genomic"
                  G3 =1.0
                  add_genotypes(model1,genofile,G3,header=true,separator=',');
            end
            outputMCMCsamples(model1,"x2")

            if single_step == false && test_method!="non_genomic"
                  out1=runMCMC(model1,phenotypes,methods=test_method,estimatePi=test_estimatePi,chain_length=100,output_samples_frequency=10,printout_frequency=50);
            elseif single_step == true && test_method!="non_genomic"
                  out1=runMCMC(model1,phenotypes_ssbr,methods=test_method,estimatePi=test_estimatePi,chain_length=100,output_samples_frequency=10,printout_frequency=50,
                              single_step_analysis=true,pedigree=pedigree);
            elseif test_method=="non_genomic"
                  out1=runMCMC(model1,phenotypes,chain_length=100,output_samples_frequency=10,printout_frequency=50);
            end

            printstyled("\n\n\n\n\n\n\n\nTest multi-trait Bayesian Alphabet analysis using complete genomic data\n\n\n",bold=true,color=:red)

            model_equation2 ="y1 = intercept + x1 + x3 + ID + dam
                              y2 = intercept + x1 + x2 + x3 + ID
                              y3 = intercept + x1 + x1*x3 + x2 + ID";

            R      = [1.0 0.5 0.5
                      0.5 1.0 0.5
                      0.5 0.5 1.0]
            model2 = build_model(model_equation2,R);

            set_covariate(model2,"x1");
            G1 = [1.0 0.5
                  0.5 1.0]
            G2 = [1.0 0.5 0.5 0.5
                  0.5 1.0 0.5 0.5
                  0.5 0.5 1.0 0.5
                  0.5 0.5 0.5 1.0]
            set_random(model2,"x2",G1);
            set_random(model2,"ID dam",pedigree,G2);

            if test_method != "non_genomic"
                  G3 = [1.0 0.5 0.5
                        0.5 1.0 0.5
                        0.5 0.5 1.0]
                  add_genotypes(model2,genofile,G3,separator=',');
            end
            outputMCMCsamples(model2,"x2")



            if single_step == false && test_method!="non_genomic" && test_method!="GBLUP"
                  out2=runMCMC(model2,phenotypes,methods=test_method,estimatePi=test_estimatePi,chain_length=100,output_samples_frequency=10,printout_frequency=50);
            elseif single_step == true && test_method!="non_genomic" && test_method!="GBLUP"
                  out2=runMCMC(model2,phenotypes_ssbr,methods=test_method,estimatePi=test_estimatePi,chain_length=100,output_samples_frequency=10,printout_frequency=50,
                              single_step_analysis=true,pedigree=pedigree);
            elseif test_method=="non_genomic"
                  out2=runMCMC(model2,phenotypes,chain_length=100,output_samples_frequency=10,printout_frequency=50);
            end
      end
end
