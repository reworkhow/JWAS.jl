## Data encryption with random effects


# Goal: obtain the number of levels of each random effect, then generate a fake column in df.
#       This step is to make sure JWAS can generate a (fake) incidence matrix with correct size.
# Note: the user must name the encrypted indidence matrix as "name_i", e.g., "batch_1", "batch_2" for "batch" with 2 levels.
function add_encrypted_random_variable_column!(mme,df)
      name_all = names(df)
      for rndTrm in mme.rndTrmVec
            rnd_variable_name = split(rndTrm.term_array[1],":")[2] #e.g., "batch"
            nlevels,i = 0,1
            while rnd_variable_name*"_"*string(i) ∈ name_all
                  i+=1
                  nlevels+=1
            end
            df[:,rnd_variable_name]         = rand(1:nlevels,size(df,1))
            df[1:nlevels,rnd_variable_name] = collect(1:nlevels) #make sure all levels are generated
            if length(unique(df[:,rnd_variable_name]))==nlevels
                  printstyled("Encryption: a random column has been generated for $rnd_variable_name with $nlevels levels. \n",bold=false,color=:green)
            else
                  error("wrong number of levels!")
            end
      end
end

# Goal: replace the fake incidence matrix with user's encrypted incidence matrix.
# Note: the user must name the encrypted indidence matrix as "name_i", e.g., "batch_1", "batch_2" for "batch" with 2 levels.
function use_encrypted_X!(mme,df)
      for modelTerm in mme.modelTerms
            if modelTerm.random_type == "I"
                  nlevel=modelTerm.nLevels
                  rnd_variable_name = split(modelTerm.trmStr,":")[2] #e.g., "batch"
                  cols_to_select    = [rnd_variable_name * "_" * string(i) for i in 1:modelTerm.nLevels] #e.g., "batch_1","batch_2",...
                  modelTerm.X       = sparse(Matrix(select(df, cols_to_select)))
                  printstyled("Encryption: the enctypted incidece matrix for $rnd_variable_name is obtained from input data ($nlevel levels). \n",bold=false,color=:green)
            end
      end
end

