
#Devloped by Zigui Wang (zigwang97@gmail.com) and Hao Cheng (qtlcheng@ucdavis.edu)
#Code based on
#Gianola, D., & Sorensen, D. (2004). Quantitative Genetic Models for Describing
#Simultaneous and Recursive Relationships Between Phenotypes. Genetics, 167(3),
#1407–1424. http://doi.org/10.1534/genetics.103.025734

# for one individual, with fully simultaneous model (this not fully recursive model)
# Λy =
# [y1 - λ_12⋅y2 - λ_13⋅y3
#  y2 - λ_21⋅y1 - λ_23⋅y3
#  y3 - λ_31⋅y1 - λ_32⋅y2]
# = y - Yλ
# =
# [y1    [y2 y3 0  0  0  0
#  y2  -  0  0  y1 y3 0  0   * [λ_12;λ_13;λ_21;λ_23;λ_31;λ_32]
#  y3]    0  0  0  0  y1 y2]

# Example Y (3 traits, 4 inds) from function `get_sparse_Y_FSM`
#  yi2 yi3    0   0    0   0
#  yj2 yj3    0   0    0   0
#  yk2 yk3    0   0    0   0
#  yl2 yl3    0   0    0   0
#    0   0  yi1 yi3    0   0
#    0   0  yj1 yj3    0   0
#    0   0  yk1 yk3    0   0
#    0   0  yl1 yl3    0   0
#    0   0    0   0  yi1 yi2
#    0   0    0   0  yj1 yj2
#    0   0    0   0  yk1 yk2
#    0   0    0   0  yl1 yl2

#  yi2 yi3    0   0    0   0
#    0   0  yi1 yi3    0   0
#    0   0    0   0  yi1 yi2
#  yj2 yj3    0   0    0   0
#    0   0  yj1 yj3    0   0
#    0   0    0   0  yj1 yj2
#  yk2 yk3    0   0    0   0
#    0   0  yk1 yk3    0   0
#    0   0    0   0  yk1 yk2
#  yl2 yl3    0   0    0   0
#    0   0  yl1 yl3    0   0
#    0   0    0   0  yl1 yl2


#1)starting values for structural coefficient λij (i≂̸j) is zero,
#2)the variable "ycorr" is used for "Λycorr" for coding convinience,
#and starting values for ycorr (i.e., Λycorr) is the ycorr obtained above
#3)no missing phenotypes in SEM
#4)Λ is (I-Λ) in Wang et al. 2020 BioRxiv.
function SEM_setup(wArray,causal_structure,mme)
    ntraits = length(wArray)
    nobs    = length(wArray[1])
    Y  = get_sparse_Y_FRM(wArray,causal_structure) #here wArray is for phenotypes (before corrected)
    Λ  = Matrix{Float64}(I,ntraits,ntraits) #structural coefficient λij (i≂̸j) is zero (starting values)
    Λy = kron(Λ,sparse(1.0I,nobs,nobs))*mme.ySparse
    causal_structure_filename = "structure_coefficient_MCMC_samples.txt"
    causal_structure_outfile  = open(causal_structure_filename,"w")   #write MCMC samples for Λ to a txt file



    return Y,Λy,causal_structure_outfile
end
# Get Y for all individuals ordered as individuals within traits (fully simultaneous model)
function get_sparse_Y_FSM(wArray)
    ntraits = length(wArray)
    data = wArray[1]
    for i in 2:length(wArray)
        data = [data wArray[i]]
    end

    Y = kron(sparse(1.0I, ntraits, ntraits),sparse(data))

    keep_col = Array{Bool}(undef,0)

    for i in 1:ntraits
        for j in 1:ntraits
            if i == j
               push!(keep_col,false)
            else
               push!(keep_col,true)
            end
        end

    end
    Y = Y[:,keep_col]  #size: #row: ntraits*nind #col: ntraits*(ntraits-1)
    return Y
end

# Get Y for all individuals ordered as individuals within traits (fully recursive model)
function get_sparse_Y_FRM(wArray,causal_structure)
    Y        = get_sparse_Y_FSM(wArray)
    ntraits  = length(wArray)
    keep_col = Array{Bool}(undef,0)
    for i in 1:ntraits
        for j in 1:ntraits
            if i != j
                if  causal_structure[i,j]==1
                    push!(keep_col,true)
                else
                    push!(keep_col,false)
                end
            end
        end
    end
    Y = Y[:,keep_col]
    return Y
end

#sample Λ and update Λycorr
function get_Λ(Y,R,Λycorr,Λy,y,causal_structure)

    ntraits = size(R,1)
    nind    = div(length(Λycorr),ntraits)
    #Define residual matrix
    bigR = kron(sparse(1:ntraits,1:ntraits,diag(R).^(-1)),sparse(1.0I,nind,nind))

    #formula calculation
    Y_R_product = Y'bigR

    #λ~N(1μ,Iσ2)
    λ_σ2     = 1.0
    λ_μ      = 0.0
    first    = Y_R_product*Y
    first   += 1/λ_σ2*sparse(I,size(first))

    #ycor   is y  - Xβ - Zu
    #Λycorr is Λy - Xβ - Zu
    #Λycorr is used in smaplings for all other parameters. However,
    #ycor, instead of Λycorr,is used for sampling Λ.
    ycorr  = Λycorr - Λy + y  #subtract old Λy
    second = Y_R_product*ycorr
    second += ones(length(second))*(λ_μ/λ_σ2)

    first_inv = inv(Matrix(first))
    mu        = vec(first_inv*second)
    var       = Symmetric(first_inv)

    # λ =  [λ12;λ21]
    # Λ =  [1 -λ12
    #      -λ21 1]
    λ = rand(MvNormal(mu,var))

    causal_matrix = tranform_lambda(λ,causal_structure)
    Λ             = I - causal_matrix
    λ_vec         = vec(causal_matrix)

    Λy[:]      = kron(Λ,sparse(1.0I,nind,nind))*y
    Λycorr[:]  = ycorr - y + Λy #add new Λy
    return λ,λ_vec
end

function tranform_lambda(lambda,causal_structure)
   causal_structure = sparse(causal_structure)
   row_index = findnz(causal_structure)[1]
   col_index = findnz(causal_structure)[2]
   Lambda    = zeros(size(causal_structure))
   for j in 1:length(lambda)
       Lambda[row_index[j],col_index[j]] = lambda[j]
   end
   return Lambda
end


# generate the MCMC samples for indirect marker effect
# In each iteration, we want to form the following equation
#
# [0    0     0         [marker1_y1 marker2_y1 ...       [0    0
#  λ12    0     0   *    marker1_y2 marker2_y2 ...  =     λ12*marker1_y1 λ12*marker2_y1
#  λ13    0     0]       marker1_y3 marker2_y3 ...]       λ13*marker1_y1 λ12*marker2_y1 ]
# = Λ * [α1,α2,...] = indirect effect

function generate_indirect_marker_effect_sample(phenotypes,output_folder,causal_structure,structure_coefficient_path)


    number_traits = size(causal_structure,1) # the row number of causal structure matrix is number of traits
    trait_vec     = names(phenotypes)[2:number_traits+1] # vector of trait name, starting from 2 for skipping the "ID"

    λ_file        = CSV.read(structure_coefficient_path, DataFrame,header = false)

    direct_effect_sample = Dict()
    io_diction  = Dict()
    # read the direct effect sample file and create indirect sample file for each trait
    for i in 1:number_traits
        direct_file_name                   = output_folder*"/MCMC_samples_marker_effects_genotypes_"*trait_vec[i]*".txt"
        direct_effect_sample[trait_vec[i]] = CSV.read(direct_file_name,DataFrame,header = true)

        indirect_file_name = output_folder*"/MCMC_samples_indirect_marker_effects_genotypes_"*trait_vec[i]*".txt"
        io_diction[trait_vec[i]] = indirect_file_name
    end

    number_sample = size(direct_effect_sample[trait_vec[1]],1)
    number_marker = size(direct_effect_sample[trait_vec[1]],2)

    marker_header = permutedims(names(direct_effect_sample[trait_vec[1]]))

    # write marker header for each indirect effect file
    for i in 1:number_traits
        open(io_diction[trait_vec[i]], "a") do io
                   writedlm(io, marker_header ,", ")
               end;

    end


    number_sample = 1

    # compute the indirect effect in each sample and write to the target file
    for i in 1:number_sample
        direct_effect_current_sample = zeros(number_traits,number_marker)
        λ_sample                     = Vector(λ_file[i,1:end])
        for j in 1:number_traits
            current_effect = direct_effect_sample[trait_vec[j]][i,1:end]
            current_effect = Vector(current_effect)

            direct_effect_current_sample[j,1:end] = current_effect'

        end

        causal_matrix = reshape(λ_sample ,number_traits,number_traits)
        indirect_effect_current_sample = compute_indirect_effect(causal_matrix,direct_effect_current_sample )

        for i in 1:number_traits
            indirect_effect_this_trait = indirect_effect_current_sample[i,1:end]
            open(io_diction[trait_vec[i]], "a") do io
                       writedlm(io, indirect_effect_this_trait' ,", ")
                   end;

        end
    end


end

# compute the indirect marker effect based on the formula in the SEM paper
function compute_indirect_effect(Λ,marker_effects )
    number_traits,number_marker = size(marker_effects)
    result                      = zeros(number_traits,number_marker)
    for i in 1:number_traits-1
        result += Λ^i * marker_effects
    end
    return (result)
end
