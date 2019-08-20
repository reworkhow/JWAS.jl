
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
    Λ = I - tranform_lambda(λ,causal_structure)


    Λy[:]      = kron(Λ,sparse(1.0I,nind,nind))*y
    Λycorr[:]  = ycorr - y + Λy #add new Λy
    return λ
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
