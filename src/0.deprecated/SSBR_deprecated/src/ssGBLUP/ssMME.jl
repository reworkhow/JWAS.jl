function ssMME(all_M,all_y,all_J,all_Z,all_X,all_W,all_A,all_num,vRes,vG)

    vAlpha = vG/all_num.num_markers
    λ1 = vRes/vG
    λ2 = vRes/vAlpha
    num_markers = all_num.num_markers

    X  =all_X.X
    M  =all_M.M
    J  =all_J.J
    W  =all_W.W
    X1 =all_X.X1
    Z1 =all_Z.Z1
    W1 =all_W.W1
    y  =all_y.y
    y1 =all_y.y1
    Ai11=all_A.Ai11
    num_nongeno=all_num.num_g1


    dimX = size(X,2)#add small numbers to diagonal of X'X in case X'X is singular
    ###solve MME
    lhs = [hcat(X'X+eye(dimX)*0.001,    X'W,                     X1'Z1);
           hcat(W'X,    W'W+eye(num_markers)λ2,  W1'Z1);
           hcat(Z1'X1,  Z1'W1,                  Z1'Z1+Ai11*λ1 )]

    rhs = [X'y; W'y; Z1'y1]

    sol =lhs\rhs

    #get EBVs
    end_fixed     =size(X,2)
    start_marker  =end_fixed+1
    end_marker    =start_marker+num_markers-1
    start_epsi    =end_marker+1
    end_epsi      =start_epsi+num_nongeno-1

    beta_hat = sol[1:end_fixed]
    mu_g = sol[end_fixed]
    alpha_hat=sol[start_marker:end_marker]
    epsi_hat = sol[start_epsi:end_epsi]

    mu = sol[1]
    aHat = mu+J*mu_g+M*alpha_hat
    aHat[1:num_nongeno,:] += epsi_hat
    return aHat,alpha_hat,beta_hat,epsi_hat
end



# df = readtable("bv.txt", eltypes =[String, Float64], separator = ' ',header=false)
# a  = Array(Float64,num_ped)
# for (i,ID) in enumerate(df[:,1])
#      j = ped.idMap[ID].seqID
#      a[j] = df[i,2]
# end
# cor(a,aHat)

# #keep (overall correlation)
# #70 0.854
# #50 0.828
# #30 0.822
# #10 0.819
# #5  0.817
# #2  0.799

# #PBLUP
function PBLUP(all_y,all_Z,all_A,all_num,vRes,vG)
    num_obs=length(all_y.y)
    pX = ones(num_obs)
    λ1 = vRes/vG
    y  = all_y.y
    Z  = all_Z.Z
    Ai = all_A.Ai

    lhs = [hcat(pX'pX,    pX'Z );
            hcat(Z'pX,    Z'Z+Ai*λ1)]
    rhs = [pX'y; Z'y]
    sol =lhs\rhs
    aHat = sol[2:length(sol)]
    return aHat
end

# cor(a,aHat)
# #PBLUP 0.813
