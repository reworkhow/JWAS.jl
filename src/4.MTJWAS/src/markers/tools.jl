function get_column(X,j)
    nrow,ncol = size(X)
    if j>ncol||j<0
        error("column number is wrong!")
    end
    indx = 1 + (j-1)*nrow
    ptr = pointer(X,indx)
    pointer_to_array(ptr,nrow)
end

function get_column_ref(X)
    ncol = size(X)[2]
    xArray = Array(Array{Float64,1},ncol)
    for i=1:ncol
        xArray[i] = get_column(X,i)
    end
    return xArray
end

function center!(X)
    nrow,ncol = size(X)
    colMeans = mean(X,1)
    BLAS.axpy!(-1,ones(nrow)*colMeans,X)
    return colMeans
end

function getXpRinvX(X, Rinv)
    ncol = size(X)[2]
    XpRinvX = [((X[:,i].*Rinv)'X[:,i])[1]::Float64 for i=1:ncol]
    return XpRinvX
end

function getXpRinvX(X)
    XpRinvX = [dot(X[:,i],X[:,i]) for i=1:size(X,2)]
    return XpRinvX
end

type GibbsMats
    X::Array{Float64,2}
    nrows::Int64
    ncols::Int64
    xArray::Array{Array{Float64,1},1}
    xpx::Array{Float64,1}
    function GibbsMats(X::Array{Float64,2}) ###More
        nrows,ncols = size(X)
        xArray = get_column_ref(X)
        XpX = getXpRinvX(X)
        new(X,nrows,ncols,xArray,XpX)
    end
end

function make_valpha(G,pi;option="theory") #from marker effect variance with π=0
    if option=="theory"
        offdiag(a,b)=a*b
    elseif option=="max_qtl" #guanrantee posdef
        offdiag(a,b)=max(a,b)
    elseif option=="min_qtl"
        offdiag(a,b)=min(a,b)
    elseif option=="sqrt"  #guarante posdef
        offdiag(a,b)=sqrt(a*b)
        elseif option=="zero" #make off-dignal of maker vaiance matrix zero
        offdiag(a,b)= Inf
    else
        error("unknown option")
    end

    pi_c        = 1 - pi
    mat_pi_c    = diagm(vec(pi_c))
    for i=1:size(mat_pi_c,2)
        for j=1:size(mat_pi_c,2)
            if i!=j
                mat_pi_c[i,j]=offdiag(mat_pi_c[i,i],mat_pi_c[j,j])
            end
        end
    end
    Valpha    = G./mat_pi_c
    println("marker effect variance with π=$pi is positive definite??: ", isposdef(Valpha))
    Valpha
end
