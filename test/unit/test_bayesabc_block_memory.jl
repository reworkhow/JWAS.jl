using Test
using JWAS
using Random
using LinearAlgebra

function make_block_inputs(nobs::Int,nmarkers::Int,block_size::Int)
    X = randn(Float64,nobs,nmarkers)
    Rinv = ones(Float64,nobs)
    starts = collect(1:block_size:nmarkers)
    XArray = [view(X,:,s:min(s+block_size-1,nmarkers)) for s in starts]
    XRinvArray = [Xb' * Diagonal(Rinv) for Xb in XArray]
    XpRinvX = [XRinvArray[i] * XArray[i] for i in eachindex(XArray)]
    xpRinvx = [dot(view(X,:,j),view(X,:,j)) for j in 1:nmarkers]
    yCorr = randn(Float64,nobs)
    α = randn(Float64,nmarkers)
    β = copy(α)
    δ = ones(Float64,nmarkers)
    return (XArray,XRinvArray,xpRinvx,X,XpRinvX,yCorr,α,β,δ,1.0,fill(1.0,nmarkers),0.95)
end

function measure_block_alloc(nobs::Int,nmarkers::Int,block_size::Int)
    args = make_block_inputs(nobs,nmarkers,block_size)
    Random.seed!(123)
    JWAS.BayesABC_block!(args...)
    args = make_block_inputs(nobs,nmarkers,block_size)
    Random.seed!(123)
    return @allocated JWAS.BayesABC_block!(args...)
end

@testset "BayesABC_block allocation scaling" begin
    nobs = 128
    block_size = 8
    alloc_small = measure_block_alloc(nobs,64,block_size)
    alloc_large = measure_block_alloc(nobs,512,block_size)
    @test alloc_large <= alloc_small * 12
end
