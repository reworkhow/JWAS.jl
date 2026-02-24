# Unit tests for solver algorithms (Jacobi, Gauss-Seidel, Gibbs)
using Test, JWAS, DataFrames, CSV, JWAS.Datasets
using LinearAlgebra, SparseArrays

# Test the internal solver algorithms directly with a small linear system
# Ax = b where A is positive definite
A = Float32[4.0 1.0; 1.0 3.0]
b = Float32[1.0; 2.0]
x_true = A \ b  # exact solution

@testset "Jacobi solver" begin
    x = zeros(Float32, 2)
    x_sol = JWAS.Jacobi(A, x, b, tolerance=1e-6, maxiter=1000, printout_frequency=10000)
    @test all(isapprox.(x_sol, x_true, atol=1e-3))
end

@testset "Gauss-Seidel solver" begin
    x = zeros(Float32, 2)
    x_sol = JWAS.GaussSeidel(A, x, b, tolerance=1e-6, maxiter=1000, printout_frequency=10000)
    @test all(isapprox.(x_sol, x_true, atol=1e-3))
end

@testset "Gibbs sampler (single-trait lambda)" begin
    x = zeros(Float64, 2)
    A64 = Float64.(A)
    b64 = Float64.(b)
    x_sol = JWAS.Gibbs(A64, x, b64, 1.0, 5000, printout_frequency=100000)
    # Gibbs is stochastic, so use looser tolerance
    @test all(isapprox.(x_sol, Float64.(x_true), atol=0.5))
end

@testset "Gibbs sampler (multi-trait general)" begin
    x = zeros(Float64, 2)
    A64 = Float64.(A)
    b64 = Float64.(b)
    x_sol = JWAS.Gibbs(A64, x, b64, 5000, printout_frequency=100000)
    @test all(isapprox.(x_sol, Float64.(x_true), atol=0.5))
end

@testset "Solvers converge to same solution" begin
    x1 = zeros(Float32, 2)
    sol_jacobi = JWAS.Jacobi(A, x1, b, tolerance=1e-6, maxiter=5000, printout_frequency=10000)

    x2 = zeros(Float32, 2)
    sol_gs = JWAS.GaussSeidel(A, x2, b, tolerance=1e-6, maxiter=5000, printout_frequency=10000)

    @test all(isapprox.(sol_jacobi, sol_gs, atol=1e-3))
end

@testset "Larger system" begin
    n = 10
    M = randn(Float32, n, n)
    A_large = M' * M + Float32(n) * I  # positive definite
    b_large = randn(Float32, n)
    x_exact = A_large \ b_large

    x = zeros(Float32, n)
    sol = JWAS.Jacobi(A_large, x, b_large, tolerance=1e-6, maxiter=10000, printout_frequency=100000)
    @test all(isapprox.(sol, x_exact, atol=1e-2))

    x = zeros(Float32, n)
    sol = JWAS.GaussSeidel(A_large, x, b_large, tolerance=1e-6, maxiter=10000, printout_frequency=100000)
    @test all(isapprox.(sol, x_exact, atol=1e-2))
end
