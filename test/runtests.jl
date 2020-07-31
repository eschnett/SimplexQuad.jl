using SimplexQuad

using Test

function integrate(f, X::AbstractMatrix, W::AbstractVector)
    np = size(X, 1)
    @assert size(W, 1) == np
    return sum(W[i] * f(X[i, :]) for i in 1:np)
end

const Types = [Float32, Float64]
@testset "Integrate polynomials T=$T D=$D N=$N" for T in Types,
D in 1:5,
N in 1:(11 - D)

    X, W = simplexquad(T, N, D)

    # Check returned types
    X::AbstractMatrix{T}
    W::AbstractVector{T}
    # Check returned sizes
    np = N^D
    @test size(X, 1) == np
    @test size(X, 2) == D
    @test size(W, 1) == size(X, 1)
    # Ensure all points lie strictly inside the simplex
    function inside(x)
        all(>(0), x) || return false
        all(<(1), x) || return false
        λ = 1 - sum(x)
        return 0 < λ < 1
    end
    @test all(inside(X[i, :]) for i in 1:np)
    # Ensure weights are strictly positive
    @test all(>(0), W)

    # Expect polynomials of orders P<= to be integrated exactly
    P = N
    pmin = CartesianIndex(ntuple(d -> 0, D))
    pmax = CartesianIndex(ntuple(d -> P, D))
    for p in pmin:pmax
        sum(p.I) <= P || continue
        f(x) = prod(x[i]^p[i] for i in 1:D)
        computed = integrate(f, X, W)
        computed::T
        # Mathematica suggests this neat result
        expected = prod(T.(factorial.(p.I))) / factorial(D + sum(p.I))
        expected::T
        @test abs(computed - expected) <= 10 * eps(T)
    end
end
