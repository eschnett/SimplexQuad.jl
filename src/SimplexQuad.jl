module SimplexQuad

using LinearAlgebra



export simplexquad

"""
    simplexquad([::Type{T},] N::Int, n::Int)
    simplexquad(N::Int, vert::AbstractMatrix)

Construct Gauss points and weights for an `n`-dimensional simplex
domain with vertices specified by the `n*(n-1)` matrix `vert`, where
each row contains the coordinates `(x1,...,xn)` for a vertex. The
order of the quadrature scheme in each dimension must be the same in
this implementation.

    N: Number of quadrature points
    T: Return type (defaults to Float64)
    n: Number of dimensions
    vert: matrix of simplex vertex coordinates

Returns:
    X: Matrix of quadrature points, where each row contains the
       coordinates for a point
    W: Array of quadrature weights

Sample Usage:

    X, W = simplexquad(n, vert)     # Specify the vertices
    X, W = simplexquad(T, n, dim)   # Specify the dimension and use unit simplex

`X` will be a matrix for which the `j`th column are the grid points in
each coordinate `xj`.

To integrate a function `f`, use e.g.

    sum(W[i] * f(X[i,:]) for i in 1:length(W))

Note: The standard `n`-dimensional simplex has vertices specified
`vert = diagm(n+1, n, ones(n))
"""
simplexquad

function simplexquad(::Type{T}, N::Int, n::Int) where {T}
    m = n + 1
    vert = diagm(m, n, ones(T, n))
    simplexquad(N, vert)
end
simplexquad(N::Int, n::Int) = simplexquad(float(Int), N, n)

function simplexquad(N::Int, vert::AbstractMatrix{T}) where {T}
    N <= 0 && @error "Number of Gauss points must be a natural number"
    
    m, n = size(vert)
    m != n+1 && @error "The matrix of vertices must have n+1 rows and n columns"
    
    Nn = N^n
    if n == 1
        # The 1-D simplex is only an interval
        q, w = rquad(T, N, 0)
        len = diff(vert; dims=1)
        X = vert[1] .+ len[1] * q
        W = abs.(len[1]) * w
        X = reshape(X, :, 1)
    else
        # Find quadrature rules for higher dimensional domains
        q = Array{Vector}(undef, n)
        w = Array{Vector}(undef, n)
        for k in 1:n 
            q[k], w[k] = rquad(T, N, n - k)
        end
        Q = ndgrid(q...)
        q = reshape(cat(Q...; dims=n), Nn, n)
        W = ndgrid(w...)
        w = reshape(cat(W...; dims=n), Nn, n)
        map = zeros(T, m, m) + I
        map[2:m, 1] .= -1
        c = map * vert
        W = abs(det(c[2:m, :])) * prod(w; dims=2)
        qp = cumprod(q; dims=2)
        e = ones(T, Nn, 1)
        X = [e [(1 .- q[:, 1:n-1]) e] .* [e qp[:, 1:n-2] qp[:, n]]] * c
        @assert size(W, 2) == 1
        W = reshape(W, :)
    end

    X::Matrix{T}
    @assert size(X, 2) == n
    W::Vector{T}
    @assert length(W) == size(X, 1)
    X, W
end

function rquad(::Type{T}, N::Int, k::Int) where {T}
    k1 = k+1
    k2 = k+2

    n = collect(1:N)'
    nnk = 2*n .+ k
    A = [T(k)/k2 repeat([T(k^2)], 1, N) ./ (nnk .* (nnk .+ 2))]

    n = collect(2:N)'
    nnk = nnk[n]
    B1 = T(4)*k1 / (k2*k2*(k+3))

    nk = n .+ k
    nnk2 = nnk .* nnk
    B = T(4)*(n .* nk).^2 ./ (nnk2 .* nnk2 - nnk2)

    ab = [A' [T(2^k1)/k1; B1; B']]
    s = sqrt.(ab[2:N, 2])
    X, V = eigen(SymTridiagonal(ab[1:N, 1], s))
    I = sortperm(X)
    X = X[I]
    x = (X .+ 1)/2
    w = (T(1)/2)^k1 * ab[1,2] * V[1,I].^2

    x::Vector{T}
    @assert length(x) == N
    w::Vector{T}
    @assert length(w) == N
    x, w
end



# These functions were formerly a part of Julia. License is MIT:
# https://julialang.org/license

ndgrid(v::AbstractVector) = copy(v)

function ndgrid(v1::AbstractVector{T}, v2::AbstractVector{T}) where T
    m, n = length(v1), length(v2)
    v1 = reshape(v1, m, 1)
    v2 = reshape(v2, 1, n)
    (repeat(v1, 1, n), repeat(v2, m, 1))
end

function ndgrid_fill(a, v, s, snext)
    for j = 1:length(a)
        a[j] = v[div(rem(j-1, snext), s)+1]
    end
end

function ndgrid(vs::AbstractVector{T}...) where T
    n = length(vs)
    sz = map(length, vs)
    out = ntuple(i->Array{T}(undef, sz), n)
    s = 1
    for i=1:n
        a = out[i]::Array
        v = vs[i]
        snext = s*size(a,i)
        ndgrid_fill(a, v, s, snext)
        s = snext
    end
    out
end



integrate(f, X, W) = sum(W[i] * f(X[i,:]) for i in 1:length(W))

end
