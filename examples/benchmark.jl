using BenchmarkTools
using SimplexQuad

# @fastmath function integrate2d(f, X, W)
#     if size(X, 2) == 2
#         @inbounds sum(W[i] * f(X[i,1], X[i,2]) for i in 1:length(W))
#     elseif size(X, 2) == 2
#         @inbounds sum(W[i] * f(X[i,1], X[i,2], X[i,3]) for i in 1:length(W))
#     else
#         @error "not implemented"
#     end
# end

@fastmath function integrate(f::F, X, W) where {F}
    @inbounds begin
        if size(X, 2) == 2
            s = zero(W[1] * f(X[1, 1], X[1, 2]))
            for i = 1:length(W)
                s += W[i] * f(X[i, 1], X[i, 2])
            end
            s
        elseif size(X, 2) == 3
            s = zero(W[1] * f(X[1, 1], X[1, 2], X[1, 3]))
            for i = 1:length(W)
                s += W[i] * f(X[i, 1], X[i, 2], X[i, 3])
            end
            s
        else
            @error "not implemented"
        end
    end
end



function bench_integrate(outfile, D, N, f, nflop_f)
    X, W = simplexquad(N, D)
    t = @benchmark integrate($f, $X, $W)
    nflop = (nflop_f + 2) * length(W)
    nflop_sec = nflop / minimum(t).time
    println(outfile,
            "D=$D N=N   $nflop Flop   $(round(nflop_sec; digits=1)) GFlop/s")
end

open("figures/benchmarks.txt", "w") do outfile
    f2(x, y) = @fastmath 1 + x + y^2
    nflop_f2 = 3
    bench_integrate(outfile, 2, 4, f2, nflop_f2)
    bench_integrate(outfile, 2, 8, f2, nflop_f2)

    f3(x, y, z) = @fastmath 1 + x + y^2 + z^3
    nflop_f3 = 6
    bench_integrate(outfile, 3, 4, f3, nflop_f3)
    bench_integrate(outfile, 3, 8, f3, nflop_f3)
end
