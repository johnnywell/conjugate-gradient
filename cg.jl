trans = transpose

function cg(A, b, x; ϵ::Float64=1e-8, imax::Int=100000)
    i = 1
    r = b - (A * x)
    d = r
    δ = dot(r, r)
    δ₀ = δ
    while i < imax && δ > (ϵ^2 * δ₀)
        q = A * d
        α = δ / dot(d, q)
        x = x + (α * d)
        if i % 50 == 0
            r = b - (A * x)
        else
            r = r - (α * q)
        end
        δₒ = δ
        δ = dot(r, r)
        β = δ/δₒ
        d = r + (β * d)
        i += 1
    end
    return x, i
end
