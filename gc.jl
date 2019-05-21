trans = transpose

function gc(A, b, x; ϵ::Float64=10e-15, imax::Int=100000)
    i = 1
    r = b - (A * x)
    d = r
    δ = trans(r) * r
    δ₀ = δ
    while i < imax && δ > (ϵ^2 * δ₀)
        q = A * d
        α = δ / (trans(d) * q)
        x = x + (α * d)
        if i % 50 == 0
            r = b - (A * x)
        else
            r = r - (α * q)
        end
        δₒ = δ
        δ = trans(r) * r
        β = δ/δₒ
        d = r + (β * d)
        i += 1
    end
    return x, i
end
