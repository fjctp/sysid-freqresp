function chirp(f0::Real, f1::Real, duration::Real, dt::Real; shape::Symbol = :log)
    t = 0:dt:duration
    if shape == :linear
        k = (f1 - f0) / duration
        phase = f0 .* t .+ 0.5 .* k .* t.^2
    elseif shape == :quad
        a2 = 1
        a1 = (f1 - f0 - a2*duration^2) / duration
        phase = f0 .* t .+ 0.5 .* a1 .* t.^2 .+ (1/3) .* a2 .* t.^3

    elseif shape == :log
        ω0, ω1 = f0, f1
        β = log(ω1 / ω0) / duration
        phase = (ω0 / β) .* (exp.(β .* t) .- 1)
    else
        error("Unsupported chirp shape. Use :linear, :quad, or :exp")
    end
    return t, sin.(phase)
end
