function chirp(w0::Real, w1::Real, duration::Real, dt::Real; shape::Symbol = :log)
    t = 0:dt:duration

    if shape == :linear
        k = (w1 - w0) / duration
        phase = w0 .* t .+ 0.5 .* k .* t.^2

    elseif shape == :log
        β = log(w1 / w0) / duration
        phase = (w0 / β) .* (exp.(β .* t) .- 1)

    else
        error("Unsupported chirp shape. Use :linear, :quad, or :exp")

    end

    return t, sin.(phase)
end
