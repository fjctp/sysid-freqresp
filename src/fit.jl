
# Helper: Fit Transfer Function to Estimated FRF
function fit_tf(f::Vector, H::Vector, order::Int)
    ω = 2pi .* f
    mag = abs.(H)
    phase = angle.(H)

    # Convert to complex frequency domain data
    s = im .* ω
    Hc = H

    # Fit TF using ControlSystems.jl
    # Use vector fitting or other method for robust fitting
    # Here we do simple least-squares fit to 1st or 2nd order TF
    if order == 1
        # Fit H(s) = k / (τs + 1)
        function model(s, p)
            k, τ = p
            return k ./ (τ .* s .+ 1)
        end
        function loss(p)
            sum(abs2, model(s, p) .- Hc)
        end
        p_opt = optimize(loss, [1.0, 1.0]).minimizer
        return tf([p_opt[1]], [p_opt[2], 1.0])
    elseif order == 2
        # Fit H(s) = kω² / (s² + 2ζωs + ω²)
        function model(s, p)
            k, ζ, ωn = p
            return k .* ωn^2 ./ (s.^2 .+ 2ζ .* ωn .* s .+ ωn^2)
        end
        function loss(p)
            sum(abs2, model(s, p) .- Hc)
        end
        p_opt = optimize(loss, [1.0, 0.7, 2.0]).minimizer
        num = [p_opt[1]*p_opt[3]^2]
        den = [1.0, 2*p_opt[2]*p_opt[3], p_opt[3]^2]
        return tf(num, den)
    else
        error("Only order 1 or 2 supported")
    end
end
