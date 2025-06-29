# System Identification in Julia using Frequency Response
# SISO, MISO, and MIMO examples with TF/State-space estimation

using ControlSystems, DSP, Plots, LinearAlgebra, Optim, FFTW

# Manual logarithmic chirp generator (no dependencies)
function manual_chirp(t::AbstractVector{<:Real}, f0::Real, f1::Real)
    T = t[end]
    β = log(f1 / f0) / T
    return sin.(2π * f0 * (exp.(β .* t) .- 1) ./ β)
end

# Helper: Estimate Frequency Response Function (FRF)
function estimate_frf(y::Vector, u::Vector, fs::Float64; nfft=1024)
    Y = dspfft(y, nfft)
    U = dspfft(u, nfft)
    Pxy = Y .* conj.(U)
    Pxx = U .* conj.(U)
    H = Pxy ./ Pxx
    f = fs .* (0:(nfft ÷ 2)) ./ nfft
    return f[1:length(H)], H[1:length(H)]
end

function dspfft(x::Vector, nfft::Int)
    X = fft(x, nfft)[1:(nfft ÷ 2 + 1)]
    return X ./ length(x)
end

# Helper: Fit Transfer Function to Estimated FRF
function fit_tf(f::Vector, H::Vector, order::Int)
    ω = 2pi .* f
    s = im .* ω
    Hc = H

    if order == 1
        model1(s, p) = p[1] ./ (p[2] .* s .+ 1)
        loss1(p) = sum(abs2, model1(s, p) .- Hc)
        result = optimize(loss1, [1.0, 1.0])
        p_opt = Optim.minimizer(result)
        return tf([p_opt[1]], [p_opt[2], 1.0])
    elseif order == 2
        model2(s, p) = p[1] .* p[3]^2 ./ (s.^2 .+ 2 .* p[2] .* p[3] .* s .+ p[3]^2)
        loss2(p) = sum(abs2, model2(s, p) .- Hc)
        result = optimize(loss2, [1.0, 0.7, 2.0])
        p_opt = Optim.minimizer(result)
        num = [p_opt[1]*p_opt[3]^2]
        den = [1.0, 2*p_opt[2]*p_opt[3], p_opt[3]^2]
        return tf(num, den)
    else
        error("Only order 1 or 2 supported")
    end
end

# SISO Example: Roll Dynamics
function siso_example()
    fs, T = 100.0, 30.0
    t = 0:1/fs:T
    u = manual_chirp(t, 0.1, 5.0)
    u = reshape(u, 1, :)
    sys_true = tf([5.0], [0.2, 1.0])
    y = lsim(sys_true, u, t).y[1, :]

    f, H = estimate_frf(y, vec(u), fs)
    sys_est = fit_tf(f, H, 1)

    bodeplot(sys_est, label="Estimated TF")
    bodeplot!(sys_true, label="True TF")
end

# MISO Example: Two Inputs → One Output
function miso_example()
    fs, T = 100.0, 30.0
    t = 0:1/fs:T
    u1 = manual_chirp(t, 0.1, 3.0)
    u2 = manual_chirp(t, 1.0, 5.0)

    u1 = reshape(u1, 1, :)
    u2 = reshape(u2, 1, :)

    sys1 = tf([2.0], [0.5, 1.0])
    sys2 = tf([1.0], [1.0, 1.0])
    y = lsim(sys1, u1, t).y[1, :] .+ lsim(sys2, u2, t).y[1, :] + 0.01*randn(length(t))

    f1, H1 = estimate_frf(y, vec(u1), fs)
    f2, H2 = estimate_frf(y, vec(u2), fs)

    sys1_est = fit_tf(f1, H1, 1)
    sys2_est = fit_tf(f2, H2, 1)

    return sys1_est, sys2_est
end

# MIMO Example: Two Inputs → Two Outputs
function mimo_example()
    fs, T = 100.0, 30.0
    t = 0:1/fs:T
    u1 = manual_chirp(t, 0.1, 3.0)
    u2 = manual_chirp(t, 1.0, 5.0)

    u1 = reshape(u1, 1, :)
    u2 = reshape(u2, 1, :)

    sys11 = tf([2.0], [0.5, 1.0])
    sys12 = tf([1.0], [1.0, 1.0])
    sys21 = tf([1.5], [0.4, 1.0])
    sys22 = tf([2.0], [0.8, 1.0])

    y1 = lsim(sys11, u1, t).y[1, :] .+ lsim(sys12, u2, t).y[1, :]
    y2 = lsim(sys21, u1, t).y[1, :] .+ lsim(sys22, u2, t).y[1, :]

    f, H11 = estimate_frf(y1, vec(u1), fs)
    _, H12 = estimate_frf(y1, vec(u2), fs)
    _, H21 = estimate_frf(y2, vec(u1), fs)
    _, H22 = estimate_frf(y2, vec(u2), fs)

    return H11, H12, H21, H22
end

# Run SISO Example
siso_example()
