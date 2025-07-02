using SysID

using Plots
using DSP: mt_coherence
using ControlSystems: tf, bodeplot, bodeplot!, setPlotScale
using ControlSystemIdentification: FRD

include("utils/chirp.jl")
include("utils/generateData.jl")

# Generate test data
tEnd = 30
dt = 1/100
t = 0:dt:tEnd
sysTrue = tf([2.0], [0.5, 1.0])

chirpFreq = [0.1, 10]

t, u = chirp(chirpFreq[1], chirpFreq[2], tEnd, dt, shape=:log)
u = reshape(u, 1, :) # convert `u` to [1, length(t)], a matrix
y, ~, ~ = generateData(sysTrue, u, t; noise=false)

# Estimate frequency response with time domain data.
w, H = estimate_frqrsp(vec(y), vec(u), dt)
sysFrd = FRD(w, H)

# Estimate transfer function from estimated frequency response.
K0, tau0 = rand(2)
x0 = [K0, tau0]
sysEst = fit_tf(w, H, x0, x -> tf(x[1], [x[2], 1.]))

# Plot bode
# - Truth system
ws = range(chirpFreq[1],stop=chirpFreq[2],length=200)
setPlotScale("dB")
pltBode = bodeplot(sysTrue, ws, label="True TF")

# - Estimated frequency response
plot!(pltBode[1], sysFrd.w, 20 .* log10.(abs.(sysFrd.r)), 
        label="Estimated Frd", 
        xscale=:log10, 
        xlims=chirpFreq,
        grid=true)
plot!(pltBode[2], sysFrd.w, angle.(sysFrd.r) .* 180 ./ π, 
        label="Estimated Frd", 
        xscale=:log10, 
        xlims=chirpFreq,
        grid=true)

# - Estimated system
bodeplot!(sysEst, ws, label="Estimated TF")

# Plot coherence
data = [u; y]
coh = mt_coherence(data, fs=1/dt)
pltCoh = plot(layout=(1,1), size=(800, 600))
plot!(pltCoh[1], coh.freq[2:end] .* 2pi, coh.coherence[2,1,2:end], 
        xscale=:log10, 
        xlims=chirpFreq,
        xlabel="Frequency (rad/s)", 
        ylabel="Coherence",
        legend=false,
        grid=true)

# Plot time-domain
pltTime = plot(t, vec(u))

display(pltBode)
display(pltCoh)
