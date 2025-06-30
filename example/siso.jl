using SysID

using Plots
using DSP: mt_coherence
using ControlSystems: tf, bodeplot, setPlotScale
using ControlSystemIdentification: FRD

include("utils/chirp.jl")
include("utils/generateData.jl")

tEnd = 30
dt = 1/100
t = 0:dt:tEnd
sysTrue = tf([5.0], [0.2, 1.0])

chirpFreq = [0.1, 10]

t, u = chirp(chirpFreq[1], chirpFreq[2], tEnd, dt, shape=:log)
u = reshape(u, 1, :) # convert `u` to [1, length(t)], a matrix
y, ~, ~ = generateData(sysTrue, u, t; noise=false)

f, H = estimate_frqrsp(vec(y), vec(u), dt)
sysFrd = FRD(f, H)
#sysEst = fit_tf(f, H, 1)

ws = range(chirpFreq[1],stop=chirpFreq[2],length=200)
bplt = bodeplot(sysTrue, ws, label="True TF")
#bodeplot!(sysEst, ws, label="Estimated TF")
setPlotScale("dB")
display(bplt)

plt = plot(layout=(3,1), size=(800, 600))
plot!(plt[1], sysFrd.w, 20 .* log10.(abs.(sysFrd.r)), 
        xscale=:log10, 
        xlims=chirpFreq,
        ylabel="Magnitude (dB)",
        title="Bode Plot",
        legend=false,
        grid=true)
plot!(plt[2], sysFrd.w, angle.(sysFrd.r) .* 180 ./ π, 
        xscale=:log10, 
        xlims=chirpFreq,
        xlabel="Frequency (rad/s)", 
        ylabel="Phase (deg)",
        legend=false,
        grid=true)
data = [u; y]
coh = mt_coherence(data, fs=1/dt)
plot!(plt[3], coh.freq[2:end] .* 2pi, coh.coherence[2,1,2:end], 
        xscale=:log10, 
        xlims=chirpFreq,
        xlabel="Frequency (rad/s)", 
        ylabel="Coherence",
        legend=false,
        grid=true)
display(plt)

plot(t, vec(u))
