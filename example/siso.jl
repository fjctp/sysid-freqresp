using SysID

using Plots
using ControlSystems: tf, bodeplot, setPlotScale
using ControlSystemIdentification: FRD

include("utils/chirp.jl")
include("utils/generateData.jl")

tEnd = 30
dt = 1/100
t = 0:dt:tEnd
sysTrue = tf([5.0], [0.2, 1.0])

chirpFreq = [0.1, 10]

t, u = chirp(chirpFreq[1], chirpFreq[2], tEnd, dt, shape=:quad)
u = reshape(u, 1, :) # convert `u` to [1, length(t)], a matrix
res = lsim(sysTrue, u, t) 
y = res.y # Size: [1, length(t)], a matrix

f, H = estimate_frqrsp(vec(y), vec(u), dt)
sysFrd = FRD(f, H)
#sysEst = fit_tf(f, H, 1)

ws = range(chirpFreq[1],stop=chirpFreq[2],length=200)
bodeplot(sysTrue, ws, label="True TF")
#plot(sysFrd, hz=false)
#bodeplot!(sysEst, ws, label="Estimated TF")
setPlotScale("dB")
