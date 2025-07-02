using Plots
using ControlSystems
using ControlSystemIdentification

include("utils/chirp.jl")
include("utils/generateData.jl")

# Generate test data
tEnd = 30
dt = 1/100
t = 0:dt:tEnd
sysTrue = tf([2.0], [0.5, 1.0])

chirpFreq = [0.1, 10] # rad/s

t, u = chirp(chirpFreq[1], chirpFreq[2], tEnd, dt, shape=:log)
u = reshape(u, 1, :) # convert `u` to [1, length(t)], a matrix
y, ~, ~ = generateData(sysTrue, u, t; noise=false)

# Estimate frequency response with time domain data.
idd = iddata(y, u, dt)
coh = coherence(idd)

pltTime = plot(idd)
pltCoh = plot(coh, yscale=:identity, grid=true)

display(pltTime)

display(pltCoh)
