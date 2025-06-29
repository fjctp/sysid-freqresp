using SysID
using ControlSystems: bodeplot

include("chirp.jl")
include("generateData.jl")

tEnd = 30
dt = 1/100
t = 0:dt:tEnd
sysTrue = tf([5.0], [0.2, 1.0])

u = chirp(0.1, 10, tEnd, dt, shape=:quad)
u = reshape(u, 1, :) # convert to [1, length(t)]
res = lsim(sys_true, u, t)
y = res.y[1, :]

f, H = estimate_frqrsp(y, vec(u), fs)
sys_est = fit_tf(f, H, 1)

bodeplot(sys_est, label="Estimated TF")
bodeplot!(sys_true, label="True TF")
