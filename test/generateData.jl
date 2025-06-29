using ControlSystems: tf, lsim

function generateData(sys, t, u; noise::Bool=false)
    y = lsim(sys, u, t)
    return y, u, t
end
