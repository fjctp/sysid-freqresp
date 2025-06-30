using ControlSystems: tf, lsim

function generateData(sys, u, t; noise::Bool=false)
    res = lsim(sys, u, t)
    return res.y, u, t
end
