using ControlSystems: tf, lsim

function generateData(sys, u, t; noise::Bool=false)
    res = lsim(sys, u, t)
    y = res.y
    if noise
        y .+= rand(size(y))
    end
    return y, u, t
end
