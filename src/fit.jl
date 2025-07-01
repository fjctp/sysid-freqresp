using Optim
using ControlSystems: freqresp

# Helper: Fit Transfer Function to Estimated FRF
function fit_tf(w::AbstractVector{<:Real}, H::AbstractVector{<:Complex}, 
    x0::AbstractVector{<:Real}, model)

    function loss(x)
        frqrsp = freqresp(model(x), w)
        frqrsp = frqrsp[1,1,:]
        frqrspErr = frqrsp .- H
        
        return 0.5 * sum(abs, frqrspErr .* conj(frqrspErr))
    end

    res = optimize(loss, [0.1, 0.1], [Inf, Inf], x0)
    display(res)

    return model(res.minimizer)
end
