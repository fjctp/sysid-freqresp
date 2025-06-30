using DSP

# Helper: Estimate Frequency Response Function (FRF)
function estimate_frqrsp(y::Vector, u::Vector, dt::Real; nfft=1024)
    fs = 1/dt
    nfft = nextpow(2, length(u))

    # PSD using Multitaper
    psd_uu = mt_pgram(u; fs=fs, nfft=nfft, nw=4)
    freq_uu = freq(psd_uu) * 2 * pi
    p_uu = power(psd_uu)

    # Cross-PSD using Multitaper
    # Stack signals in matrix: channels × samples
    UY = [reshape(u, 1, :); reshape(y, 1, :)]
    cpsd_uy = mt_cross_power_spectra(UY; fs=fs, nfft=nfft, nw=4)
    p_uy = cpsd_uy.power[2,1,:] # cross-spectrum U and Y

    # Estimate FRF: H = S_yu / S_uu
    H = p_uy ./ p_uu

    return freq_uu[2:end], H[2:end]
end
