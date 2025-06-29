using DSP

# Helper: Estimate Frequency Response Function (FRF)
function estimate_frqrsp(y::Vector, u::Vector, fs::Float64; nfft=1024)

    # PSD using Welch
    Pxx = DSP.Periodograms.welch_pgram(u; fs=fs, nfft=512, window=hanning)
    freq = DSP.Periodograms.freq(Pxx)
    Puu = DSP.Periodograms.power(Pxx)

    # Cross-PSD using multitaper
    # Stack signals in matrix: channels × samples
    XY = [y; u]
    cpsd_mt = DSP.Periodograms.mt_cross_power_spectra(XY; fs=fs, nfft=512, nw=4)
    C = DSP.Periodograms.power(cpsd_mt)   # dims: 2×2×length(freq)
    Pyu = C[1,2,:]                       # cross-spectrum Y and U

    # Estimate FRF: H = S_yu / S_uu
    H = Pyu ./ Puu

    return freq, H
end
