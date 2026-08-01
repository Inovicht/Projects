function AmpFrec_analysis(t,y)
    [frec, amp] = amplitudeSpectrumOneSided(t,y);
    semilogy(frec,amp);
end