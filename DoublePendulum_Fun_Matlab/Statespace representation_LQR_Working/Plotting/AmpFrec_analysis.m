function AmpFrec_analysis(t, y)
    % Calculate the amplitude spectrum
    [frec, amp] = amplitudeSpectrumOneSided(t, y);

    % Create the plot
    figure;
    semilogy(frec, amp);

    % Add title and labels
    title('Amplitude Spectrum Analysis');
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    
    % Add grid
    grid on;
end