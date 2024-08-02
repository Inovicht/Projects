function [frec, amp] = amplitudeSpectrumOneSided(t,x)

F = fft(x);

n = length(t);

Af = 1/(t(end)-t(1));

for i = 1:n
    f(i) = (i-1)*Af;
    C(i) = F(i)/n;
end

n2 = floor(n/2);

frec = f(1:n2);
amp = [abs(C(1)) 2*abs(C(2:n2))];

% figure;
semilogy(frec,amp);
title('AMPLITUDE SPECTRUM');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

