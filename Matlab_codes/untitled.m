figure
Fs = T/dt;  % Sampling rate in Hz
N = length(mean(sensor_data,1));
fft_r = zeros(1,N);
for i = 1:Nx
    disp(i)
    fft_r(i,:) = exp(-abs(i-250)*1i)*fft(sensor_data(i,:));
end
fft_result = mean(fft_r);
frequencies = (0:N-1) * (Fs / N);
frequencies = frequencies(1:N/2);
fft_result = fft_result(1:N/2);
fft_result = abs(fft_result);
loglog(frequencies, fft_result);
title('Frequency Domain Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude(a.u.)');
grid on
grid minor