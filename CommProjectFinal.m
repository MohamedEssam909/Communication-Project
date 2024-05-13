num_bits = 64;  % Number of bits
bits = randi([0, 1], 1, num_bits);  % Generate random bit stream
bit_rate = 1;

T = length(bits)/bit_rate; % Full time of bit sequence
n = 200;
N = n*length(bits);
dt = T/N;
t = 0:dt:T;
x = zeros(1,length(t)); % Output signal
for i = 0:length(bits)-1
  if bits(i+1) == 1
    x(i*n+1:(i+1)*n) = 1;
  else
    x(i*n+1:(i+1)*n) = 0;
  end
end

graphics_toolkit('gnuplot');

figure;
plot(t, x, 'LineWidth', 2);
xlabel('Time');
ylabel('Amplitude');
title('UNRZ Encoded Signal');

ylim([0, 1.1]);  % Set the y-axis limit from 0 to 1.1
xlim([0, num_bits]);



% Plot the spectral domain representation (Fourier transform)
sampling_freq = 1 / dt;
N = length(x);  % Length of the signal
frequencies = (-sampling_freq/2):(sampling_freq/N):(sampling_freq/2 - sampling_freq/N);  % Frequency axis
spectrum = abs(fftshift(fft(x))) / N;  % FFT
figure;
plot(abs(frequencies), spectrum * 2, 'r', 'LineWidth', 3);
xlabel('Frequency (kHz)');
ylabel('Magnitude');
title('Spectral Domain Representation of Unipolar NRZ Encoded Signal');

max_magnitude = max(spectrum);
ylim([0, 1.1 * max_magnitude]);
xlim([-10, 10]);



