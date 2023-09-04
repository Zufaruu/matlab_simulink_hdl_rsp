%% Wave Generator
%golay coding
golay_reg = golay(256);
golay_symbol_index = 1;
golay_code_value = golay_reg(golay_symbol_index, :);

n_subcarrier = 256;
f_start = 20e6;
f_end = 22.56e6;
T_symbol = 1e-4;
symbol_delay = 12;
freq_sampling = n_subcarrier / T_symbol;
f = (linspace(f_start, f_end, n_subcarrier))'; %subcarrier frequencies
t = (linspace(0, T_symbol, n_subcarrier))'; %signal time
pri_t = (linspace(0, T_symbol*symbol_delay, n_subcarrier*symbol_delay))';

%phase shifting
% alpha = deg2rad(90);
% deltat = -alpha ./ (2*pi*f);
% phase_shift = exp(1j * 2 * pi * f .* deltat);
% signal = golay_code_value .* phase_shift;

%ifft
iq_signal = ifft(golay_code_value, n_subcarrier)';

% %offset-freq signal
% freq_shift = 5000
% iq_signal = frequencyOffset(iq_signal,freq_sampling,freq_shift);

%IQ modulator
i_signal = real(iq_signal);
q_signal = imag(iq_signal);
tx_signal = (i_signal + q_signal);

% periodogram(tx_signal, fs)
% periodogram(tx_signal,[],n_subcarrier,freq_sampling,'centered')

% Showing time-domain symbol
% figure(1);
% plot(t, tx_signal);
% xlabel('Time (s)');
% ylabel('Amplitude (V)');
% title('OFDM Symbol in Time Domain');

% % showing freq-domain symbol
% fft_tx_signal = fft(tx_signal);
% freq_axis = (0:n_subcarrier-1)*(freq_sampling/n_subcarrier);
% power = abs(fft_tx_signal).^2/n_subcarrier;
% figure(2);
% plot(freq_axis,power);
% xlabel('Frequency');
% ylabel('Power');
% title('OFDM Symbol in Frequency Domain');

%% Channel
% zero padding and delaying signal
delay_lag = 1;
delayed_symbol_signal = circshift(tx_signal, mod(delay_lag, n_subcarrier));
zero_padded_signal = [tx_signal; zeros((n_subcarrier*(symbol_delay-1)), 1)];
delayed_signal = circshift(zero_padded_signal, delay_lag);

if delay_lag > 256*11
    delayed_signal(1:255) = 0;
end

% noising signal
corrupted_signal = awgn(delayed_signal, 50);

% Showing time-domain channeled signal
figure(3);
plot(pri_t, corrupted_signal);
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Received Signal Along PRI');
%% RSP
rx_signal = corrupted_signal;

%cross-correlation
[c,lags] = xcorr(rx_signal, zero_padded_signal);

%range processing
size_lags =  size(lags, 2);
lags_abs = lags(1, round(size_lags/2):size_lags);
c_abs = c(round(size_lags/2):size_lags, 1);

rxy = c_abs/max(c_abs);
tau = lags_abs * T_symbol / max(lags_abs);
lag_of_max_c = lags_abs(c_abs == max(c_abs));
phase_of_max_c = (lag_of_max_c / size(lags_abs, 2)) * 2 * pi;

figure(4);
stem(tau, rxy);
xlabel('\tau (s)');
ylabel('R_x_y');
title('Cross Correlation of Transmitted and Received Signal');

lag_delay_time = (1/freq_sampling) * lag_of_max_c;
range_target = 3e8 * lag_delay_time / 2;
lag_of_max_c
range_target