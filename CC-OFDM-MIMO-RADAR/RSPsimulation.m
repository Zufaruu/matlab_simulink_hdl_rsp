%% Wave Generator
%golay coding
golay_reg = golay(256);
golay_symbol_index = 1;
golay_code_value = golay_reg(golay_symbol_index, :);

n_subcarrier = 256;
f_start = 20e6;
f_end = 22.56e6;
T_symbol = 1e-4;
f = (linspace(f_start, f_end, n_subcarrier))'; %subcarrier frequencies
t = (linspace(0, T_symbol, n_subcarrier))'; %signal time

%phase shifting
% alpha = deg2rad(90);
% deltat = -alpha ./ (2*pi*f);
% phase_shift = exp(1j * 2 * pi * f .* deltat);
% signal = golay_code_value .* phase_shift;

%ifft
iq_signal = ifft(golay_code_value, n_subcarrier);
i_signal = real(iq_signal);
q_signal = imag(iq_signal);
tx_signal = i_signal + q_signal;

% plot(t, tx_signal);

%% Channel
delayed_signal = circshift(tx_signal, 128);
corrupted_signal = awgn(delayed_signal, 30);

%% RSP
rx_signal = corrupted_signal;

%cross-correlation
[c,lags] = xcorr(rx_signal, tx_signal);

%range processing
size_lags =  size(lags, 2);
lags_abs = lags(1, round(size_lags/2):size_lags);
c_abs = c(1, round(size_lags/2):size_lags);

lag_of_max_c = lags_abs(c_abs == max(c_abs));
phase_of_max_c = (lag_of_max_c / max(lags_abs)) * 2 * pi;
stem(lags_abs,c_abs)

r_ua = 3e8 * T_symbol / 2;
range_target = (phase_of_max_c / (2*pi)) * r_ua;
range_target



