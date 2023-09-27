clear;
%% Params
n_subcarrier = 256;
freq_carrier = 3e9;
freq_IF = 2e7;
freq_shift = 0;
lambda = 3e8 / freq_carrier;
T_symbol = 1e-4;
N_sample = 8192; %8192
symbols_per_frame = 2;
pri = T_symbol*symbols_per_frame;
freq_sampling = N_sample / T_symbol;
range_per_sampling_period = 3e8 * (1/freq_sampling) / 2;
rcs = 1;
Gt = 2936;
Gr = 2936;
Bn = (1/T_symbol)/symbols_per_frame;                        % BW post doppler
t = (linspace(0, T_symbol, n_subcarrier))';                 %symbol time
t2 = (linspace(0, T_symbol, N_sample))';                    %symbol time interpolated
pri_t = (linspace(0, pri, N_sample*symbols_per_frame))';    %pri time

%% Wave Generator
%golay coding
golay_reg = golay(n_subcarrier);
golay_symbol_index = 1;
golay_code_value = golay_reg(golay_symbol_index, :);

%phase shifting
% alpha = deg2rad(90);
% deltat = -alpha ./ (2*pi*f);
% phase_shift = exp(1j * 2 * pi * f .* deltat);
% signal = golay_code_value .* phase_shift;

%ifft
iq_signal = ifft(golay_code_value, n_subcarrier)';

if N_sample > n_subcarrier
    a = 1:1:n_subcarrier;
    b = linspace(1, n_subcarrier, N_sample);
    iq_signal = interp1(a,iq_signal,b)';
end

% %offset-freq and IF Freq signal 
carrier = exp(-1i*2*pi*(freq_IF+freq_shift)*t2);
iq_signal = iq_signal .* carrier;

%IQ modulator
i_signal = real(iq_signal);
q_signal = imag(iq_signal);
tx_signal = (i_signal + q_signal);

% periodogram(tx_signal,[],512,freq_sampling,'centered')

% Showing time-domain symbol
% figure(1);
% plot(t, tx_signal);
% xlabel('Time (s)');
% ylabel('Amplitude (V)');
% title('OFDM Symbol in Time Domain');
% 
% figure(9);
% plot(t2, z);
% xlabel('Time (s)');
% ylabel('Amplitude (V)');
% title('OFDM Symbol Interpolated in Time Domain');

% showing freq-domain symbol
% fft_tx_signal = fft(tx_signal);
% freq_axis = (0:n_subcarrier-1)*(freq_sampling/n_subcarrier);
% power = abs(fft_tx_signal).^2/n_subcarrier;
% figure(2);
% plot(freq_axis,power);
% xlabel('Frequency');
% ylabel('Power');
% title('OFDM Symbol in Frequency Domain');

% fft_z_signal = fft(tx_signal);
% freq_axis = (0:N_sample-1)*(freq_sampling/N_sample);
% power = abs(fft_z_signal).^2/N_sample;
% fft_z_signal_center = fftshift(fft_z_signal);
% freq_axis_centered = (-N_sample/2:N_sample/2-1)*(freq_sampling/N_sample); 
% power_centered = abs(fft_z_signal_center).^2/N_sample;
% 
% figure(2);
% plot(freq_axis_centered,power_centered);
% xlabel('Frequency');
% ylabel('Power');
% title('OFDM Symbol in Frequency Domain');

% scope = spectrumAnalyzer(SampleRate=freq_sampling,AveragingMethod="exponential",...
%     PlotAsTwoSidedSpectrum=false,...
%     RBWSource="auto",SpectrumUnits="dBW");
% scope(tx_signal)

% powerbw(tx_signal,freq_sampling)
%% Channel
% zero padding and delaying signal
delay_lag = 500;
% delayed_symbol_signal = circshift(tx_signal, mod(delay_lag, N_sample));
ref_signal = [tx_signal; zeros((N_sample*(symbols_per_frame-1)), 1)];
zero_padded_signal = [tx_signal; zeros((N_sample*(symbols_per_frame-1)), 1)];
delayed_signal = circshift(zero_padded_signal, delay_lag);
if delay_lag > N_sample*(symbols_per_frame-1) && symbols_per_frame ~= 1
    delayed_signal(1:N_sample-1) = 0;
end

% attenuating signal (free-space loss in Voltage)
Lfs = sqrt(((4*pi)^3) * ((range_per_sampling_period*delay_lag)^4) / ((lambda^2) * rcs * Gt * Gr));
attenuated_signal = delayed_signal * Lfs / Lfs;

% noising signal
N0 = 1.38e-23 * 290 * Bn; %temperature noise floor
corrupted_signal = attenuated_signal; %awgn(attenuated_signal, 100);

% Showing time-domain channeled signal
% figure(3);
% plot(pri_t, corrupted_signal);
% xlabel('Time (s)');
% ylabel('Amplitude (V)');
% title('Received Signal Across PRI');

% figure(10);
% plot(pri_t, ref_signal);
% xlabel('Time (s)');
% ylabel('Amplitude (V)');
% title('Reference Signal');
%% RSP
rx_signal = corrupted_signal;

%cross-correlation
[c,lags] = xcorr(rx_signal, ref_signal);

%range processing
size_lags =  size(lags, 2);
lags_abs = lags(1, round(size_lags/2):size_lags);
c_abs = c(round(size_lags/2):size_lags, 1);

% rxy = c_abs/max(c_abs);
rxy = c_abs;
tau = lags_abs * T_symbol / max(lags_abs);
lag_of_max_c = lags_abs(c_abs == max(c_abs))
% phase_of_max_c = (lag_of_max_c / size(lags_abs, 2)) * 2 * pi;

lag_delay_time = (1/freq_sampling) * lag_of_max_c;
range_target = 3e8 * lag_delay_time / 2

% figure(4);
% stem(tau, rxy);
% xlabel('\tau (s)');
% ylabel('R_x_y');
% title('Cross Correlation of Transmitted and Received Signal');

