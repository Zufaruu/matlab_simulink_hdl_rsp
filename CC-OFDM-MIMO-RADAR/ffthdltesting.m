% obs_samples = [0.0019
%     0.0041
%     0.0006
%    -0.0036
%    -0.0028
%     0.0018
%     0.0039
%     0.0005
%    -0.0034
%    -0.0026
%    0
%    0
%    0
%    0
%    0
%    0
%    0
%    0
%    0
%    0
%    0
%    0
%    0
%    0
%    0
%    0
%    0
%    0
%    0
%    0
%    0
%    0];

obs_samples = [0.02246
    0.04883
    0.006836
   -0.04688
   -0.03613
    0.02344
    0.05078
    0.006836
   -0.04785
   -0.03711
   0
   0
   0
   0
   0
   0
   0
   0
   0
   0
   0
   0
   0
   0
   0
   0
   0
   0
   0
   0
   0
   0];

N = 32;
V = 32;
Fs = 5000;
y_fixed = fi(obs_samples,1,12,11);

tempfft = dsphdl.FFT;
loopCount = getLatency(tempfft,N,V)+N/V;
Yf = zeros(V,loopCount);
validOut = false(V,loopCount);
for loop = 1:1:loopCount
    if ( mod(loop,N/V) == 0 )
        i = N/V;
    else
        i = mod(loop,N/V);
    end
    [Yf(:,loop),validOut(loop)] = HDLFFT32((y_fixed(:,i)),(loop<=N/V));
end

C = Yf(:,validOut==1);
Yf_flat = C(:);

Yr =  bitrevorder(Yf_flat);
plot(Fs/2*linspace(0,1,N/2),2*abs(Yr(1:N/2)/N))
plot(2*abs(Yf_flat(1:N/2)/N))
title('Single-Sided Amplitude Spectrum of Noisy Signal y(t)')
xlabel('Frequency (Hz)')
ylabel('Output of FFT (f)')

function [yOut,validOut] = HDLFFT32(yIn,validIn)
%HDLFFT128 
% Processes one sample of FFT data using the dsphdl.FFT System object(TM)
% yIn is a fixed-point scalar or column vector. 
% validIn is a logical scalar value.
% You can generate HDL code from this function.

  persistent fft32;
  if isempty(fft32)
    fft32 = dsphdl.FFT(FFTLength=32,BitReversedOutput=false);
  end    
  [yOut,validOut] = fft32(yIn,validIn);
end