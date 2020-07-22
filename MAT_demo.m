% Input : acoustic guitar tone, open A
[x_,SR] = wavread('45.wav');

% Discard transient part
x_ = x_(1000:end);

% Coefficients for cosine window (2nd-order, continuous, cf. Nuttall 81)
a_ = [1;4/3;1/3];

% Order of the window
P = length(a_)-1;

% The MAT needs an a priori, rough estimate of the Fundamental Frequency.
FF_Hz = 110;

% Determine window length according to period length and window order
T_samps = round(SR/FF_Hz);
N = T_samps*2*(P+1);

% Synthesize window
w_ = zeros(N,1);
for p=0:P
    w_ = w_ + a_(p+1)*cos(p*2*pi/N*(-N/2:N/2-1)');
end

% For the freq. evaluation of the partials, I prefer to use the CSPE
% method (cf. Short and Garcia 2006)
% Simple, and no need for extra zero-padding. Hence FFT length M is simply
% power-of-two above N.
% Simply use a higher power-of-two for smoother plotting.
M = 2^ceil(log2(N));

% For CSPE : y_ is one-sample delayed frame.
y_ = x_(2:N+1);
x_ = x_(1:N);

X_ = fft(x_.*w_,M);
Y_ = fft(y_.*w_,M);
% CSPE spectrum
C_ = conj(X_).*Y_;

% Median-Adjustive Trajectories. In this implementation, the a priori FF is
% expected in frequency bins.
[FF,IC,K] = MAT(C_,FF_Hz/SR*M);

% Returned FF is in bins. Convert to Hz.
FF = FF/M*SR;

plot((0:M-1)/M*SR,db(1/N*X_),'k')
k_ = 1:K;
xtick = k_*FF.*sqrt(1+IC*k_.^2);
set(gca,'xtick',xtick,'ylim',[-96 0],'xlim',[0 SR/2],'xticklabel',k_)
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')