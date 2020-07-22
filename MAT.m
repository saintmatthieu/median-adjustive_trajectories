% function [FF,IC,k] = MAT(C_,FF_in)

% Takes CSPE spectrum C_ (Short and Garcia, 2006) and a priori estimate
% FF_in of fundamental frequency of analysed string tone.

% Returns medians of all possible Fundamental Frequency (FF) and
% Inharmonicity Coefficient (IC) estimates upon equations (11) and
% (12) in DAFx-10 paper (Hodgkinson et al.).

% Let X_ = fft(x_(1:M)) and Y_ = fft(x_(2:M+1)). Then, C_ = conj(X_).*Y_.
% (Preferably use window before taking the FFTs.)

% ¡¡¡ Important !!!
% FF_in is expected in ¡frequency bins! E.g. for a fundamental frequency in
% herz FF_hz, a sampling rate SR and an FFT length M, then FF_in =
% M*FF_hz/SR. Likewise, FF is returned in frequency bins.

function [FF,IC,k] = MAT(C_,FF_in)

M = length(C_);

% A frequency-scaling constant, used recurrently...
M2pi = .5*M/pi;
    
% Make vector of peak indices.
mag_ = abs(C_);
Cp_ = find(mag_(1:end-2)<mag_(2:end-1)...
                      &mag_(2:end-1)>mag_(3:end))+1;

% Maximal frequency for search on 60dB SNR basis
f_max = find(mag_(1:.5*M)>max(mag_(Cp_))*1e-6, 1, 'last' )-1;

% Estimate maximal number of partials within [0,f_max]
K = ceil(f_max/FF_in);
f_ = zeros(K,1);   % Partial frequencies

% (K^2-K)/2 is number of potential FF and IC estimates
FF_ = zeros(.5*(K^2-K),1);
IC_ = FF_;

% Initiate peak search with a priori FF estimate FF_in
Cpi = closest(Cp_,1,FF_in+1);

f_(1) = M2pi*angle(C_(Cp_(Cpi)));

% Estimate of partial C_ for next index of k
fk = 2*f_(1);

for k=2:K
    Cpi = closest(Cp_,Cpi,fk+1);
    f_(k) = M2pi*angle(C_(Cp_(Cpi)));
    
    I = .5*(k^2-k);     % Current number of FF & IC estimates
    i_ = I+(1-k:-1)+1;
    l = (1:k-1)';
    
    % Refer to paper in DAFx10 for formulae
    fk2 = f_(k)^2;
    fl2 = f_(l).^2;
    k2 = k*k;
    l2 = l.*l;
    FF_(i_) = (l2.*l2*fk2-k2*k2*fl2)./(l2*k2.*(l2-k2));
    IC_(i_) = -(l2*fk2-k2*fl2)./(l2.*l2*fk2-k2*k2*fl2);
    
    % Get median ; computationally more efficient this way
    FF_(1:I) = sort(FF_(1:I));
    IC_(1:I) = sort(IC_(1:I));
    if ~mod(I,2)
        FF = mean(sqrt(FF_(.5*I+[0 1])));
        IC = mean(IC_(.5*I+[0 1]));
    else
        FF = sqrt(FF_(.5*(I+1)));
        IC = IC_(.5*(I+1));
    end
    
    % For indices k<5, median estimate might not be too reliable yet.
    % Until then, just assume harmonic model.
    if k<5, fk = (k+1)*FF;
    else fk = (k+1)*FF*sqrt(1+IC*(k+1)^2); end
    
    if fk>f_max, break, end

end

% Just to make sure that no impossible value is returned.
IC = max(IC,0);

function Cpi = closest(Cp_,Cpi,i)

P = length(Cp_);
while Cp_(Cpi)<i
    Cpi = Cpi+1;
    if Cpi>P, Cpi = Cpi-1; break, end
end
if Cpi>1, Cpi = Cpi - (i < mean(Cp_(Cpi+[-1 0]))-1); end