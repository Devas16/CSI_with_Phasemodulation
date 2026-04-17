function Rates = my_compute_rate_modified(qbit, Nsub, L, TRX, fc, BW_des, Ptx, showfig, sigma_e)
% =========================================================================
% my_compute_rate_modified.m
% MODIFIED version of my_compute_rate.m from Yu & Dai (IEEE Trans. Commun. 2025)
%
% CHANGES FROM ORIGINAL:
%   1. Added input argument: sigma_e (CSI error std dev)
%   2. Calls Beamforming_modified() instead of Beamforming()
%      — all imperfect CSI and quantization logic is inside Beamforming_modified
%   3. Rate computation loop is identical to original
%
% INPUTS:
%   qbit    - quantization bits (use large number e.g. 100 for continuous)
%   Nsub    - number of OFDM subcarriers
%   L       - RIS side length (m)
%   TRX     - [x1;y1;z1;x2;y2;z2] (6x1)
%   fc      - carrier frequency (Hz)
%   BW_des  - bandwidth (Hz)
%   Ptx     - transmit power (dBm)
%   showfig - plot intermediate figures flag
%   sigma_e - CSI noise std dev (0 = perfect CSI)
%
% OUTPUTS:
%   Rates(1) - Rate with classical beamforming (QBFmat)
%   Rates(2) - Rate with SPM-based beamforming (ModBFmat)  [FZ-SPM]
%   Rates(3) - Theoretical upper bound rate
% =========================================================================

if nargin < 9
    sigma_e = 0;    % default: perfect CSI
end

Ngrid = 1000;

x1 = TRX(1); y1 = TRX(2); z1 = TRX(3);
x2 = TRX(4); y2 = TRX(5); z2 = TRX(6);

c = 3e8;
Ptx_lin = 10.^((Ptx-30)/10);
Lb      = c/fc;
spac    = Lb/2;
hp      = 6.625e-34;
kb      = 1.38e-23;
Tk      = 290;

N = round(L/spac);
xcoord = linspace(-L/2, L/2, N);
ycoord = xcoord;
[xx, yy] = meshgrid(xcoord, ycoord);

a       = (sqrt((xx-x1).^2+(yy-y1).^2+z1^2) + sqrt((xx-x2).^2+(yy-y2).^2+z2^2))/2;
a_vec   = linspace(min(a(:)), max(a(:)), Ngrid)';
Distmat  = a*2;
distgrid = a_vec*2;
Distmat1 = sqrt((xx-x1).^2+(yy-y1).^2+z1^2);
Distmat2 = sqrt((xx-x2).^2+(yy-y2).^2+z2^2);

% =========================================================================
% Call MODIFIED Beamforming — passes sigma_e for imperfect CSI
% =========================================================================
[ModBFmat, QBFmat, modfn] = Beamforming_modified(qbit, L, TRX, fc, BW_des, Ptx, showfig, sigma_e);

% Frequency grid (identical to original)
freqvec     = 2*pi * linspace(fc-BW_des/2, fc+BW_des/2, Nsub);
antt_delays = Distmat(:)/c;
AF_product  = antt_delays * freqvec;
h_0         = (Distmat1(:)*freqvec/c) .* (Distmat2(:)*freqvec/c);

% Channel matrix (true channel — used for received power)
Wb_mch = exp(-1i*AF_product) ./ h_0;

% Received power spectra
SCpowers_SPM     = abs(ModBFmat(:).' * Wb_mch).^2;   % FZ-SPM power
SCpowers_classic = abs(QBFmat(:).'  * Wb_mch).^2;   % Classical power

% Upper bound power (flat in-band)
SCpowers_ideal = zeros(size(freqvec));
SCpowers_ideal(abs(freqvec/2/pi - fc) <= BW_des/2) = 1;
Pow = sum(modfn.^2) * spac^-4 * (a_vec(end)-a_vec(1))/Ngrid * c / (freqvec(2)-freqvec(1)) * pi ./ h_0(end/2,end/2)^2;
SCpowers_ideal_norm = SCpowers_ideal * Pow / sum(SCpowers_ideal);

% Rate computation (identical formula to original)
Delta_f  = freqvec(2) - freqvec(1);
Rates    = zeros(3, 1);
phin_nr  = hp * freqvec;
phin_dr  = exp(hp*freqvec/(kb*Tk)) - 1;
phin     = phin_nr ./ phin_dr;
powalloc = Ptx_lin / BW_des;

% Rate 1: Classical beamforming
N_eff    = phin ./ SCpowers_classic;
SNR_eff  = powalloc ./ N_eff;
Rates(1) = Delta_f * sum(log2(1 + SNR_eff));

% Rate 2: FZ-SPM beamforming
N_eff    = phin ./ SCpowers_SPM;
SNR_eff  = powalloc ./ N_eff;
Rates(2) = Delta_f * sum(log2(1 + SNR_eff));

% Rate 3: Upper bound
N_eff    = phin ./ SCpowers_ideal_norm;
SNR_eff  = powalloc ./ N_eff;
Rates(3) = Delta_f * sum(log2(1 + SNR_eff));

end
