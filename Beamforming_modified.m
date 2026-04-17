function [ModBFmat, QBFmat, modfn, phase_mod] = Beamforming_modified(qbit, L, TRX, fc, BW_des, Ptx, showfig, sigma_e)
% =========================================================================
% Beamforming_modified.m
% MODIFIED version of Beamforming.m from Yu & Dai (IEEE Trans. Commun. 2025)
%
% CHANGES FROM ORIGINAL:
%   1. Added input argument: sigma_e  (CSI error standard deviation)
%   2. Added imperfect CSI: h_est = h + e, where e ~ CN(0, sigma_e^2)
%   3. Phase alignment and SPM phase are computed on h_est, not h
%   4. Quantization (quantz) is unchanged — operates on estimated phases
%
% INPUTS (same as original + sigma_e):
%   qbit    - number of quantization bits (Inf = continuous)
%   L       - RIS side length (m)
%   TRX     - [x1;y1;z1;x2;y2;z2] TX and RX positions
%   fc      - carrier frequency (Hz)
%   BW_des  - system bandwidth (Hz)
%   Ptx     - transmit power (dBm)
%   showfig - flag to show intermediate figures
%   sigma_e - CSI error std dev (0 = perfect CSI, default = 0)
%
% OUTPUTS (identical to original):
%   ModBFmat  - SPM-based beamforming weight matrix (FZ-SPM)
%   QBFmat    - Classical beamforming weight matrix
%   modfn     - Fresnel zone reflective intensity v(a)
%   phase_mod - designed phase across Fresnel zones psi(a)
% =========================================================================

% --- Default: perfect CSI if sigma_e not provided ---
if nargin < 8
    sigma_e = 0;
end

Ngrid = 1000;

x1 = TRX(1); y1 = TRX(2); z1 = TRX(3);
x2 = TRX(4); y2 = TRX(5); z2 = TRX(6);
z12 = z1^2; z22 = z2^2;

c = 3e8;
Ptx_lin = 10.^((Ptx-30)/10);
Lb = c/fc;
spac = Lb/2;
hp = 6.625e-34;
kb = 1.38e-23;
Tk = 290;

N = round(L/spac);

xcoord = linspace(-L/2, L/2, N);
ycoord = xcoord;
[xx, yy] = meshgrid(xcoord, ycoord);

alpha = atan((y1-y2)/(x1-x2));
xc = (x1+x2)/2;
yc = (y1+y2)/2;
xxp = (xx-xc)*cos(alpha) + (yy-yc)*sin(alpha);
yyp = -(xx-xc)*sin(alpha) + (yy-yc)*cos(alpha);

LxiBF = [xxp(1,1), xxp(1,end), xxp(end,end), xxp(end,1)];
LyiBF = [yyp(1,1), yyp(1,end), yyp(end,end), yyp(end,1)];
Ldx = circshift(LxiBF,-1) - LxiBF;
Ldy = circshift(LyiBF,-1) - LyiBF;
LMat = [LxiBF; Ldx; LyiBF; Ldy];

u = (x1-xc)*cos(alpha) + (y1-yc)*sin(alpha);

a = (sqrt((xx-x1).^2+(yy-y1).^2+z1^2) + sqrt((xx-x2).^2+(yy-y2).^2+z2^2))/2;
a_vec = linspace(min(a(:)), max(a(:)), Ngrid)';
da_vec = a_vec(2) - a_vec(1);

b_vec = sqrt(a_vec.^2 - u^2);
x0_vec = 1/4*u*(z22-z12)./b_vec.^2;
k_vec = sqrt(1 + (x0_vec/u).^2 - 1/2*(z12+z22)./b_vec.^2);

A_e  = b_vec.^2.*Ldx.^2 + a_vec.^2.*Ldy.^2;
B_e  = b_vec.^2.*Ldx.*(LxiBF-x0_vec) + a_vec.^2.*Ldy.*LyiBF;
C_e  = b_vec.^2.*(LxiBF-x0_vec).^2 + a_vec.^2.*LyiBF.^2 - k_vec.^2.*a_vec.^2.*b_vec.^2;

t1 = (-B_e - sqrt(B_e.^2 - A_e.*C_e))./A_e;
t2 = (-B_e + sqrt(B_e.^2 - A_e.*C_e))./A_e;

isReal   = (imag(t1) == 0);
INrange1 = (t1<=1 & t1>=0) & isReal;
INrange2 = (t2<=1 & t2>=0) & isReal;
INrange  = INrange1 | INrange2;

cos_t1 = (LxiBF-x0_vec+Ldx.*(t1))./k_vec./a_vec;
sin_t1 = (LyiBF+Ldy.*t1)./k_vec./b_vec;
theta1 = mod(imag(log(cos_t1+1i*sin_t1)), 2*pi);

cos_t2 = (LxiBF-x0_vec+Ldx.*(t2))./k_vec./a_vec;
sin_t2 = (LyiBF+Ldy.*t2)./k_vec./b_vec;
theta2 = mod(imag(log(cos_t2+1i*sin_t2)), 2*pi);

dtheta1 = mod(theta1-circshift(theta2,1,2),2*pi).*INrange1.*circshift(INrange2,1,2);
dtheta2 = mod(theta1-circshift(theta2,2,2),2*pi).*INrange1.*circshift(INrange2,2,2).*circshift(~INrange,1,2);
dtheta3 = mod(theta1-circshift(theta2,3,2),2*pi).*INrange1.*circshift(INrange2,3,2).*circshift(~INrange,1,2).*circshift(~INrange,2,2);
dtheta4 = mod(theta1-circshift(theta2,4,2),2*pi).*INrange1.*circshift(INrange2,4,2).*circshift(~INrange,1,2).*circshift(~INrange,2,2).*circshift(~INrange,3,2);

t_proj  = ((x0_vec-LxiBF).*Ldx+(-LyiBF).*Ldy)./(Ldx.^2+Ldy.^2);
inRIS   = (sum(isReal,2)==0).*all(t_proj>0 & t_proj<1, 2);

modfn1  = sum(dtheta1,2)+sum(dtheta2,2)+sum(dtheta3,2)+sum(dtheta4,2)+inRIS*2*pi;
modfn1  = (a_vec.^2/2.*(-(z22-z12).^2/4./b_vec.^5+(z12+z22)./b_vec.^3)+k_vec.^2.*(b_vec.^2+1/2*u*u)./b_vec).*modfn1;

dstheta1 = (sin(theta1)-circshift(sin(theta2),1,2)).*INrange1.*circshift(INrange2,1,2);
dstheta2 = (sin(theta1)-circshift(sin(theta2),2,2)).*INrange1.*circshift(INrange2,2,2).*circshift(~INrange,1,2);
dstheta3 = (sin(theta1)-circshift(sin(theta2),3,2)).*INrange1.*circshift(INrange2,3,2).*circshift(~INrange,1,2).*circshift(~INrange,2,2);
dstheta4 = (sin(theta1)-circshift(sin(theta2),4,2)).*INrange1.*circshift(INrange2,4,2).*circshift(~INrange,1,2).*circshift(~INrange,2,2).*circshift(~INrange,3,2);
modfn2   = sum(dstheta1,2)+sum(dstheta2,2)+sum(dstheta3,2)+sum(dstheta4,2);
modfn2   = -a_vec*u.*(z22-z12)./(b_vec.^3).*k_vec.*modfn2;

dstheta1 = (sin(2*theta1)-circshift(sin(2*theta2),1,2)).*INrange1.*circshift(INrange2,1,2);
dstheta2 = (sin(2*theta1)-circshift(sin(2*theta2),2,2)).*INrange1.*circshift(INrange2,2,2).*circshift(~INrange,1,2);
dstheta3 = (sin(2*theta1)-circshift(sin(2*theta2),3,2)).*INrange1.*circshift(INrange2,3,2).*circshift(~INrange,1,2).*circshift(~INrange,2,2);
dstheta4 = (sin(2*theta1)-circshift(sin(2*theta2),4,2)).*INrange1.*circshift(INrange2,4,2).*circshift(~INrange,1,2).*circshift(~INrange,2,2).*circshift(~INrange,3,2);
modfn3   = sum(dstheta1,2)+sum(dstheta2,2)+sum(dstheta3,2)+sum(dstheta4,2);
modfn3   = -1/4.*k_vec.^2*u^2./b_vec.*modfn3;

modfn = modfn1 + modfn2 + modfn3;

Distmat  = a*2;
distgrid = a_vec*2;
Distmat1 = sqrt((xx-x1).^2+(yy-y1).^2+z1^2);
Distmat2 = sqrt((xx-x2).^2+(yy-y2).^2+z2^2);

% =========================================================================
% TRUE channel (perfect knowledge) — used for actual received signal power
% =========================================================================
Hmat = (Lb/(2*pi)) * exp(-1i*2*pi*fc*Distmat/c) ./ Distmat1 ./ Distmat2;

% =========================================================================
% MODIFICATION 1 — IMPERFECT CSI
% Add complex Gaussian noise to the channel to simulate estimation error.
% h_est = h + e,  e ~ CN(0, sigma_e^2)
% sigma_e = 0 reduces to perfect CSI (original paper behaviour).
% =========================================================================
noise_real = sigma_e/sqrt(2) * randn(size(Hmat));
noise_imag = sigma_e/sqrt(2) * randn(size(Hmat));
Hmat_est   = Hmat + noise_real + 1i*noise_imag;   % ESTIMATED channel

% =========================================================================
% MODIFICATION 2 — PHASE DESIGN USES ESTIMATED CHANNEL
% Original: BFmat = exp(-1i*angle(Hmat))
% Modified: BFmat = exp(-1i*angle(Hmat_est))
% This is the key change — phase alignment is done on imperfect estimate.
% =========================================================================
BFmat  = exp(-1i*angle(Hmat_est));    % <-- uses estimated channel
QBFmat = quantz(BFmat, qbit);         % quantize the estimated-channel phase

% SPM phase design (unchanged — operates on modfn which is geometry-based)
P         = cumsum(modfn.^2);
E         = P(end);
phase_mod = 2*pi*BW_des*cumsum(P)/E*2*(a_vec(2)-a_vec(1))/c - pi*BW_des*2*(a_vec-a_vec(1))/c;

indx_valid = ones(size(a)) > 0;
GMod_off   = exp(1i*compute_phase(Distmat, phase_mod, distgrid, indx_valid));

% =========================================================================
% MODIFICATION 3 — QUANTIZATION APPLIED TO COMBINED SPM + ESTIMATED PHASE
% ModBFmat = quantize( SPM_phase * estimated_classical_phase )
% When sigma_e=0 and qbit=Inf, this equals the original paper exactly.
% =========================================================================
ModBFmat = quantz(GMod_off .* BFmat, qbit);   % <-- BFmat from h_est
ModBFmat = ModBFmat .* indx_valid;
QBFmat   = QBFmat   .* indx_valid;

end
