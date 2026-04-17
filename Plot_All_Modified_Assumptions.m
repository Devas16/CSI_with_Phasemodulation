% =========================================================================
% Plot_All_Modified_Assumptions.m
%
% MAIN SIMULATION SCRIPT — generates all 3 required plots:
%   Plot 1 — Average Rate vs SNR/Transmit Power
%             (Ideal | Imperfect CSI | Quantization only | Both combined)
%   Plot 2 — Average Rate vs CSI error variance (sigma_e^2)
%   Plot 3 — Average Rate vs Number of quantization bits (b)
%
% Based on: Yu & Dai, "Near-Field Wideband Beamforming for RIS Based on
%           Fresnel Zones", IEEE Trans. Commun., Vol. 73, No. 7, Jul 2025.
%
% HOW TO RUN:
%   1. Place this file in the same folder as the original paper code
%   2. Also place Beamforming_modified.m and my_compute_rate_modified.m
%      in the same folder
%   3. Run: Plot_All_Modified_Assumptions
%
% PARAMETERS YOU CAN CHANGE:
%   Ntrials  - Monte Carlo trials (increase for smoother curves, ~50 is fast)
%   sigma_e  - CSI error std dev values for Plot 2
%   b_vals   - quantization bits for Plot 3
%   Ptx_vec  - transmit power range (dBm) for Plot 1
% =========================================================================

clear; clc; close all;

%% ── SHARED SYSTEM PARAMETERS (identical to original paper) ──────────────
Nsub   = 256;       % Number of OFDM subcarriers
L      = 1;         % RIS side length (m) — paper uses 1m for main results
fc     = 30e9;      % Carrier frequency: 30 GHz
BW_des = 1.5e9;     % Bandwidth: 1.5 GHz
Ptx_0  = 10;        % Fixed transmit power for Plots 2 & 3 (dBm)
Ntrials = 20;       % Monte Carlo trials (increase for smoother results)

% Fixed TX/RX positions (same as paper Fig. 7 / Fig. 8 geometry)
% BS at (6.4, 5, 14.4) m, UE at (-4.8, 5, 6.4) m
TRX = [6.4; 5; 14.4; -4.8; 5; 6.4];

fprintf('System setup: RIS = %.1fm x %.1fm, fc = %.0fGHz, BW = %.1fGHz\n',...
    L, L, fc/1e9, BW_des/1e9);
fprintf('Running %d Monte Carlo trials per point...\n\n', Ntrials);

%% ========================================================================
%  PLOT 1 — Rate vs Transmit Power (SNR)
%  Four curves: Ideal | Imperfect CSI | Quantization only | Both combined
% =========================================================================
fprintf('Computing Plot 1: Rate vs Transmit Power...\n');

Ptx_vec  = 0:5:30;      % Transmit power sweep (dBm)
qbit_inf = 100;          % Large number = effectively continuous phase
qbit_q   = 2;            % Fixed quantization bits for "quantization only" curve
sigma_csi = 0.05;        % Fixed CSI error std dev for "imperfect CSI" curve

% Pre-allocate [4 curves x length(Ptx_vec)]
R_ideal    = zeros(1, length(Ptx_vec));   % Curve 1: Ideal (no errors)
R_csi      = zeros(1, length(Ptx_vec));   % Curve 2: Imperfect CSI only
R_quant    = zeros(1, length(Ptx_vec));   % Curve 3: Quantization only
R_both     = zeros(1, length(Ptx_vec));   % Curve 4: CSI error + Quantization
R_upper    = zeros(1, length(Ptx_vec));   % Theoretical upper bound

for pi_idx = 1:length(Ptx_vec)
    ptx = Ptx_vec(pi_idx);
    
    r_ideal = zeros(Ntrials,1);
    r_csi   = zeros(Ntrials,1);
    r_quant = zeros(Ntrials,1);
    r_both  = zeros(Ntrials,1);
    r_upper = zeros(Ntrials,1);
    
    for t = 1:Ntrials
        % Curve 1: Ideal — perfect CSI, continuous phase
        tmp       = my_compute_rate_modified(qbit_inf, Nsub, L, TRX, fc, BW_des, ptx, 0, 0);
        r_ideal(t) = tmp(2);    % FZ-SPM rate
        r_upper(t) = tmp(3);    % upper bound
        
        % Curve 2: Imperfect CSI only — continuous phase, sigma_e > 0
        tmp       = my_compute_rate_modified(qbit_inf, Nsub, L, TRX, fc, BW_des, ptx, 0, sigma_csi);
        r_csi(t)   = tmp(2);
        
        % Curve 3: Quantization only — perfect CSI, b-bit phase
        tmp        = my_compute_rate_modified(qbit_q, Nsub, L, TRX, fc, BW_des, ptx, 0, 0);
        r_quant(t) = tmp(2);
        
        % Curve 4: Both — imperfect CSI + quantization
        tmp       = my_compute_rate_modified(qbit_q, Nsub, L, TRX, fc, BW_des, ptx, 0, sigma_csi);
        r_both(t)  = tmp(2);
    end
    
    R_ideal(pi_idx)  = mean(r_ideal);
    R_csi(pi_idx)    = mean(r_csi);
    R_quant(pi_idx)  = mean(r_quant);
    R_both(pi_idx)   = mean(r_both);
    R_upper(pi_idx)  = mean(r_upper);
    
    fprintf('  Ptx = %2d dBm done\n', ptx);
end

% ── Figure 1 ──────────────────────────────────────────────────────────────
figure('Position', [100 100 800 600]);
plot(Ptx_vec, R_upper/1e9,  'k--',  'LineWidth', 2.0); hold on; grid on;
plot(Ptx_vec, R_ideal/1e9,  'b-o',  'LineWidth', 2.0, 'MarkerSize', 7);
plot(Ptx_vec, R_csi/1e9,    'r-s',  'LineWidth', 2.0, 'MarkerSize', 7);
plot(Ptx_vec, R_quant/1e9,  'g-^',  'LineWidth', 2.0, 'MarkerSize', 7);
plot(Ptx_vec, R_both/1e9,   'm-d',  'LineWidth', 2.0, 'MarkerSize', 7);
hold off;

xlabel('Transmit Power (dBm)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Average Achievable Rate (Gbit/s)', 'FontSize', 14, 'FontName', 'Times New Roman');
title('Rate vs. Transmit Power under Modified Assumptions', 'FontSize', 14, 'FontName', 'Times New Roman');
legend(sprintf('Upper Bound (Ideal)'), ...
       sprintf('Ideal: Perfect CSI, Continuous Phase'), ...
       sprintf('Imperfect CSI only (\\sigma_e = %.2f)', sigma_csi), ...
       sprintf('Quantization only (%d-bit)', qbit_q), ...
       sprintf('Both: Imperfect CSI + %d-bit Quant.', qbit_q), ...
       'Location', 'northwest', 'FontSize', 12, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 13, 'FontName', 'Times New Roman');
set(get(gca,'Legend'), 'Box', 'off');

fprintf('Plot 1 done.\n\n');

%% ========================================================================
%  PLOT 2 — Rate vs CSI Error Variance (sigma_e^2)
%  Shows rate degradation as CSI quality worsens
% =========================================================================
fprintf('Computing Plot 2: Rate vs CSI Error Variance...\n');

% sigma_e values: from perfect (0) to heavily corrupted (0.3)
sigma_vec = [0, 0.01, 0.02, 0.05, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3];
var_vec   = sigma_vec.^2;   % variance for x-axis label

R_vs_sigma_SPM     = zeros(1, length(sigma_vec));
R_vs_sigma_classic = zeros(1, length(sigma_vec));
R_vs_sigma_upper   = zeros(1, length(sigma_vec));

for si = 1:length(sigma_vec)
    se = sigma_vec(si);
    
    r_spm = zeros(Ntrials, 1);
    r_cls = zeros(Ntrials, 1);
    r_ub  = zeros(Ntrials, 1);
    
    for t = 1:Ntrials
        % Continuous phase (qbit_inf), varying CSI error
        tmp      = my_compute_rate_modified(qbit_inf, Nsub, L, TRX, fc, BW_des, Ptx_0, 0, se);
        r_spm(t) = tmp(2);    % FZ-SPM
        r_cls(t) = tmp(1);    % Classical
        r_ub(t)  = tmp(3);    % Upper bound
    end
    
    R_vs_sigma_SPM(si)     = mean(r_spm);
    R_vs_sigma_classic(si) = mean(r_cls);
    R_vs_sigma_upper(si)   = mean(r_ub);
    
    fprintf('  sigma_e = %.3f (var = %.5f) done\n', se, se^2);
end

% ── Figure 2 ──────────────────────────────────────────────────────────────
figure('Position', [200 100 800 600]);
plot(var_vec, R_vs_sigma_upper/1e9,   'k--',  'LineWidth', 2.0); hold on; grid on;
plot(var_vec, R_vs_sigma_SPM/1e9,     'b-o',  'LineWidth', 2.0, 'MarkerSize', 8);
plot(var_vec, R_vs_sigma_classic/1e9, 'r-s',  'LineWidth', 2.0, 'MarkerSize', 8);
hold off;

xlabel('CSI Error Variance \sigma_e^2', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Average Achievable Rate (Gbit/s)', 'FontSize', 14, 'FontName', 'Times New Roman');
title(sprintf('Rate vs. CSI Error Variance (P_{tx} = %d dBm, Continuous Phase)', Ptx_0), ...
    'FontSize', 14, 'FontName', 'Times New Roman');
legend('Upper Bound', ...
       'FZ-SPM (Proposed)', ...
       'Classical Beamforming', ...
       'Location', 'southwest', 'FontSize', 12, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 13, 'FontName', 'Times New Roman');
set(get(gca,'Legend'), 'Box', 'off');

% Annotate the perfect-CSI point
text(var_vec(1)+0.0005, R_vs_sigma_SPM(1)/1e9 + 0.05, 'Perfect CSI', ...
    'FontSize', 11, 'FontName', 'Times New Roman', 'Color', 'blue');

fprintf('Plot 2 done.\n\n');

%% ========================================================================
%  PLOT 3 — Rate vs Number of Quantization Bits (b)
%  Shows convergence to continuous-phase performance as bits increase
% =========================================================================
fprintf('Computing Plot 3: Rate vs Quantization Bits...\n');

b_vals   = 1:8;          % 1-bit to 8-bit quantization
sigma_e_vals = [0, 0.05, 0.1];   % Perfect CSI, medium error, high error

R_vs_bits = zeros(length(sigma_e_vals), length(b_vals));

for bi = 1:length(b_vals)
    b = b_vals(bi);
    
    for si = 1:length(sigma_e_vals)
        se = sigma_e_vals(si);
        
        r_tmp = zeros(Ntrials, 1);
        for t = 1:Ntrials
            tmp       = my_compute_rate_modified(b, Nsub, L, TRX, fc, BW_des, Ptx_0, 0, se);
            r_tmp(t)  = tmp(2);   % FZ-SPM rate
        end
        R_vs_bits(si, bi) = mean(r_tmp);
    end
    fprintf('  b = %d bits done\n', b);
end

% Continuous-phase reference lines (b -> infinity limit)
R_continuous = zeros(1, length(sigma_e_vals));
for si = 1:length(sigma_e_vals)
    r_tmp = zeros(Ntrials, 1);
    for t = 1:Ntrials
        tmp       = my_compute_rate_modified(qbit_inf, Nsub, L, TRX, fc, BW_des, Ptx_0, 0, sigma_e_vals(si));
        r_tmp(t)  = tmp(2);
    end
    R_continuous(si) = mean(r_tmp);
end

% ── Figure 3 ──────────────────────────────────────────────────────────────
colors  = {'b', 'r', [0.6 0.1 0.8]};
markers = {'-o', '-s', '-^'};
labels  = {'\sigma_e = 0 (Perfect CSI)', ...
           sprintf('\\sigma_e = %.2f', sigma_e_vals(2)), ...
           sprintf('\\sigma_e = %.2f', sigma_e_vals(3))};

figure('Position', [300 100 800 600]);
hold on; grid on;

for si = 1:length(sigma_e_vals)
    plot(b_vals, R_vs_bits(si,:)/1e9, markers{si}, ...
        'Color', colors{si}, 'LineWidth', 2.0, 'MarkerSize', 8);
end

% Draw horizontal dashed lines for continuous-phase limit
for si = 1:length(sigma_e_vals)
    yline(R_continuous(si)/1e9, '--', 'Color', colors{si}, 'LineWidth', 1.2, ...
        'Label', sprintf('Continuous (\\sigma_e=%.2f)', sigma_e_vals(si)), ...
        'LabelHorizontalAlignment', 'right', 'FontSize', 10);
end

hold off;
xlabel('Number of Quantization Bits (b)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Average Achievable Rate (Gbit/s)', 'FontSize', 14, 'FontName', 'Times New Roman');
title(sprintf('Rate vs. Quantization Bits (P_{tx} = %d dBm)', Ptx_0), ...
    'FontSize', 14, 'FontName', 'Times New Roman');
legend(labels, 'Location', 'southeast', 'FontSize', 12, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 13, 'FontName', 'Times New Roman', ...
         'XTick', b_vals, 'XLim', [0.8 8.2]);
set(get(gca,'Legend'), 'Box', 'off');

fprintf('Plot 3 done.\n\n');
fprintf('All 3 plots generated successfully.\n');
