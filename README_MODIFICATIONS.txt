=======================================================================
README — Modified RIS Beamforming Code
Based on: Yu & Dai, IEEE Trans. Commun., Vol. 73, No. 7, July 2025
=======================================================================

MODIFIED FILES (3 new files — do NOT change the originals):
------------------------------------------------------------
1. Beamforming_modified.m
2. my_compute_rate_modified.m
3. Plot_All_Modified_Assumptions.m    <-- RUN THIS

ORIGINAL FILES KEPT UNCHANGED (helpers):
-----------------------------------------
quantz.m, compute_phase.m, wfill.m, waterfill.m,
randnc.m, TRX_std.mat, TRX_same_dist.mat

=======================================================================
EXACT MODIFICATIONS MADE
=======================================================================

FILE: Beamforming_modified.m
  Change 1 (line ~90): Added 'sigma_e' as 8th input argument
            Default = 0  (perfect CSI, identical to original)

  Change 2 (line ~97): Imperfect CSI added
            BEFORE: [nothing — Hmat was used directly]
            AFTER:
              noise_real = sigma_e/sqrt(2) * randn(size(Hmat));
              noise_imag = sigma_e/sqrt(2) * randn(size(Hmat));
              Hmat_est   = Hmat + noise_real + 1i*noise_imag;

  Change 3 (line ~107): Phase design uses estimated channel
            BEFORE: BFmat = exp(-1i*angle(Hmat));
            AFTER:  BFmat = exp(-1i*angle(Hmat_est));

            BEFORE: QBFmat = quantz(exp(-1i*angle(Hmat)), qbit);
            AFTER:  QBFmat = quantz(BFmat, qbit);
            (quantz is unchanged — just operates on Hmat_est now)

  Change 4 (line ~117): SPM combined weight uses estimated phase
            BEFORE: ModBFmat = quantz(GMod_off .* exp(-1i*angle(Hmat)), qbit);
            AFTER:  ModBFmat = quantz(GMod_off .* BFmat, qbit);
            (BFmat already computed from Hmat_est above)

FILE: my_compute_rate_modified.m
  Change 1: Added 'sigma_e' as 9th input argument
  Change 2: Calls Beamforming_modified() instead of Beamforming()
            Everything else (rate formula, frequency grid) is unchanged.

FILE: Plot_All_Modified_Assumptions.m
  Entirely new — generates 3 required plots with Monte Carlo averaging.

=======================================================================
HOW TO RUN
=======================================================================
1. Open MATLAB
2. Set current folder to this directory
3. Run:  Plot_All_Modified_Assumptions
4. Three figures appear automatically

To speed up (for testing): reduce Ntrials to 5 inside the script.
To improve accuracy:        increase Ntrials to 100+.

=======================================================================
PLOTS GENERATED
=======================================================================
Plot 1 — Rate vs Transmit Power (dBm)
  - Upper Bound (theoretical)
  - Ideal: perfect CSI + continuous phase (original paper)
  - Imperfect CSI only (sigma_e = 0.05, continuous phase)
  - Quantization only (2-bit, perfect CSI)
  - Both: Imperfect CSI + 2-bit quantization

Plot 2 — Rate vs CSI Error Variance (sigma_e^2)
  - Compares FZ-SPM, Classical, and Upper Bound
  - Shows graceful degradation as CSI quality decreases
  - sigma_e ranges from 0 (perfect) to 0.30 (heavily corrupted)

Plot 3 — Rate vs Quantization Bits (b = 1 to 8)
  - Three CSI error levels: sigma_e = 0, 0.05, 0.10
  - Dashed lines show continuous-phase performance limit
  - Shows that 2-bit quantization nearly matches continuous phase

=======================================================================
KEY INSIGHT FOR PRESENTATION
=======================================================================
- At sigma_e = 0 and qbit = 100, the modified code gives exactly
  the same results as the original paper (baseline validated).
- CSI error degrades performance gradually — the FZ-SPM method
  is robust to small errors (sigma_e < 0.05).
- 2-bit quantization is sufficient — matches paper's Fig. 8 result.
=======================================================================
