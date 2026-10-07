# V5.0 Verification and Reconstruction Code (2025-09-06)

This folder collects code written while preparing TRDEGM V5.0
(`V5_0_초안.tex`). It covers results that were recorded in the master
log but whose original scripts had been lost. Where the original script
could be found in the project backup, it was re-run. Where it could not,
the code was re-implemented from the written description of the method
and used for independent verification. The tables below list which
claim of V5.0 each file verifies, and whether the match is complete or
partial.

**Classification criteria**: "Complete" means the result matches the
original down to the decimal places. "Partial" means the qualitative
direction or pattern is right but the exact numbers differ; V5.0 itself
honestly downgrades these with `\PartReflected` or `\UnresolvedTag`.
"Failed (not reflected)" means the result could not be reproduced and
was not carried into the V5.0 text at all.

## Complete match (identical to the original to the decimal place; reflected in V5.0 as `\Verified`)

| File | What it verifies (location in V5.0) | Key result |
|---|---|---|
| `analyze_hop_result.py` | §6.1, **the key result cited directly in the abstract** | p=0.00806 (document: p=0.008) |
| `angle_data_builder.py` | §7.5 preprocessing (generates the input data for track_ab) | angle_data.npz generated |
| `track_ab_stability_jacobian.py` | §7.5 Track A/B stability | unstable modes 0/768 over the full range K0=0.05–8.0 |
| `A1_omega_density.py` + `A1_shuffle_check.py` | §7.6 codec 1 (local density) | 0.53σ (consistent with the artifact conclusion) |
| `A2_velocity_coupling.py` | §7.6 codec 2 (velocity squared) | 3.50σ |
| `A3_inertia_modulation.py` | §7.6 codec 3 (inertia) | 0.70σ (consistent with the non-significance conclusion) |
| `A4_recheck.py` (+ `A4_angle_term_원본.txt`) | §7.6 codec 4 (mean bond angle) | delta=+0.08914±0.03187 |
| `A5_topological_coupling.py` | §7.6 codec 5 (coordination number) | 3.13σ |
| `model_ii_rotation_stiffness.py` | §7.3 rotational stiffness ratio | Model I=1.6176, Model II=1.0822 |
| `weyl_law_resolution.py` | §4.3 Weyl law | ratio 1.000 (over the full range k=0.5–3.0) |
| `hysteresis_check.py` | §3.3 hysteresis / first-order phase transition | maximum forward/backward difference 0.35 at K=0.08–0.11 |
| `mutual_info_fragmentation.py` | §11.12 mutual information vs. fragmentation | maximum at low K (0.06), minimum in the fragmentation region |

## Partial match (supports only the qualitative direction; honestly downgraded in V5.0 with `\PartReflected` / `\UnresolvedTag`)

| File | What it verifies | What matches | What does not match |
|---|---|---|---|
| `co_convergence_sphericity.py` | §4.1 convergence to CO | monotonic approach toward the target value | precise convergence rate ("halved each time depth doubles") |
| `susceptibility_scan.py` | §11.2 susceptibility plateau | a "plateau, not a narrow point" shape | exact range (0.05–0.07 vs 0.06–0.08) |
| `fusion_binding_test.py` | §11.5 nuclear-fusion binding energy | peak pattern at depth=1 | absolute numerical scale |
| `hexadecapole_kspace_check.py` | §5.2 hexadecapole symmetry | hexadecapole is nonzero, R²=1.00000 | the document's "8×" figure (a different metric, so not comparable) |

## Failed; not reflected in the V5.0 text at all (recorded honestly only in the direction-summary document)

| File | What was attempted | Why it failed |
|---|---|---|
| `model_iii_666only.py` | §7.3 666-only stage 3 (0.732) | the triplet classification criterion differs from the original (the definition of the 6-12-12 group does not match); possibly a coincidence, so not reflected |
| `dla_reconstruct.py` | §11.12 DLA fractal dimension | off-lattice DLA is sensitive to parameters; D_f does not match even with a standard implementation |
| `probabilistic_life_test.py` | §11.8 threshold-type curve for the biotechnology analogy | all 0% because of a difference in the definition of the local order variable |
| `task3_multihop_robust.py` | §6.1 supporting qualitative evidence | kept only as supplementary material, since analyze_hop_result.py already verifies the result completely |

## Data integrity

Six `.npy` binary data files in the project were found to have been
corrupted during a text-encoding conversion. Data regenerated from
scratch with `build_basis.py` was confirmed to be **completely
identical** to the intact originals that the user later located on
their desktop, which finally settled the integrity of the data.

## Usage

Most scripts require shared infrastructure modules such as
`network_base.py` and `swing_eom.py`. These are already in the project
repository (`/mnt/project` or the GitHub repository
`rd-energy-grid-model`), so they are not included here. See the import
statements at the top of each script.
