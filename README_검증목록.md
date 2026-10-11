# V5.0 Verification Code and Data / V5.0 검증 코드·데이터 목록

**The Rhombic Dodecahedral Energy Grid Model V5.0 (Final 2)**
**마름모십이면체 에너지 그리드 모델 V5.0 (최종본2)**

This folder collects the code and data that were rerun to verify the results
reported in the paper, during the V5.0 consolidation session (2026-09-19 to 20)
and the Final 2 review (2026-10-08). Each row lists where in the paper the
result appears. Section numbers follow the paper's Korean and English editions,
which share the same numbering (appendices A–F).

이 폴더는 V5.0 종합작업 세션(2026-09-19~20)과 최종본2 점검(2026-10-08)에서
논문의 결과를 검증하기 위해 직접 재실행한 코드와 데이터를 모은 것입니다. 각 행에
그 결과가 반영된 본문 위치를 적었습니다. 절 번호는 국문판과 영문판이 같습니다
(부록 A–F).

## How to run / 실행 방법

- Python 3 with numpy, scipy, networkx. Run each script from inside this folder.
  Python 3와 numpy, scipy, networkx가 필요합니다. 각 스크립트는 이 폴더 안에서 실행합니다.
- Some files ending in `.txt` are Python scripts; run them as `python3 <file>.txt`.
  `.txt`로 끝나는 파일 중 일부는 Python 스크립트이며 `python3 <파일>.txt`로 실행합니다.
- `supercell_builder.py` is imported by `q1_trackAB_full.txt`, `q1_trackAB_cross.txt`
  and `trackAB_hessian_curvature_check.py`, so it is stored as `.py`.
  `supercell_builder.py`는 `q1_trackAB_full.txt`, `q1_trackAB_cross.txt`,
  `trackAB_hessian_curvature_check.py`가 불러오므로 `.py`로 저장되어 있습니다.

## Changes in Final 2 / 최종본2 변경 사항

- `supercell_builder.txt` → `supercell_builder.py` (renamed so that it can be imported).
  다른 스크립트가 불러올 수 있도록 확장자를 변경했습니다.
- Added `trackAB_hessian_curvature_check.py` and `gravity_exponent_check.py`, cited in the paper (Appendix C.49, §9.2.1).
  본문(부록 C.49, §9.2.1)에 인용된 새 검증 스크립트 2개를 추가했습니다.
- Added `build_hub_graph_mid150.py` and the regenerated `hub_graph_mid150.npz`, which
  `analyze_hop_result.py`, `arrival_time_test.py` and `mid150_matched_worker.py` require.
  `analyze_hop_result.py`, `arrival_time_test.py`, `mid150_matched_worker.py`에 필요한
  `hub_graph_mid150.npz`를 재생성하고 생성 스크립트를 추가했습니다.
- Corrected three rows: kepler_test (circular check), embedding_sensitivity_check (√2 claim),
  gravity_v2_final (−0.964 not in the data).
  세 행을 정정했습니다: kepler_test(순환 검증), embedding_sensitivity_check(√2 주장),
  gravity_v2_final(−0.964가 데이터에 없음).
- Corrected the gravity rows again: the meta-extrapolation uses N_k upper limits 32–88 (8 points) out of
  N_k = 16–88, and its infinite-lattice limit depends on the extrapolation form (−0.975 to −0.992).
  `gravity_exponent_check.py` now recomputes this from `gravity_v2_raw.npz`.
  중력 관련 행을 다시 정정했습니다: 메타외삽은 계산한 N_k=16–88 가운데 상한 32–88의 8점을 쓰며,
  무한격자 극한은 외삽 방식에 따라 −0.975~−0.992로 달라집니다. `gravity_exponent_check.py`가
  `gravity_v2_raw.npz`로 이를 다시 계산합니다.

## Verification list / 검증 목록

| File<br>파일 | Section<br>본문 위치 | Result<br>확인한 결과 | Source<br>원본 출처 |
|---|---|---|---|
| co_convergence_sphericity.py | §4.4.3 (RD and CO)<br>§4.4.3(RD와 CO) | Sphericity at depth 0–4: 0.97343 → 0.93747, exactly as in the paper<br>depth 0–4 구형도 0.97343→0.93747, 본문과 정확히 일치 | Recovered lost code (소실코드_복구분.zip)<br>소실코드_복구분.zip/recovered_output/ |
| network_base.py | §4.2; shared dependency<br>§4.2, 여러 곳의 공통 의존성 | Reproduces S(d)/C(d) for depth 0–4<br>S(d)/C(d) depth 0–4 재현 | Project archive<br>프로젝트 보관본 |
| bond_stiffness.py | §5.5 (Atomic·Molecular Structure)<br>§5.5(원자·분자 구조) | Bond force constants; normalized frequencies 1.000:0.950:0.420; N6/N12 total stiffness 9.724/9.4517<br>결합력상수, 정규화진동수 1.000:0.950:0.420, N6/N12 총강성 9.724/9.4517 | Project archive<br>프로젝트 보관본 |
| bond_curve_test.txt | §5.5 | Equilibrium-distance error 0.45%; dissociation-energy order 6-6 > 6-12 > 12-12<br>평형거리 오차 0.45%, 해리에너지 순서(6-6>6-12>12-12) | Project archive<br>프로젝트 보관본 |
| synthetic_bipartite_sweep.txt | §4.3.5 (heterogeneous-connection amplification)<br>§4.3.5(이종연결 증폭) | Non-monotonic pattern peaking at 6 cross-neighbors<br>교차이웃수 6에서 정점을 찍는 비단조 패턴 | Project archive<br>프로젝트 보관본 |
| frequency_as_temperature.txt | §4.3.5 (threshold phenomena)<br>§4.3.5(문턱현상) | Sharp jump at frequency ratio 1.1–1.2; √2 lies in the saturated range<br>진동수비 1.1–1.2 구간 급격한 도약, √2는 포화구간 | Project archive<br>프로젝트 보관본 |
| swing_eom.py | Shared dependency of the two scripts above<br>위 두 코드의 공통 의존성 | – | Project archive<br>프로젝트 보관본 |
| hexadecapole_kspace_check.py | §4.4.1 (Distortion and Effective Gravity; hexadecapole)<br>§4.4.1(왜곡과 유효 중력, 육극자) | c4 = K4 fit, R² = 1.00000<br>c4=K4 피팅 R²=1.00000 | Requires basis_data.npy<br>basis_data.npy 필요 |
| basis_data.npy | Input data for the script above<br>위 코드의 입력 데이터 | – (the project copy was corrupted; this is the intact one)<br>–(프로젝트 버전은 손상됨, 이것이 정상본) | Paper code backup (논문코드백업.zip)<br>논문코드백업.zip |
| gravity_v2_final.npz | §4.4.1 (gravity exponent from original data)<br>§4.4.1(중력지수 원본 재현) | slope = −0.99578, R² = 0.99982. The cited value −0.964 does not appear at any N_k in the original data (range −1.026 to −0.995; Final 2 correction)<br>slope=−0.99578, R²=0.99982. 본문 인용값 −0.964는 원본 데이터 어느 N_k에서도 나오지 않음(범위 −1.026~−0.995, 최종본2 정정) | Paper code backup (논문코드백업.zip)<br>논문코드백업.zip |
| gravity_v2_meta.npz | §4.4.1 | Exponent per N_k upper limit 32 → 88 (8 points). Infinite-lattice limit depends on the extrapolation form: 1/N_k −0.975 (R² 0.958), 1/N_k² −0.989 (R² 0.995), free power law −0.992<br>N_k 상한 32→88별 지수(8점). 무한격자 극한은 외삽 방식에 따라 다름: 1/N_k −0.975(R² 0.958), 1/N_k² −0.989(R² 0.995), 자유 거듭제곱 −0.992 | Paper code backup (논문코드백업.zip)<br>논문코드백업.zip |
| gravity_v2_raw.npz | §4.4.1 | Raw Green's-function table, N_k = 16–88 (12 values)<br>원시 Green 함수 테이블, N_k=16–88(12개 값) | Paper code backup (논문코드백업.zip)<br>논문코드백업.zip |
| rd_ft_activation.txt | §5.6 (Two-Path Interference)<br>§5.6(두 갈래 경로와 간섭) | Activation probabilities 13.35 / 21.21 / 0.00 / 65.44%<br>활성화확률 13.35/21.21/0.00/65.44% | Project archive<br>프로젝트 보관본 |
| analyze_hop_result.py | §6.1.1 (representative methodological result, p = 0.008)<br>§6.1.1(대표 방법론적 성과, p=0.008) | hop1 = 0.1907, hop3 = 0.1586, 10/12, p = 0.00806 (requires hub_graph_mid150.npz)<br>홉1=0.1907, 홉3=0.1586, 10/12, p=0.00806 (hub_graph_mid150.npz 필요) | handoff_2a.zip |
| mid150_matched_ckpt.npz | Input data for the script above (174-hub raw data)<br>위 코드의 입력 데이터(174허브 원자료) | – | handoff_2a.zip |
| mid150_matched_queue.json | Input data for the script above<br>위 코드의 입력 데이터 | – | handoff_2a.zip |
| mid150_matched_worker.py | Script that generated the 174-hub raw data (reference)<br>174허브 원자료 생성 스크립트(참고용) | – | handoff_2a.zip |
| build_hub_graph_mid150.py, hub_graph_mid150.npz | Input graph for the three scripts above and arrival_time_test.py (**new in Final 2**)<br>위 세 코드와 arrival_time_test.py의 입력 그래프(**최종본2 신규**) | The original file was lost and is regenerated from network_base.py (depth = 1, buffer_radius = 5, r² ≤ 8.75, 174 hubs, ω = 2π·(1, 1/√2)). Source selection matches the original queue exactly, one no-kick trial matches the checkpoint with zero difference, and p = 0.00806 is reproduced<br>원본 파일이 유실되어 network_base.py로 재생성(depth=1, buffer_radius=5, r²≤8.75, 174허브, ω=2π·(1, 1/√2)). 소스 선택이 원본 큐와 정확히 일치하고, 무킥 시행 하나가 체크포인트와 차이 0이며, p=0.00806이 그대로 재현됨 | Newly written<br>신규 작성 |
| hysteresis_check.txt | §6.1 (hysteresis)<br>§6.1(이력현상) | At K = 0.10: forward 0.5849 / backward 0.9383<br>K=0.10에서 forward 0.5849/backward 0.9383 | Project archive<br>프로젝트 보관본 |
| R_function_fit.txt | §6.1.2 (renormalization map)<br>§6.1.2(재규격화 사상) | r = 0.8316, R² = 0.99868 (b_n); p = 1.3789, R² = 0.99808 (a_n − 1)<br>r=0.8316, R²=0.99868(b_n), p=1.3789, R²=0.99808((a_n−1)) | Project archive<br>프로젝트 보관본 |
| reconstruct_erepr.py | §6.3.1 (ER=EPR analogy; **independent reconstruction**, not the original)<br>§6.3.1(ER=EPR 유비, **독립 재구성**, 원본 아님) | Flattening at T = 30 (near-uniform); T = 2–10 undecidable; qualitative agreement<br>T=30 평평화(균등분포 근사), T=2–10 판단불가, 정성적 일치 | Newly written (based on supercell_builder.py)<br>신규 작성(supercell_builder.py 기반) |
| finnis_sinclair.py, fs_hub_corr.py | §7.2 (Geometric Variables and Energy; **reconstruction attempt**)<br>§7.2(기하 변수와 에너지, **재구성 시도**) | Coordination–displacement correlation in the opposite direction (r = −0.73 to −0.84); left undetermined<br>배위수–변위 상관 반대 방향(r=−0.73~−0.84), 미확정으로 남김 | Project archive + Newly written<br>프로젝트 보관본 + 신규 작성 |
| model_ii_rotation_stiffness.txt | §7.3 (Three-Body Angular Term)<br>§7.3(3체 각도항) | Model I = 1.6176, Model II uniform weight = 1.0822; exact match with original<br>Model I=1.6176, Model II 균일가중=1.0822, 원본과 완전 일치 | Project archive<br>프로젝트 보관본 |
| model_ii_scan2d.txt | §7.4 (Negative Angular Terms)<br>§7.4(음의 각도항) | (0,0) = 1.6176, (0.8,0) = 1.0822 reconfirmed; systematic worsening as λ12 increases<br>(0,0)=1.6176, (0.8,0)=1.0822 재확인, λ12 증가 시 악화되는 체계적 패턴 | Project archive<br>프로젝트 보관본 |
| q1_trackAB_full.txt, supercell_builder.py | §9.1.2 (decoherence; **partial verification**)<br>§9.1.2(결어긋남, **부분 검증**) | Not reproduced at default settings; qualitative agreement with larger noise parameters (λ2.0 → 0.88, λ5.0 → 0.16–0.51); Track A divergence consistent with §7.5<br>기본 설정으론 재현 안 됨, 노이즈 파라미터를 키우면 정성적 부합(λ2.0→0.88, λ5.0→0.16~0.51), Track A 발산은 §7.5 불안정성과 일치 | Project archive + Paper code backup (논문코드백업.zip)<br>프로젝트 보관본 + 논문코드백업.zip |
| goal1_orbit_setup.txt, part11_nbody.txt, kepler_test.py | §9.2.1 (Gravity and Geometry; re-collapse time; **Final 2 correction**)<br>§9.2.1(중력과 기하, 재응집시간, **최종본2 정정**) | t ~ r0^1.4831 (R² = 1.000000). For F ∝ r^−n the free-fall exponent is (n+1)/2 by dimensional analysis, so 1.482 from the input n = 1.964 is not an independent check; the paper's t ~ r0^1.5 is inconsistent with n = 1.964 (see gravity_exponent_check.py)<br>t~r0^1.4831(R²=1.000000). F∝r^−n에서 자유낙하 지수는 차원분석으로 (n+1)/2이므로 입력 n=1.964에서 1.482가 나오는 것은 독립 검증이 아님. 본문의 t~r0^1.5는 n=1.964와 맞지 않음(gravity_exponent_check.py 참고) | Project archive + Newly written<br>프로젝트 보관본 + 신규 작성 |
| k1_scan_unified.py | §9.3.1 (renaming of the "temperature axis")<br>§9.3.1("온도축" 재명명) | Peak v_std = 0.424 at K ≈ 0.22, then a sharp drop; exact match with original<br>K≈0.22에서 v_std=0.424 정점 후 급락, 원본과 완전 일치 | Project archive<br>프로젝트 보관본 |
| q1_trackAB_cross.txt | §10.5 (decoherence pre-check; **partial verification**)<br>§10.5(결어긋남 사전 체크, **부분 검증**) | Track A diverges (consistent with §7.5); the follow-up quantitative comparison code was not found<br>Track A 발산(§7.5와 일치), 정량 비교 후속 코드는 찾지 못함 | Project archive<br>프로젝트 보관본 |
| translate_exp_power_greensfn.txt | §10.5 (Bernstein's theorem; **partial verification**)<br>§10.5(베른슈타인 정리, **부분 검증**) | The theorem itself is a mathematical fact; the measured Green's-function R² is low and no clearly dominant pattern was confirmed<br>정리 자체는 수학적 사실, 그린함수 실측 R²는 낮고 뚜렷한 우세 패턴 미확인 | Project archive<br>프로젝트 보관본 |
| nbody_recombination_v2.txt | §11.3.2 / §11.4 (N-body re-collapse; **contradiction found and resolved**)<br>§11.3.2/§11.4(N체 재응집, **모순 발견·해소**) | Bound fraction 1.1–4.4%; final spread 5.5–7× the initial (dispersal); contradiction with the earlier §11.4 text found and corrected<br>속박비율 1.1~4.4%, 최종 퍼짐이 초기의 5.5~7배(흩어짐), §11.4 기존 서술과의 모순을 발견해 수정 | Project archive<br>프로젝트 보관본 |
| arrival_time_test.py | **§5.3 / §4.5 / §10.5 / §10.6 / §11.12 / Appendix C (all 9 places; large-scale reliability review)**<br>**§5.3/§4.5/§10.5/§10.6/§11.12/부록 C(9곳 전부, 대규모 신뢰성 재검토)** | Hop distance vs. arrival time reconstructed: R² = 0.0025 (no relation), opposite to the paper's R² ≈ 0.90/0.74. Consistent with the "provenance not established + null result (r = −0.039, p = 0.297)" warning of propagation_recheck.py; all 9 places downgraded from verified to open<br>홉거리–도착시간 재구성: R²=0.0025(무관계), 본문의 R²≈0.90/0.74와 반대. propagation_recheck.py(출처 감사 도구)의 "출처 미확립 + null result(r=−0.039, p=0.297)" 경고와 정합. 9곳 전부 검증→미해결로 하향 | Newly written (requires hub_graph_mid150.npz)<br>신규 작성(hub_graph_mid150.npz 필요) |
| mutual_info_fragmentation.txt | §11.12 (Subjective Dynamics Entity; one of 4 failed "integration" attempts)<br>§11.12(주관 역학 개체, "통합" 실패 시도 4회 중 1회) | K = 0.06 → 3.6447, K = 0.105 → 0.0450; exact match with original to the decimal<br>K=0.06→3.6447, K=0.105→0.0450, 원본과 소수점까지 완전 일치 | Project archive<br>프로젝트 보관본 |
| spectrumC_multiseed.txt | §11.12 (global criticality vs. local fragmentation)<br>§11.12(전역 임계성 vs 국소 파편화) | Susceptibility peak K = 0.06, fragmentation peak K = 0.105, mismatch = False; exact match with original<br>감수율 피크 K=0.06, 파편화 피크 K=0.105, 불일치=False, 원본과 완전 일치 | Project archive<br>프로젝트 보관본 |
| depth2_core_fraction.txt | §11.12 (ninth component, core fraction)<br>§11.12(아홉 번째 구성요소, core fraction) | core_frac = 0.5448 at K = 0.110, exactly the paper's "initially 0.545"<br>K=0.110에서 core_frac=0.5448, 본문의 "처음엔 0.545"와 정확히 일치 | Project archive<br>프로젝트 보관본 |
| n6n12_closedform.py | §11.12 (closed-form derivation, reclassification; **pure combinatorics**)<br>§11.12(폐형식 유도, 재분류 사례, **순수 조합론 계산**) | 4th difference = 0; leading terms N6 = 20/3, N12 = 10/3; limit = 2/3; all exactly as in the paper<br>4차 차분=0, N6 최고차항=20/3, N12 최고차항=10/3, 극한=2/3, 전부 본문과 정확히 일치 | Newly written (pure computation)<br>신규 작성(순수 계산) |
| escape_vs_return_scan.txt | §11.12 (escape and return; critical slowing)<br>§11.12(이탈·복귀, 임계 둔화) | K = 0.105 → 1000 steps, K = 0.108 → 2500 steps (critical-slowing peak), K = 0.100/0.110 → 500 steps; exact match with original<br>K=0.105→1000스텝, K=0.108→2500스텝(임계 둔화 정점), K=0.100/0.110→500스텝, 원본과 정확히 일치 | Project archive<br>프로젝트 보관본 |
| embedding_sensitivity_check.py | §11.12.3 (root cause of N6 preferential incorporation) / §7.4 (**pure algebra; Final 2 correction**)<br>§11.12.3(N6 우선편입 근본원인)/§7.4(**순수 대수 계산, 최종본2 정정**) | Sensitivity ratio of F(ρ) = −ε√ρ is exactly √2. It equals the RD rhombic diagonal ratio numerically, but comes from the coordination ratio 20/10 and the choice of embedding function, so it is not confirmed as "the same fundamental constant"<br>F(ρ)=−ε√ρ의 민감도비=√2(정확). RD 마름모 대각선비와 수치는 같지만, 배위수비 20/10과 매립함수 선택에서 나오는 값이므로 '같은 근본상수'로 확정하지 않음 | Newly written (pure computation)<br>신규 작성(순수 계산) |
| bond_stiffness_used.txt | §13.20 (geometric derivation of coupling strength)<br>§13.20(결합 강도의 기하학적 유도) | Relative ratios 1.000:0.902:0.177; N6 = 9.724 / N12 = 9.452 (ratio 0.972); exact match with original<br>상대비 1.000:0.902:0.177, N6=9.724/N12=9.452(비율 0.972), 원본과 완전 일치 | Project archive<br>프로젝트 보관본 |
| integration_as_multifacet_convergence.txt | §11.12 (component 2, integration; **citation error found and corrected**)<br>§11.12(구성요소 2 통합, **인용 오류 발견·수정**) | N6 = −0.1112/0.9761, random mean = 0.1083/0.6244; exact match with original. The metric cited in the paper was corrected to the one actually computed<br>N6=−0.1112/0.9761, 무작위 평균=0.1083/0.6244, 원본과 정확히 일치. 본문에서 잘못 인용된 지표를 실제 계산 지표로 정정 | Project archive<br>프로젝트 보관본 |
| lyapunov_chaos_test.txt | §11.12 (quantitative confirmation of the chaos window)<br>§11.12(카오스 창 정량 확정) | K = 0.05/0.07 (negative) → K ≈ 0.085 (zero crossing) → K = 0.095–0.115 (positive) → K = 0.12 (negative again); exact match with original<br>K=0.05/0.07(음수)→K≈0.085(0 교차)→K=0.095–0.115(양수)→K=0.12(다시 음수), 원본과 정확히 일치 | Project archive<br>프로젝트 보관본 |
| deltaomega_exponent_verify.txt | §11.12 (chaos-window position vs. Δω; relation check)<br>§11.12(카오스 창 위치와 Δω의 인과관계, 관계식 검증) | At scale 0.75, chaos exactly at the predicted K = 0.0732 (Lyapunov 0.0902), stable on both sides; exact match with original<br>배율 0.75에서 예측 K=0.0732 지점에 정확히 카오스(리아푸노프 0.0902), 양옆은 안정, 원본과 정확히 일치 | Project archive<br>프로젝트 보관본 |
| K_cycling_test.txt | §11.12 (natural-selection analogy; K-cycling re-ignition test)<br>§11.12(자연선택 유비, K 순환 재점화 검정) | 7 clusters in the first LOW phase (fragmentation) → fixed at 2 clusters through 3 cycles after HIGH (no re-ignition); exact match with original<br>첫 LOW 국면 클러스터 7개(파편화)→HIGH 이후 3회 순환 내내 클러스터 2개 고정(재점화 안 됨), 원본과 정확히 일치 | Project archive<br>프로젝트 보관본 |
| leading_process_test.txt | §11.12 (leading subjective process; multi-stage hypothesis)<br>§11.12(주도하는 주관 프로세스, 다단계 가설) | Eigenvector-centrality–core-membership correlation for N12 only = 0.892; exact match with original<br>N12만의 고유벡터중심성–핵심소속 상관=0.892, 원본과 정확히 일치 | Project archive<br>프로젝트 보관본 |
| fusion_binding_test.txt | §11.5 (nuclear-fusion analogy)<br>§11.5(핵융합 유비) | Binding energy per new shell 0.0507 → 0.0521 → 0.0509 (peak at depth 1); exact match with original<br>신규 껍질당 결합에너지 0.0507→0.0521→0.0509(depth 1 정점), 원본과 정확히 일치 | Project archive<br>프로젝트 보관본 |
| track2_3_ER_EPR.txt | §6.3.1 exploration (**not used**)<br>§6.3.1 탐색 과정(**미사용**) | 2-qubit simplified version; differs from the paper's 96-hub setting, results did not match; replaced by reconstruct_erepr.py<br>2큐빗 단순화 버전, 본문의 96허브 설정과 달라 결과 불일치 확인 후 reconstruct_erepr.py로 대체 | Paper code backup (논문코드백업.zip)<br>논문코드백업.zip |
| a3_b7_spectral_decimation.py | §6.1.2 exploration (**not used**)<br>§6.1.2 탐색 과정(**미사용**) | Per-depth spectrum vs. infinite-lattice KS test; a different method from the paper's "30.8% RG composition" computation (gives different p values)<br>depth별 스펙트럼–무한격자 KS 검정, 본문의 "R 합성 30.8%" 계산과 다른 방식(다른 p값)임을 확인 | Project archive<br>프로젝트 보관본 |
| A1_omega_density.txt, A1_shuffle_check.txt | §7.6 (position–phase coupling)<br>§7.6(위치–위상 결합) | Local-density modulation 0.53σ (artifact); exact match with original<br>국소밀도 변조 0.53σ(인공물), 원본과 정확히 일치 | Project archive + Paper code backup (논문코드백업.zip)<br>프로젝트 보관본 + 논문코드백업.zip |
| A2_velocity_coupling.txt | §7.6 | Neighbor squared-velocity modulation 3.50σ (significant); exact match<br>이웃 속도제곱 변조 3.50σ(유의), 원본과 정확히 일치 | Project archive<br>프로젝트 보관본 |
| A3_inertia_modulation.txt | §7.6 | Inertia modulation 0.70σ (not significant); exact match<br>관성 변조 0.70σ(비유의), 원본과 정확히 일치 | Project archive<br>프로젝트 보관본 |
| A4_angle_term.txt | §7.6 | Mean bond angle 1.73σ (15 seeds, early stage); matches original<br>결합각 평균 1.73σ(15시드, 초기 단계), 원본과 일치 | Project archive<br>프로젝트 보관본 |
| A5_topological_coupling.txt | §7.6 | Coordination-number modulation 3.13σ (significant); exact match<br>배위수 변조 3.13σ(유의), 원본과 정확히 일치 | Project archive<br>프로젝트 보관본 |
| trackAB_hessian_curvature_check.py | §7.5 (stability-cause correction) / §13.3 / Appendix C.49 (**new in Final 2**)<br>§7.5(안정성 원인 정정)/§13.3/부록 C.49(**최종본2 신규**) | E_FS-only Hessian: analytic rank 86, 202 zero modes (cancellation at p = 2q); Track A min eigenvalue −4597.79 with 285 negative modes; Track B min positive eigenvalue 2.0076; in 6 sign tests of K′ and K″, negative modes appear only when K″ > 0<br>E_FS 단독 헤시안 해석적 랭크 86·영모드 202개(p=2q 상쇄), Track A 최소 고유값 −4597.79·음의 모드 285개, Track B 최소 양의 고유값 2.0076, K′·K″ 부호 시험 6종에서 음의 모드는 K″>0일 때만 발생 | Newly written (requires supercell_builder.py, basis_data.npy)<br>신규 작성(supercell_builder.py, basis_data.npy 필요) |
| gravity_exponent_check.py | §4.4.1 / §9.2.1 / Appendix A.33 (**new in Final 2**)<br>§4.4.1/§9.2.1/부록 A.33(**최종본2 신규**) | Exponent per N_k upper limit in the original data: −1.026 to −0.995 (no −0.964); meta-extrapolation recomputed from the raw table and compared across extrapolation forms (−0.975 to −0.992); force exponent from the potential ≈ 1.996 (r^−2.15 output not preserved); free-fall exponent (n+1)/2 confirmed for n = 1.5 / 1.964 / 1.996 / 2 / 2.15<br>원본 데이터 N_k 상한별 지수 −1.026~−0.995(−0.964 없음), 원시 테이블로 메타외삽 재계산 및 외삽 방식별 비교(−0.975~−0.992), 퍼텐셜에서 얻는 힘 지수 ≈1.996(r^−2.15 출력 미보존), 자유낙하 지수 (n+1)/2 확인(n=1.5/1.964/1.996/2/2.15) | Newly written (requires gravity_v2_*.npz)<br>신규 작성(gravity_v2_*.npz 필요) |

## Files cited in the paper but not included / 본문에 인용되었으나 포함되지 않은 파일

These files are named in the paper but their originals were not available when
this folder was assembled.
아래 파일은 본문에 이름이 나오지만, 이 폴더를 정리할 때 원본을 확보하지 못했습니다.

- **Appendix C, "The V4 Computational Lineage" / 부록 C "V4 계산 계보"**:
  `angle_data_builder.py`, `build_basis_data.py`, `orbital_hexadecapole_fibonacci.py`,
  `rd_fcc_voronoi_check.py`, `renormalization_R_function.py`, `swing_eom_dynamics.py`,
  `track_a_stability_jacobian.py`, `track_ab_stability_jacobian.py`, `weyl_law_resolution.py`
- **Provenance-audit tools / 출처 감사 도구** (§5.3.1): `propagation_recheck.py`, `PROVENANCE_AUDIT.md`
- **Follow-up verification pack / 후속 검증 묶음** (Appendix C): `06_rg_recursive_check.py`

## License / 라이선스

See `LICENSE.md` (CC BY-NC-SA 4.0, the same license as the paper).
`LICENSE.md`를 참고하세요(논문 본문과 같은 CC BY-NC-SA 4.0).
