# V5.0 종합작업 세션 — 재현 검증 코드/데이터 전체 목록

이번 최종본종합 작업 세션(2026-09-19~20)에서 직접 재실행하여
검증한 코드와 데이터입니다. 본문에는 결과만 반영돼 있고, 아직
부록C(계산 방법)에는 파일명이 기록되지 않았습니다 — 다음 작업
때 부록C 보완이 필요합니다. GitHub 공개 전 정리용으로 미리
모아둡니다.

| 파일 | 검증한 본문 위치 | 확인한 결과 | 원본 출처 |
|---|---|---|---|
| co_convergence_sphericity.py | §4.4.3(RD와CO) | depth0-4 구형도 0.97343→0.93747, 초안과 정확히 일치 | 소실코드_복구분.zip/recovered_output/ |
| network_base.py | §4.2, 여러 곳의 공통 의존성 | S(d)/C(d) depth0-4 재현 | /mnt/project/ |
| bond_stiffness.py | §5.5(원자·분자구조) | 결합력상수, 정규화진동수 1.000:0.950:0.420, N6/N12 총강성 9.724/9.4517 | /mnt/project/ |
| bond_curve_test.txt | §5.5 | 평형거리오차 0.45%, 해리에너지순서(6-6>6-12>12-12) | /mnt/project/ |
| synthetic_bipartite_sweep.txt | §4.3.5(이종연결증폭) | 교차이웃수6에서 정점 찍는 비단조 패턴 | /mnt/project/ |
| frequency_as_temperature.txt | §4.3.5(문턱현상) | 진동수비1.1-1.2 구간 급격한 도약, √2는 포화구간 | /mnt/project/ |
| swing_eom.py | 위 두 코드의 공통 의존성 | - | /mnt/project/ |
| hexadecapole_kspace_check.py | §4.4.1(왜곡과유효중력, 육극자) | c4=K4 피팅 R²=1.00000 | basis_data.npy 필요 |
| basis_data.npy | 위 코드의 입력데이터(프로젝트 버전은 손상됨) | - | 논문코드백업.zip(정상본) |
| gravity_v2_final.npz | §4.4.1(중력지수 원본재현) | slope=-0.99578, R²=0.99982 | 논문코드백업.zip |
| gravity_v2_meta.npz | §4.4.1 | N_k=32→88 수렴시퀀스 | 논문코드백업.zip |
| gravity_v2_raw.npz | §4.4.1 | 원시 그린함수 테이블 | 논문코드백업.zip |
| rd_ft_activation.txt | §5.6(두갈래경로와간섭) | 활성화확률 13.35/21.21/0.00/65.44% | /mnt/project/ |
| analyze_hop_result.py | §6.1(초록 핵심성과, p=0.008) | 홉1=0.1907,홉3=0.1586,10/12,p=0.00806 | handoff_2a.zip |
| mid150_matched_ckpt.npz | 위 코드의 입력데이터(174허브 원자료) | - | handoff_2a.zip |
| mid150_matched_queue.json | 위 코드의 입력데이터 | - | handoff_2a.zip |
| mid150_matched_worker.py | 174허브 원자료 생성 스크립트(참고용) | - | handoff_2a.zip |
| hysteresis_check.txt | §6.1(이력현상) | K=0.10에서 forward0.5849/backward0.9383 | /mnt/project/ |
| R_function_fit.txt | §6.1.2(재규격화사상) | r=0.8316, R²=0.99868(b_n), p=1.3789, R²=0.99808((a_n-1)) | /mnt/project/ |
| reconstruct_erepr.py | §6.3.1(ER=EPR유비, **독립 재구성**, 원본 아님) | T=30 평평화(균등분포 근사), T=2-10 판단불가, 정성적 일치 | 신규 작성(supercell_builder.txt 기반) |
| finnis_sinclair.py, fs_hub_corr.py | §7.2(기하변수와에너지, **재구성시도**) | 배위수-변위상관 반대방향(r=-0.73~-0.84), 정직하게 미확정 | /mnt/project/ + 신규작성 |
| model_ii_rotation_stiffness.txt | §7.3(3체각도항) | Model I=1.6176, Model II균일가중=1.0822, 원본과 완전일치 | /mnt/project/ |
| model_ii_scan2d.txt | §7.4(음의각도항) | (0,0)=1.6176,(0.8,0)=1.0822 재확인, λ12증가시 악화되는 체계적패턴 확인 | /mnt/project/ |
| q1_trackAB_full.txt, supercell_builder.txt | §9.1.2(결어긋남, **부분검증**) | 기본설정으론 재현안됨, 노이즈 파라미터 키우면 정성적 부합(λ2.0→0.88, λ5.0→0.16~0.51), TrackA발산은 §7.5 불안정성과 일치 | /mnt/project/ + 논문코드백업.zip |
| goal1_orbit_setup.txt, part11_nbody.txt, kepler_test.py | §9.2.1(중력과기하, 재응집시간, **독립재구성**) | t~r0^1.4831(R²=1.000000), 초안의 t~r0^1.5(R²=0.999995)와 정성적·정량적 부합 | /mnt/project/ + 신규작성 |
| k1_scan_unified.py | §9.3.1("온도축"재명명) | K≈0.22에서 v_std=0.424 정점, 이후 급락, 원본과 완전일치 | /mnt/project/ |
| q1_trackAB_cross.txt | §10.5(결어긋남 사전체크, **부분검증**) | TrackA 발산(§7.5와 일치), 정량비교 후속코드는 못찾음 | /mnt/project/ |
| translate_exp_power_greensfn.txt | §10.5(베른슈타인정리, **부분검증**) | 정리 자체는 수학적 사실, 그린함수 실측 R²는 낮고 명확한 우세패턴 미확인 | /mnt/project/ |
| nbody_recombination_v2.txt | §11.3.2/§11.4(N체재응집, **모순발견·해소**) | 속박비율1.1~4.4%, 최종퍼짐이 초기의5.5~7배(흩어짐), §11.4의 기존서술과 모순 발견해 수정 | /mnt/project/ |
| arrival_time_test.py | **§5.3/§4.5/§10.5/§10.6/§11.12/부록C(9곳 전부, 대규모 신뢰성 재검토)** | 홉거리-도착시간 재구성: R²=0.0025(무관계), 초안의 R²≈0.90/0.74와 반대. propagation_recheck.py(사용자 제공 감사도구)의 "출처미확립+null result(r=-0.039,p=0.297)" 경고와 정합. 9곳 전부 Verified→UnresolvedTag로 하향 | 신규작성(hub_graph_mid150.npz 재사용) |
| mutual_info_fragmentation.txt | §11.12(주관역학개체, "통합"4회실패시도 중 1회) | K=0.06→3.6447, K=0.105→0.0450, 원본과 소수점까지 완전일치 | /mnt/project/ |
| spectrumC_multiseed.txt | §11.12(전역임계성vs국소파편화) | 감수율피크K=0.06, 파편화피크K=0.105, 불일치=False, 원본과 완전일치 | /mnt/project/ |
| depth2_core_fraction.txt | §11.12(아홉번째구성요소, core fraction) | K=0.110에서 core_frac=0.5448, 초안의 "처음엔0.545"와 정확일치 | /mnt/project/ |
| n6n12_closedform.py | §11.12(폐형식유도, 재분류사례, **순수조합론계산**) | 4차차분=0, N6최고차항=20/3, N12최고차항=10/3, 극한=2/3, 전부 초안과 정확일치 | 신규작성(순수계산, 코드검증 불필요) |
| escape_vs_return_scan.txt | §11.12(이탈복귀, 임계둔화) | K=0.105→1000스텝, K=0.108→2500스텝(임계둔화정점), K=0.100/0.110→500스텝, 원본과 정확일치 | /mnt/project/ |
| embedding_sensitivity_check.py | §11.12(N6우선편입근본원인)/§7.4(우연근접판정 업데이트, **순수대수계산**) | F(ρ)=-ε√ρ의 민감도비=√2(정확), RD마름모대각선비와 동일 근본상수임을 확정, §7.4 서술 갱신 | 신규작성(순수계산) |
| bond_stiffness_used.txt | §13(결합강도기하학적유도, U-20) | 상대비1.000:0.902:0.177, N6=9.724/N12=9.452(비율0.972), 원본과 완전일치 | /mnt/project/ |
| integration_as_multifacet_convergence.txt | §11.12(구성요소2통합, **§120번오류 발견수정**) | N6=-0.1112/0.9761, 무작위평균=0.1083/0.6244, 원본과 정확일치. 제가 §120번에서 잘못 인용한 "상호정보량등3지표"를 정확한 지표로 정정 | /mnt/project/ |
| lyapunov_chaos_test.txt | §11.12(카오스창정량적확정) | K=0.05/0.07(음수)→K≈0.085(0교차)→K=0.095-0.115(양수)→K=0.12(다시음수), 원본과 정확일치 | /mnt/project/ |
| deltaomega_exponent_verify.txt | §11.12(카오스창위치Δω인과관계, 관계식검증) | 배율0.75에서 예측K=0.0732 지점에 정확히 카오스(리아푸노프0.0902), 양옆은 안정, 원본과 정확일치 | /mnt/project/ |
| K_cycling_test.txt | §11.12(자연선택유비, K순환재점화검정) | 첫LOW국면 클러스터7개(파편화)→HIGH이후 3회순환 내내 클러스터2개 고정(재점화안됨), 원본과 정확일치 | /mnt/project/ |
| leading_process_test.txt | §11.12(주도하는주관프로세스, 다단계가설) | N12만의 고유벡터중심성-핵심소속 상관=0.892, 원본과 정확일치 | /mnt/project/ |
| fusion_binding_test.txt | §11.5(핵융합유비) | 신규껍질당결합에너지 0.0507→0.0521→0.0509(depth1정점), 원본과 정확일치 | /mnt/project/ |
| track2_3_ER_EPR.txt | §6.3.1 탐색 과정(**미사용**) | 2큐빗 단순화 버전, 초안이 말한 96허브 설정과 달라 결과 불일치 확인 후 reconstruct_erepr.py로 대체 | 논문코드백업.zip |
| a3_b7_spectral_decimation.py | §6.1.2 탐색 과정(**미사용**) | depth별 스펙트럼-무한격자 KS검정, 초안이 말한 "R함성 30.8%" 계산과 다른 방식(다른 p값 나옴) 확인 | /mnt/project/ |
| A1_omega_density.txt, A1_shuffle_check.txt | §7.6(위치-위상결합) | 국소밀도변조 0.53σ(인공물), 원본과 정확일치 | /mnt/project/ + 논문코드백업.zip |
| A2_velocity_coupling.txt | §7.6 | 이웃속도제곱변조 3.50σ(유의), 원본과 정확일치 | /mnt/project/ |
| A3_inertia_modulation.txt | §7.6 | 관성변조 0.70σ(비유의), 원본과 정확일치 | /mnt/project/ |
| A4_angle_term.txt | §7.6 | 결합각평균 1.73σ(15시드, 초기단계), 원본과 일치 | /mnt/project/ |
| A5_topological_coupling.txt | §7.6 | 배위수변조 3.13σ(유의), 원본과 정확일치 | /mnt/project/ |

## 다음 작업 시 필요한 것
1. ~~위 표 내용을 부록C(계산 방법)의 해당 소절에 파일명과 함께 기록.~~
   **완료(2026-09-19)** — 9개 코드 전부 부록C의 해당 소절(CO수렴계산,
   Shannon entropy계산, N-slit interference계산, Propagation계산,
   장거리scaling의수렴문제, Hexadecapole계산, 신규 C.12 화학결합
   정성적재현계산, 신규 C.13 이종연결증폭·문턱현상계산)에 파일명·
   수치와 함께 기록 완료.
2. 본문에서 "원본 코드"라고만 언급한 부분에 부록C 소절 참조 추가는
   추후 여유 있을 때 보완.
3. GitHub 공개 시 이 폴더(검증코드_모음/) 구조를 그대로 리포지토리
   하위 디렉토리로 사용 가능.
