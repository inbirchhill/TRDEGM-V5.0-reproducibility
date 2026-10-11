"""중력형 장거리 지수와 자유낙하 지수 점검 (V5.0 최종본2, §4.4.1 / §9.2.1 / 부록 A.33 / 부록 C)

 (1) 원본 데이터(gravity_v2_meta.npz, gravity_v2_final.npz)에 기록된 지수를
     그대로 출력한다. meta 의 각 행은 N_k 상한(32~88)별 지수이고, 범위는
     -1.026 ~ -0.995 이다. 그동안 인용해 온 -0.964는 어느 행에서도 나오지 않는다.
 (2) 원시 테이블(gravity_v2_raw.npz, N_k=16~88의 12개 값)로 메타외삽 절차를
     다시 계산해 meta 값을 재현하고, 무한격자 극한을 외삽 방식별로 비교한다.
     1/N_k: -0.975 (R2 0.958), 1/N_k^2: -0.989 (R2 0.995),
     a+b*N_k^-p: -0.992 (p~2.6). 극한값은 외삽 방식에 따라 달라진다.
 (3) 퍼텐셜 지수 alpha 를 미분하면 힘 지수는 alpha+1 이다
     (원본 데이터 0.996 -> 1.996, 인용값 0.964 -> 1.964).
     F ~ r^-2.15 는 별도 계산의 보고값으로, 그 계산 출력은 보존돼 있지 않다.
 (4) F ∝ r^-n 에서 정지 출발 자유낙하 시간은 차원분석으로 t ∝ r0^((n+1)/2).
     kepler_test.py 의 적분값 1.4831 은 입력 n=1.964 에서 나오는 1.482 와 같다
     (독립적인 Kepler 검증이 아님). 서로 다른 n 으로 직접 확인한다.

필요 파일: gravity_v2_meta.npz, gravity_v2_final.npz, gravity_v2_raw.npz (같은 폴더)
"""
import numpy as np
from scipy.stats import linregress
from scipy.optimize import curve_fit

print("(1) 원본 데이터의 N_k 상한별 지수")
meta = np.load('gravity_v2_meta.npz')['meta']
for Nk, slope, r2, _ in meta:
    print(f"  N_k 상한={int(Nk):3d}: 지수={slope:.4f}, R2={r2:.5f}")
fin = np.load('gravity_v2_final.npz')
print(f"  최종 채택값(N_k=48~88, 1/N_k 외삽): 지수={float(fin['slope']):.4f}, R2={float(fin['r2']):.4f}")
print(f"  -0.964 와의 최소 차이: {np.abs(meta[:, 1] + 0.964).min():.4f}")

print("\n(2) 원시 테이블로 메타외삽 재계산 (gravity_channel_v2 와 같은 절차)")
raw = np.load('gravity_v2_raw.npz')
Nk_list = list(raw['Nk_list']); R = raw['R_dist']; G = raw['G_table']
inv = 1.0 / np.array(Nk_list, float)
rows = []
for cut in [24, 32, 40, 48, 56, 64, 72, 80, 88]:
    idx = [i for i, n in enumerate(Nk_list) if n <= cut]
    if len(idx) < 4:
        continue
    use = idx[-4:]
    Ginf = np.zeros(len(R)); fr = np.zeros(len(R))
    for j in range(len(R)):
        r = linregress(inv[use], np.abs(G[use, j])); Ginf[j] = r.intercept; fr[j] = r.rvalue ** 2
    v = (Ginf > 0) & (fr >= 0.90)
    rf = linregress(np.log(R[v]), np.log(Ginf[v]))
    rows.append((cut, rf.slope))
rows = np.array(rows)
print(f"  계산한 N_k: {Nk_list[0]}~{Nk_list[-1]} ({len(Nk_list)}개), 메타외삽 점: N_k 상한 "
      f"{int(rows[0, 0])}~{int(rows[-1, 0])} ({len(rows)}개)")
print(f"  meta 파일과의 최대 차이: {np.abs(rows[:, 1] - meta[:, 1]).max():.2e}")
c, s = rows[:, 0], rows[:, 1]
for name, x in [('1/N_k  ', 1 / c), ('1/N_k^2', 1 / c ** 2)]:
    r = linregress(x, s)
    print(f"  {name} 외삽: 극한 지수={r.intercept:.4f}, R2={r.rvalue ** 2:.4f}")
p, cov = curve_fit(lambda n, a, b, q: a + b * n ** (-q), c, s, p0=(-1, 10, 1.5), maxfev=20000)
print(f"  a+b*N_k^-p 외삽: 극한 지수={p[0]:.4f} (p={p[2]:.2f})")
for k in (6, 5, 4):
    a1 = linregress(1 / c[-k:], s[-k:]).intercept; a2 = linregress(1 / c[-k:] ** 2, s[-k:]).intercept
    print(f"  마지막 {k}점: 1/N_k {a1:.4f}, 1/N_k^2 {a2:.4f}")

print("\n(3) 퍼텐셜 지수 -> 힘 지수 (alpha+1)")
print(f"  원본 데이터 {-float(fin['slope']):.4f} -> {1 - float(fin['slope']):.4f}, 인용값 0.964 -> 1.964"
      "  (r^-2.15 는 출력이 보존되지 않은 별도 보고값)")

print("\n(4) 자유낙하 지수 = (n+1)/2 (차원분석)")


def freefall_time(r0, n, K=5.84, dt=1e-3):
    r, v, t = r0, 0.0, 0.0
    while r > 0.1 and t < 1e4:
        v += -K / r ** n * dt
        r += v * dt
        t += dt
    return t


r0s = np.array([3.0, 5.0, 8.0, 12.0, 18.0, 25.0])
for n in [1.5, 1.964, 1.996, 2.0, 2.15]:
    ts = np.array([freefall_time(r, n) for r in r0s])
    res = linregress(np.log(r0s), np.log(ts))
    print(f"  n={n:5.3f}: 적분 지수={res.slope:.4f}, 예측 (n+1)/2={(n+1)/2:.4f}")
