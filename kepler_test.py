"""자유낙하(방사) 재응집시간이 케플러 t~r0^1.5 스케일링을 따르는지 검증.
힘법칙(K_GRAV=5.84, GRAV_EXP=1.964)은 part11_nbody.txt의 기확립값 재사용.
"""
import numpy as np
from scipy.stats import linregress

GRAV_EXP = 1.964
K_GRAV = 5.84

def radial_freefall_time(r0, dt=0.001):
    """중심을 향해 정지상태에서 낙하, 원점 근방(r<0.1) 도달 시간"""
    r = r0
    v = 0.0
    t = 0.0
    while r > 0.1 and t < 10000:
        F = -K_GRAV / r**GRAV_EXP  # 인력(중심방향, 음수)
        v += F*dt
        r += v*dt
        t += dt
        if r <= 0:
            break
    return t

r0_list = np.array([3.0, 5.0, 8.0, 12.0, 18.0, 25.0])
times = np.array([radial_freefall_time(r0) for r0 in r0_list])
print("r0, t:")
for r0, t in zip(r0_list, times):
    print(f"  r0={r0}: t={t:.4f}")

logr = np.log(r0_list); logt = np.log(times)
res = linregress(logr, logt)
print(f"\n스케일링: t ~ r0^{res.slope:.4f} (R^2={res.rvalue**2:.6f})")
print(f"케플러 예측(1.5)과의 차이: {abs(res.slope-1.5):.4f}")
