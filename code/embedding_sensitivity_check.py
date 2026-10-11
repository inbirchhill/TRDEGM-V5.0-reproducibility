"""임베딩함수 F(rho)=-eps*sqrt(rho)의 민감도비(N6 밀도10 대 N12 밀도20)가
정확히 sqrt(2)인지 직접 대수 계산으로 검증. F(rho)=-eps*sqrt(rho)는
이 논문에서 이미 확립된 Finnis-Sinclair 임베딩함수(원자분자구조, §5.5)."""
import numpy as np
rho_N6, rho_N12 = 10, 20
# dF/drho = -eps/(2*sqrt(rho)), 민감도(|dF/drho|)는 1/sqrt(rho)에 비례
sens_ratio = np.sqrt(rho_N12/rho_N6)
print(f'민감도비(N12 대 N6) = sqrt(rho_N12/rho_N6) = {sens_ratio:.6f}')
print(f'sqrt(2) = {np.sqrt(2):.6f}')
print(f'일치여부: {np.isclose(sens_ratio, np.sqrt(2))}')
