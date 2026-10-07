"""N6/N12 허브 개수 수열의 4차 차분과 폐형식 다항식, 극한값 검증.
순수 조합론 계산(코드 실행 불필요, 뉴턴 전진차분법)."""
from fractions import Fraction as F

N6 = [56, 184, 432, 840, 1448]
N12 = [38, 116, 260, 490, 826]

def diff(seq, k):
    for _ in range(k):
        seq = [seq[i+1]-seq[i] for i in range(len(seq)-1)]
    return seq

print('N6 4차차분:', diff(N6, 4))
print('N12 4차차분:', diff(N12, 4))

N6f = [F(x) for x in N6]
d0,d1,d2,d3 = diff(N6f,0)[0], diff(N6f,1)[0], diff(N6f,2)[0], diff(N6f,3)[0]
lead6 = d3/6
N12f = [F(x) for x in N12]
e0,e1,e2,e3 = diff(N12f,0)[0], diff(N12f,1)[0], diff(N12f,2)[0], diff(N12f,3)[0]
lead12 = e3/6
print('N6 최고차항계수:', lead6, '(초안: 20/3)')
print('N12 최고차항계수:', lead12, '(초안: 10/3)')
print('극한 N6/(N6+N12):', lead6/(lead6+lead12), '(초안: 2/3)')
