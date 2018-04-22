from sympy import symbols, trigsimp, sin, cos, pi, diff, sqrt
a1, K, L, gv, ga = symbols("a1 K L g_a g_v")

a = pi / 2 + K * (L / 2 - a1)

r1 = (1 / K + ga * cos(gv * a)) * cos(a)
r2 = (1 / K + ga * cos(gv * a)) * sin(a)

dr1 = diff(r1, a1)
dr2 = diff(r2, a1)

n1 = -dr2 / sqrt(dr1**2 + dr2**2)
n2 = dr1 / sqrt(dr1**2 + dr2**2)

dn1 = diff(n1, a1)
dn2 = diff(n2, a1)

print(trigsimp(dn2))
