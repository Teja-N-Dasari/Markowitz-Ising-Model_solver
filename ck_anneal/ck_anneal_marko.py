import numpy as np
import itertools
import pandas as pd
from math import exp, log

def Marko_Ising(b, C, r, w=2, theta=(0.8,0.1,0.1)):
    theta1, theta2, theta3 = theta
    C = np.asarray(C, dtype=float)
    r = np.asarray(r, dtype=float)
    m = C.shape[0] #no. of assets
    n = m * w 
    pw = 1.0 / (2 ** (w - 1))
    factors = np.zeros(n, dtype=float) #all bits
    mapping = {} #(asset, bit_position)
    idx = 0
    for u in range(m):
        for k in range(w):
            factors[idx] = 2 ** k # it is not k-1, beacuse we have 0 index alredy fromm range(w)
            # print(u, 2**k)
            mapping[idx] = (u, k)
            idx += 1
    # print(mapping)
    ri_index = np.zeros(n, dtype=float)
    for i in range(n):
        u, _ = mapping[i]
        ri_index[i] = factors[i] * r[u]
    # print(ri_index)
    Q = np.zeros((n, n), dtype=float)
    q_lin = np.zeros(n, dtype=float)
    for i in range(n):
        u, _ = mapping[i]
        for j in range(i, n):
            v, _ = mapping[j]
            ci_j_idx = factors[i] * factors[j] * C[u, v]
            Q_ij = theta2 * ci_j_idx + theta3 * (b ** 2) * (pw ** 2) * (factors[i] * factors[j])
            Q[i, j] = Q_ij
            Q[j, i] = Q_ij
    for i in range(n):
        q_lin[i] = -theta1 * ri_index[i] - 2.0 * theta3 * (b ** 2) * (pw ** 2) * factors[i]
    J = (1/4) * Q.copy()
    J = 0.5 * (J + J.T)
    print(J)
    # np.fill_diagonal(J, 0.0)
    sum_Q_rows = np.sum(Q, axis=1)
    h = 0.5 * q_lin + 0.5 * sum_Q_rows
    gamma = theta3 * (b ** 2)
    delta = 0.25 * np.sum(Q) + 0.5 * np.sum(q_lin) + gamma
    return h, J, delta, mapping
    
def sk_anneal(J, h, T0=1.0,T_min=1e-3,tau=5000):    
    N = len(h)
    #initializing m randomly in {+1, -1}
    m = np.random.choice([-1.0, 1.0], size=N)
    print(m)
    # linear annealing schedule
    def temperature(t):
        return T0 - (T0 - T_min) * (t / tau)
    
    #here, J is matrix, m is updated matrix at each time t and so we are doing J@m each time for the feedback
    for t in range(tau):
        T = temperature(t)
        q = np.mean(m * m)        
        h_int = J @ m 
        h_eff = (h_int + h) - (1.0 - q) * m
        m = np.tanh(h_eff / T)
    return m

    
# ---------------- Assets data ----------------
r = np.array([
    0.021,  
    0.035, 
    0.048,
    0.055, 
    0.060, 
    0.072,  
    0.081, 
    0.095, 
    0.110,
    0.135
])

C = np.array([
    [0.040, 0.018, 0.012, 0.010, 0.009, 0.008, 0.007, 0.006, 0.005, 0.004],
    [0.018, 0.050, 0.020, 0.017, 0.015, 0.013, 0.012, 0.010, 0.009, 0.007],
    [0.012, 0.020, 0.060, 0.025, 0.022, 0.018, 0.016, 0.014, 0.012, 0.010],
    [0.010, 0.017, 0.025, 0.070, 0.028, 0.023, 0.020, 0.018, 0.015, 0.013],
    [0.009, 0.015, 0.022, 0.028, 0.080, 0.030, 0.025, 0.022, 0.019, 0.016],
    [0.008, 0.013, 0.018, 0.023, 0.030, 0.090, 0.035, 0.030, 0.026, 0.022],
    [0.007, 0.012, 0.016, 0.020, 0.025, 0.035, 0.100, 0.040, 0.035, 0.030],
    [0.006, 0.010, 0.014, 0.018, 0.022, 0.030, 0.040, 0.110, 0.045, 0.038],
    [0.005, 0.009, 0.012, 0.015, 0.019, 0.026, 0.035, 0.045, 0.120, 0.042],
    [0.004, 0.007, 0.010, 0.013, 0.016, 0.022, 0.030, 0.038, 0.042, 0.130]
])

# C should be positive defanite, strictly

b=10
w = 3 #Number of bits

h, J, delta, mapping = Marko_Ising(b, C, r, w=w, theta=(0.8,0.1,0.1))

m_final = sk_anneal(J, h, T0=1.0, T_min=1e-3, tau=5000)
print(m_final)
