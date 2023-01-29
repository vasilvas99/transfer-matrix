import numpy as np
from ising_tm_py import generate_ising_tm
from scipy.sparse.linalg import eigsh

def corr_len(n_spins, temp):
    transfer_matrix = generate_ising_tm(int(n_spins), -1.0, float(temp))
    eig_val, _ = eigsh(transfer_matrix, 2)
    eig_val = sorted(eig_val, reverse=True)[:2]
    cln = 1.0 / (np.log(eig_val[0] / abs(eig_val[1])))
    return cln

N_MIN = 8
N_MAX = 8

T_MIN = 2.0
T_MAX = 2.8
T_MAX_ITER= 200

t_step = (T_MAX-T_MIN)/T_MAX_ITER

print("t", end="")
for n in range(N_MIN, N_MAX+1):
    print(f"\t n = {n}", end="")
print("")

for t_iter in range(0, T_MAX_ITER+1):
    t = T_MIN+t_step*(t_iter )
    print(f"{t}", end="")
    for n in range(N_MIN, N_MAX+1):
        xi = corr_len(n, t)
        print(f"\t{xi}", end="")
    print("")

