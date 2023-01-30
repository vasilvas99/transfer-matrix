import numpy as np
from ising_tm_py import generate_ising_tm
from scipy.sparse.linalg import eigsh
import argparse
import math 
from contextlib import redirect_stdout
import os

def corr_len(n_spins, temp):
    transfer_matrix = generate_ising_tm(int(n_spins), -1.0, float(temp))
    eig_val, _ = eigsh(transfer_matrix, 2)
    eig_val = sorted(eig_val, reverse=True)[:2]
    cln = 1.0 / (np.log(eig_val[0] / abs(eig_val[1])))
    return cln


parser = argparse.ArgumentParser()
parser.add_argument("--N_MIN", required=True, type=int)
parser.add_argument("--N_MAX", required=True, type=int)

parser.add_argument("--T_MIN", required=True, type=float)
parser.add_argument("--T_MAX", required=True, type=float)
parser.add_argument("--T_MAX_ITER", required=True, type=int)

parser.add_argument("--ARRAY_POS", required=True, type=int)
parser.add_argument("--ARRAY_MAX", required=True, type=int)

parser = parser.parse_args()

N_MIN = parser.N_MIN
N_MAX = parser.N_MAX


T_MIN = parser.T_MIN
T_MAX = parser.T_MAX

deltaT = (T_MAX-T_MIN)/(parser.ARRAY_MAX+1)

T_MIN = T_MIN + deltaT*parser.ARRAY_POS
T_MAX = T_MIN + deltaT*(parser.ARRAY_POS+1)

T_MAX_ITER = parser.T_MAX_ITER
T_MAX_ITER = math.ceil(T_MAX_ITER/(parser.ARRAY_MAX+1))

t_step = (T_MAX-T_MIN)/T_MAX_ITER

try: 
    os.mkdir("OUTPUT")
except Exception as ex:
    print("Output dir exists...")

with open(f"OUTPUT/output_{parser.ARRAY_POS}.dat", "w") as f:
    with redirect_stdout(f):
        if parser.ARRAY_POS == 0:
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
