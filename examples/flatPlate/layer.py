import pandas as pd
import matplotlib.pyplot as plt
import numpy as np



def blasius_system(n, f):
    N = len(n)
    A = np.zeros((N, N))
    b = np.zeros(N)

    # Third derivarive part
    for i in range(2, N-2):
        h = (n[i+1] - n[i-1])/2
        # 2*f''' + f''f = 0
        A[i, i-2] = -0.5*2/h**3
        A[i, i-1] = 1*2/h**3
        A[i, i+1] = -1*2/h**3
        A[i, i+2] = 0.5*2/h**3
    # Second derivative part
    for i in range(1, N-1):
        h = (n[i+1] - n[i-1])/2
        # 2*f''' + f''f = 0
        A[i, i-1] += f[i]/h**2
        A[i, i] += f[i]*-2/h**2
        A[i, i+1] += f[i]/h**2
    # Boundary conditions
    # f(0) = 0
    A[0, 0] = 1
    b[0] = 0
    # f'(0) = 0
    h = n[1] - n[0]
    A[1, 0] = -1.5/h
    A[1, 1] = 2/h
    A[1, 2] = -0.5/h
    b[1] = 0
    # f'(inf) = 1
    h = n[N-1] - n[N-2]
    A[N-1, N-3] = 0.5/h
    A[N-1, N-2] = -2/h
    A[N-1, N-1] = 1.5/h
    b[N-1] = 1
    return A, b


def blasius_solution(n):
    # Solve 2*f''' + f''f = 0
    #   f'' = -2*f'''/f
    # f'(0) = 0, f(0) = 0, f'(inf) = 1
    f = n*1

    err = 1.
    while err > 1e-8:
        A, b = blasius_system(n, f)
        fp = np.linalg.solve(A, b)
        df = fp - f
        f += 0.5*df
        err = np.linalg.norm(df)

    return f



mu = 2e-5
pr = 0.72
cp = 1

x = 0.8

data = pd.read_csv("layer.csv")
rho = data["rho"]

nu = mu/rho

y_in = data["Points_1"]
u_in = data["U_0"]

U = np.max(u_in)

n = y_in * np.sqrt(U/(x*nu))
u = u_in / U


ns = np.linspace(0, 50, 1000)
fs = blasius_solution(ns)
us = np.gradient(fs, ns, edge_order=2)

plt.figure()
plt.plot(n[n<8], u[n<8])
plt.plot(ns[ns<8], us[ns<8])
plt.xlabel("$\eta$")
plt.ylabel("u/U")
plt.legend(["fvhyper", "blasius"])
plt.show()
