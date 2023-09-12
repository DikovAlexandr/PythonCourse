import matplotlib.pyplot as plt
import numpy as np

U = np.loadtxt("solve.csv")
U_0 = np.loadtxt("analytical.csv")
X = np.loadtxt("X.csv")
plt.figure(figsize=(10, 10))
plt.plot(X, U, label='Heat')
plt.plot(X, U_0, label='Analytic')
plt.legend()
plt.savefig('graph.png')
plt.show()
