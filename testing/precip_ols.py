import numpy as np

y = [0.5,0.5,0.7]
x = [1805.751, 1780.286, 1509.772]

A = np.vstack([x, np.ones(len(x))]).T


m, c = np.linalg.lstsq(A, y)[0]

print m, c
