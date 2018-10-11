import math
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

def step(d, r, p_0, p_1):
    p = np.zeros(len(d))
    for i in range(len(d)):
        if d[i] < r:
            p[i] = p_0
        else:
            p[i] = p_1

    plt.title("step function")
    return p

def normal(d, sigma):
    plt.title("normal distribution")
    return np.exp(-d**2 / (2 * sigma**2))

def lognormal(d, mu, sigma):
    plt.title("log-normal distribution")
    return np.exp(-(np.log(d) - mu)**2 / (2 * sigma**2))

d = np.linspace(0.0, 32.0, num=3201)

# step function
#p = step(d, 20.0, 0.7, 0.0)
p = normal(d, 10.0)
#p = lognormal(d, 1.5, 0.75)

plt.plot(d, p)

plt.xlim(0.0, math.sqrt(1000.0))
plt.ylim(-0.1, 1.1)
plt.xlabel("distance")
plt.ylabel("probability")

plt.show()
