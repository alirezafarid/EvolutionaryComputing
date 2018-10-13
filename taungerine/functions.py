import math
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

d_max = math.sqrt(1000)

def step(d, p_0, r_0, p_1, r_1, p_2):
    plt.title("step function")
    p = np.zeros(len(d))
    for i in range(len(d)):
        if r_1 < d[i]:
            p[i] = p_2
        elif r_0 < d[i]:
            p[i] = p_1
        else:
            p[i] = p_0

    return p

def normal(d, mu, sigma):
    plt.title("normal distribution")
    return np.exp(-(d - mu)**2 / (2 * sigma**2))

def lognormal(d, mu, sigma):
    plt.title("log-normal distribution")
    return np.exp(-(np.log(d) - mu)**2 / (2 * sigma**2))

d = np.linspace(0.0, 100.0, num=10001)

#plt.plot(d, step(d, 0.5, 1.0 / 3.0 * d_max, 0.0, 2.0 / 3.0 * d_max, 0.0))
#plt.plot(d, step(d, 0.0, 1.0 / 3.0 * d_max, 0.5, 2.0 / 3.0 * d_max, 0.0))
#plt.plot(d, step(d, 0.0, 1.0 / 3.0 * d_max, 0.0, 2.0 / 3.0 * d_max, 0.5))

#plt.plot(d, normal(d, 0.0, 0.125 * d_max))
#plt.plot(d, normal(d, 0.0, 0.25  * d_max))
#plt.plot(d, normal(d, 0.0, 0.5   * d_max))

#plt.plot(d, normal(d, d_max / 2.0, 0.125 * d_max))
#plt.plot(d, normal(d, d_max / 2.0, 0.25  * d_max))
#plt.plot(d, normal(d, d_max / 2.0, 0.5   * d_max))

#plt.plot(d, normal(d, d_max, 0.25 * d_max))
#plt.plot(d, normal(d, d_max, 0.5  * d_max))
#plt.plot(d, normal(d, d_max, 1.0  * d_max))

plt.plot(d, 1-normal(d, 0.0, 0.125 * d_max))
plt.plot(d, 1-normal(d, 0.0, 0.25  * d_max))
plt.plot(d, 1-normal(d, 0.0, 0.5   * d_max))

plt.plot(d, lognormal(d, 1.0, 1.0))
plt.plot(d, lognormal(d, 2.0, 1.0))
plt.plot(d, lognormal(d, 3.0, 1.0))

plt.xlim(0.0, d_max)
plt.ylim(-0.1, 1.1)
plt.xlabel("distance")
plt.ylabel("probability")

plt.show()
