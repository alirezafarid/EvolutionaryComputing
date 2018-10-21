import math
import matplotlib.pyplot as plt
import numpy as np

def parse(i, N):
    for j in range(N):
        scores[i][j] = float(lines[j])

    line            = lines[-2].split(" ")
    final_scores[i] = float(line[1])

    line        = lines[-1].split(" ")
    runtimes[i] = int(line[1][:-3])

# number of samples
S = 100

# initialize final scores and runtimes
final_scores = np.zeros(S)
runtimes     = np.zeros(S, dtype="uint64")

# read first sample
with open('results/storeresults0.txt') as f:
    lines = f.readlines()

# number of outputs, excluding final score
N = len(lines) - 2

# initialize scores
scores = np.zeros([S, N])

# parse first sample
parse(0, N)

# read other samples
for i in range(1, S):
    with open('results/storeresults%d.txt' % i) as f:
        lines = f.readlines()

    # parse sample
    parse(i, N)

print("Average Score:   %7g      SD: %7g   " % (final_scores.mean(), final_scores.std()), end="")
print("  [%g, %g]" % (final_scores.min(), final_scores.max()))
print("Average Runtime: %7g ms   SD: %7g ms" % (runtimes.mean(),     runtimes.std()))

losses = 10.0 - scores

for i in range(S):
    plt.plot(losses[i], color="k", alpha=0.1/math.sqrt(S))

# plot mean taken in linear space
loss_mean = losses.mean(0)
plt.plot(losses.std(0), label=r"$\sigma$")
plt.plot(loss_mean, color="r", label=r"$\mu$")

# plot std taken in log space
loss_std = np.exp(np.log(losses).std(0))
#upper    = loss_mean * loss_std
#lower    = loss_mean / loss_std
#plt.fill_between(np.linspace(0, N-1, num=N), upper, lower, color="r", alpha=0.1)
#plt.fill_between(np.linspace(0, N-1, num=N), loss_mean + losses.std(0), loss_mean - losses.std(0), color="r", alpha=0.1)

plt.xlim(0, N-1)
#plt.ylim(10**np.floor(np.log10(lower.min())-1), 10E1)
plt.ylim(10E-17, 10E1)

plt.yscale("log")

plt.xlabel("Pairings")
plt.ylabel("Loss")

plt.legend()

plt.show()
