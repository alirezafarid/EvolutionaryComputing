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

for i in range(S):
    plt.plot(scores[i], color="k", alpha=0.1/math.sqrt(S))
    
plt.plot(scores.mean(0), color="r")
plt.fill_between(np.linspace(0, N-1, num=N), scores.mean(0) + scores.std(0), scores.mean(0) - scores.std(0), color="r", alpha=0.1)
#plt.plot(scores.mean(0) - scores.std(0))
#plt.plot(scores.mean(0) + scores.std(0))

plt.xlim(0, N-1)
plt.ylim(-1, 11)

plt.xlabel("Generations")
plt.ylabel("Score")

plt.show()
