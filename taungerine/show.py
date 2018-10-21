import math
import matplotlib.pyplot as plt
import numpy as np

root = 'results/'

def parse(i, N):
    for j in range(N):
        line            = lines[j].split(" ")
        scores[i][j]    = float(line[1])
        if best_individual[0] < scores[i][j]:
            best_individual[0] = scores[i][j]
            for k in range(10):
                best_individual[k+1] = line[k+2]
    
        distances[i][j] = float(line[0])

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
with open(root + 'storeresults0.txt') as f:
    lines = f.readlines()

# number of outputs, excluding final score
N = len(lines) - 2

# initialize scores
scores = np.zeros([S, N])

# initialize distances
distances = np.zeros([S, N])

# initialize best individual
best_individual = [-math.inf]*11;

# parse first sample
parse(0, N)

# read other samples
for i in range(1, S):
    with open(root + 'storeresults%d.txt' % i) as f:
        lines = f.readlines()

    # parse sample
    parse(i, N)

### PLOT SCORES ################################################################

avg_final_loss = 10.0 - final_scores.mean()
if avg_final_loss == 0.0:
    log_avg_loss = -math.inf
else:
    log_avg_loss = math.log10(avg_final_loss)
std_final_scores = final_scores.std()
if avg_final_loss == 0.0:
    log_std_loss = -math.inf
else:
    log_std_loss = math.log10(std_final_scores)
print("Log Average Loss:    %8g      Log SD: %8g   " % (log_avg_loss,    log_std_loss), end="")
print("  [%g, %g]" % (final_scores.min(), final_scores.max()))
print("    Average Runtime: %8g ms       SD: %8g ms" % (runtimes.mean(), runtimes.std()))
print("Best Individual")
print(best_individual[0])
for i in range(10):
    print(best_individual[i+1])

losses = 10.0 - scores

for i in range(S):
    plt.plot(losses[i], color="k", alpha=0.1/math.sqrt(S))

# plot mean taken in linear space
loss_mean = losses.mean(0)
plt.plot(losses.std(0), label=r"$\sigma$")
plt.plot(loss_mean, color="r", label=r"$\mu$")

# plot std taken in log space
#loss_std = np.exp(np.log(losses).std(0))
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

### PLOT DISTANCES #############################################################

for i in range(S):
    plt.plot(distances[i], color="k", alpha=0.5/S)

# plot mean taken in linear space
distances_mean = distances.mean(0)
plt.plot(distances.std(0), color="orange", label=r"$\sigma$")
plt.plot(distances_mean, color="g", label=r"$\mu$")

plt.xlim(0, N-1)
plt.ylim(10E-11, math.sqrt(1000.0))

plt.xlabel("Pairings")
plt.ylabel("Distance")

plt.yscale("log")

plt.legend()

plt.show()
