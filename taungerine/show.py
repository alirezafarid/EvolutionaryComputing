import math

# number of samples
S = 10

scores   = [0.0] * S
runtimes = [0]   * S

for i in range(S):
    with open('storeresults%d.txt' % i) as f:
        lines = f.readlines()

    line      = lines[0].split(" ")
    scores[i] = float(line[1])

    line        = lines[1].split(" ")
    runtimes[i] = int(line[1][:-3])

score_mean = sum(scores)/S
score_SD   = math.sqrt(sum([(score - score_mean)**2 for score in scores])/S)
print("Average Score:   %7g      SD: %7g" % (score_mean, score_SD))

runtime_mean = sum(runtimes)/S
runtime_SD   = math.sqrt(sum([(runtime - runtime_mean)**2 for runtime in runtimes])/S)
print("Average Runtime: %7g ms   SD: %7g ms" % (runtime_mean, runtime_SD))
