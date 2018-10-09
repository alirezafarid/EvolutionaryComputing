import sys,os
import csv
import functools
import numpy as np
import matplotlib.pyplot as plt

path = "results/"

#---------------------------functions
def process_data(filename):
    if not os.path.exists(path + filename):
        print("Python Error: file path '"+  path + filename  +"' was not found!")
        return
    
    file_content = open(path + filename,"r").read()
    num_generations = file_content.split("ms")
    len(num_generations)
    
    data = []
    time = []; score =[]
    g = file_content.split("ms")
    for i in range(0,len(g)-1):
        row = []
        generation = g[i].split("\n")
        for j in generation:
            if ( j !=''):
                row.append(j)
        time.append(float(row[-1].split(":")[1].replace(" ","")))
        score.append(float(row[-2].split(":")[1].replace(" ","")))
        data.append(list(map(lambda x: float(x),row[:-2])))
    
    data_avg =[]
    for g in range(len(data[0])):
        temp = 0
        for i in range(len(data)):
            temp = temp + data[i][g]
        
        data_avg.append(temp / len(data))
    return score,time,data_avg

def plot_fitness(plt,data,data_label):
    plt.plot(data,label=data_label)
    plt.xlabel("Number of Generations")
    plt.ylabel("Fitness Value")
    plt.ylim(0,10)
    plt.legend()

def plot_time(plt,data,data_label):
    plt.plot(data,label=data_label)
    plt.xlabel("Experiment Number")
    plt.ylabel("Time")
    plt.legend()

def plot_score(plt,data,data_label):
    plt.plot(data,label=data_label)
    plt.xlabel("Experiment Number")
    plt.ylabel("Best Score")
    plt.legend()




#--------------------------- take the file name as the first argument

filename = sys.argv[1]
func_name = sys.argv[2]
if "Evaluation" in func_name:   
    func_name = func_name.split("Evaluation")[0] 
if "Function" in func_name:
    func_name = func_name.split("Function")[0] 
#--------------------------- returns the average best fitness value in each generation,
#list of execution time, and best scores.



result = process_data(filename)

score = result[0]
time = result[1]
data_avg = result[2]

print("average best score: " + str(np.average(score)) +"\n" )
print("STD of best score: " + str(np.asarray(score).std()) + "\n")
print("average execution time: " + str(np.average(time)) + "\n" )
print("STD of execution time: " + str(np.asarray(time).std()) + "\n")

#--------------------------- create a dictionary of result and write it in csv format

result_dict = {"score":score,"time":time,"score_per_generation":data_avg}
with open(path + filename +'_'+ func_name+'.csv', 'w') as f:
    for key in result_dict.keys():
        f.write("next line is %s:\n%s\n"%(key,','.join(str(e) for e in result_dict[key])))


#---------------------------create chart for best score
#---------------------------create chart for execution time
#---------------------------create chart for average fitness value after each generation

plt.figure()
plot_fitness(plt,data_avg,"Average Score Per Generation, " + func_name)
# for i in range(len(data)):
#     plot_fitness(plt,data[i],"Experiment "+ str(i))
plt.grid();plt.savefig(path + filename + '_' + func_name  +"_avg.png")


plt.figure()
plot_time(plt,time,"Execution Time, " + func_name)
plt.grid();plt.savefig(path + filename + '_' + func_name + "_time.png")

plt.figure()
plot_score(plt,score,"Score,"+ func_name)
plt.grid();plt.savefig(path + filename+  '_' + func_name + "_score.png")