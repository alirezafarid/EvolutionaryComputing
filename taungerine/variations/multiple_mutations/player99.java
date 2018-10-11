import org.vu.contest.ContestSubmission;
import org.vu.contest.ContestEvaluation;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.Math;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Properties;
import java.util.Random;

public class player99 implements ContestSubmission
{
	Random rnd_;
	ContestEvaluation evaluation_;
    private int evaluations_limit_;
    boolean[] function = new boolean[4];
	
	public player99()
	{
		rnd_ = new Random();
	}
	
	public void setSeed(long seed)
	{
		// Set seed of algortihms random process
		rnd_.setSeed(seed);
	}

	public void setEvaluation(ContestEvaluation evaluation)
	{
		// Set evaluation problem used in the run
		evaluation_ = evaluation;
		
		// Get evaluation properties
		Properties props = evaluation.getProperties();
        // Get evaluation limit
        evaluations_limit_ = Integer.parseInt(props.getProperty("Evaluations"));
		// Property keys depend on specific evaluation
		// E.g. double param = Double.parseDouble(props.getProperty("property_name"));
        boolean isMultimodal = Boolean.parseBoolean(props.getProperty("Multimodal"));
        boolean hasStructure = Boolean.parseBoolean(props.getProperty("Regular"));
        boolean isSeparable  = Boolean.parseBoolean(props.getProperty("Separable"));

		// Do sth with property values, e.g. specify relevant settings of your algorithm
        if(isMultimodal){
            // Do sth
        }else{
            // Do sth else
        }
        
        if(hasStructure){
            // Do sth
        }else{
            // Do sth else
        }
        
        if(isSeparable){
            // Do sth
            function[0] = true;
        }else{
            // Do sth else
        }
        
        if (!isMultimodal && !hasStructure) {
            function[1] = true;
        }
        
        if (isMultimodal && hasStructure) {
            function[2] = true;
        }
        
        if (isMultimodal && !hasStructure) {
            function[3] = true;
        }
    }
    
    private void printBest(Individual[][] pop) {
        double best = pop[0][0].fitness;
        for (int i = 1; i < pop.length; i++) {
            double f = pop[i][0].fitness;
            if (best < f) {
                best = f;
            }
        }
        System.out.println(best);
    }
    
	public void run()
	{
		// parameters
        int islands       = 4;    // number of population islands
        int pop_size      = 250; // number of individuals per island
        int exchange      = 5;   // number of individuals to exchange
        double sigma_init = 1.0;  // initial sigma
        double pressure   = 0.25; // selection pressure
        int children      = 100;  // children per pair
        
        if (function[2]) children = 100;
        if (function[3]) sigma_init = 0.01;
        
        // initialize population
        Individual[][] pop = new Individual[islands][pop_size];
        for (int i = 0; i < islands; i++) {
            for (int j = 0; j < pop_size; j++) {
                pop[i][j]  = new Individual(i*pop_size+j, sigma_init, rnd_);
            }
        }
        
        // calculate initial fitness
        int evals = 0;
        for (int i = 0; i < islands; i++) {
            for (int j = 0; j < pop_size; j++) {
                pop[i][j].fitness = (double) evaluation_.evaluate(pop[i][j].x);
                evals++;
            }
        }
        
        // sort by fitness
        for (int i = 0; i < islands; i++) {
            Arrays.sort(pop[i]);
        }
        
        while(evals<evaluations_limit_){
            for (int n = 0; n < islands; n++) {
                double r = rnd_.nextDouble();
                for (int i = 0; i < pop_size; i++) {
                    if (r < 1 - Math.exp(-(i+1)*pressure)) {
                        // individual i selected for mating
                        Individual father = pop[n][i];
                        
                        for (int j = 0; j < pop_size; j++) {
                            if (i != j) {
                                r = rnd_.nextDouble();
                                if (r < father.sdf0(pop[n][j])) {
                                    // individual j selected for mating with individual i
                                    Individual mother = pop[n][j];
                                    
                                    // print fitness of best current individual
                                    printBest(pop);

                                    // create children of individuals i and j
                                    for (int k = 0; k < children; k++) {
                                        Individual child = pop[n][pop_size-k-1];
                                        child.crossover3(father, mother);
                                        child.mutation();
                                        
                                        if (evals<evaluations_limit_) {
                                            // calculate child's fitness
                                            child.fitness = (double) evaluation_.evaluate(child.x);
                                            evals++;
                                        }
                                        double   best_fitness = child.fitness;
                                        double[] best_x       = child.x;
                                        for (int m = 0; m < 1; m++) {
                                            child.mutationB();
                                            if (evals<evaluations_limit_) {
                                                // calculate child's fitness
                                                child.fitness = (double) evaluation_.evaluate(child.x);
                                                evals++;
                                            }
                                            if (best_fitness < child.fitness) {
                                                best_fitness = child.fitness;
                                                best_x       = child.x;
                                            }
                                        }
                                        child.fitness = best_fitness;
                                        child.x       = best_x;
                                    }
                                    // sort island by fitness
                                    Arrays.sort(pop[n]);
                                    
                                    break;
                                }
                            }
                        }
                        break;
                    }
                }
            }
            if (1 < islands) {
                // exchange random individuals between islands
                for (int n = 0; n < islands; n++) {
                    for (int i = 0; i < exchange; i++) {
                        
                        // select random individual from this island
                        int r_i = rnd_.nextInt(pop_size);
                        Individual tmp = pop[n][r_i];
                        
                        int r_n;
                        do {
                            r_n = rnd_.nextInt(islands);
                        } while (r_n == n);
                        int r_j = rnd_.nextInt(pop_size);
                        
                        // exchange with random individual from other island
                        pop[n][r_i]   = pop[r_n][r_j];
                        pop[r_n][r_j] = tmp;
                    }
                }
            }
        }
    }
}
