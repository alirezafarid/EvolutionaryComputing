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
        int islands       = 1;    // number of population islands
        int pop_size      = 1000; // number of individuals per island
        int exchange      = 5;    // number of individuals to exchange
        double sigma_init = 1.0;  // initial sigma
        double pressure   = 0.25; // selection pressure
        int children      = 100;  // children per pair
        
        // bespoke parameters
        if (function[1]) {
            
        
        } else if (function[2]) {
            islands       = 5;
            pop_size      = 200;
            exchange      = 5;
            sigma_init    = 1.0;
            pressure      = 0.2;
            children      = 80;
        
        } else if (function[3]) {
            sigma_init = 0.01;
        
        }
        
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
        
        // first island
        int n = 0;
        
        while(evals < evaluations_limit_){
            
            // sort by fitness
            for (int m = 0; m < islands; m++) {
                Arrays.sort(pop[m]);
            }
            
            // search for mother
            for (int i = 0; i < pop_size; i++) {
                double r_i = rnd_.nextDouble();
                if (r_i < 1 - Math.exp(-(i+1)*pressure)) {
                    // individual i selected as mother
                    Individual mother = pop[n][i];
                    
                    // search for father
                    for (int j = 0; j < pop_size; j++) {
                        double r_j = rnd_.nextDouble();
                        if (i != j && r_j < mother.smf0(pop[n][j])) {
                            // individual j selected as father
                            Individual father = pop[n][j];
                            
                            // print fitness of best current individual
                            printBest(pop);
                            
                            // create children of mother and father
                            for (int k = 0; k < children; k++) {
                                Individual child = pop[n][pop_size-k-1];
                                child.crossover3(mother, father);
                                child.mutation();
                                
                                if (evals < evaluations_limit_) {
                                    // calculate child's fitness
                                    child.fitness = (double) evaluation_.evaluate(child.x);
                                    evals++;
                                }
                            }
                            break; // j loop
                            
                        } else {
                            // continue to next individual j
                            
                        }
                    }
                    break; // i loop
                    
                } else {
                    // continue to next individual i
                    
                }
            }
            if (1 < islands) {
                // exchange random individuals with other islands
                for (int i = 0; i < exchange; i++) {
                    
                    // select random individual from this island
                    int r_i = rnd_.nextInt(pop_size);
                    Individual tmp = pop[n][r_i];
                    
                    int r_m;
                    do {
                        r_m = rnd_.nextInt(islands);
                    } while (r_m == n);
                    int r_j = rnd_.nextInt(pop_size);
                    
                    // exchange with random individual from other island
                    pop[n][r_i]   = pop[r_m][r_j];
                    pop[r_m][r_j] = tmp;
                }
            }
            
            // next island
            n = ++n % islands;
        }
    }
}
