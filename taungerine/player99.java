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
    
    public double distance(Individual a, Individual b) {
        double distance = 0.0;
        
        for (int i = 0; i < 10; i++) {
            double d = a.x[i] - b.x[i];
            distance += d*d;
        }
        
        return Math.sqrt(distance);
    }
    
    public double mate(double distance) {
        double d = 10.0;
        return Math.exp(-(distance/d)*(distance/d));
    }
    
	public void run()
	{
		// Run your algorithm here
        double P = 1.0;
        int Q;
        if (function[2]) {
            Q = 400;
        } else {
            Q = 100;
        }
        
        int evals    = 0;
        int islands  = 1;
        int pop_size = 1000;
        
        // init population
        Individual[][] pop = new Individual[islands][pop_size];
        for (int i = 0; i < islands; i++) {
            for (int j = 0; j < pop_size; j++) {
                pop[i][j]  = new Individual(i*pop_size+j, rnd_);
            }
        }
        
        // calculate fitness
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
                    if (r < 1 - Math.exp(-(i+1)/P)) {
                        for (int j = 0; j < pop_size; j++) {
                            if (i != j) {
                                r = rnd_.nextDouble();
                                if (r < mate(distance(pop[n][i], pop[n][j]))) {
                                    for (int k = 0; k < Q; k++) {
                                        pop[n][pop_size-k-1].crossover3(pop[n][i], pop[n][j]);
                                        pop[n][pop_size-k-1].mutation();
                                        
                                        if (evals<evaluations_limit_) {
                                            pop[n][pop_size-k-1].fitness = (double) evaluation_.evaluate(pop[n][pop_size-k-1].x);
                                            evals++;
                                        }
                                    }
                                    Arrays.sort(pop[n]);
                                    
                                    break;
                                }
                            }
                        }
                        break;
                    }
                }
            }
        }
    }
}
