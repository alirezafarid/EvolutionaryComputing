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
            //System.out.println("It's multimodal!");
        }else{
            // Do sth else
        }
        
        if(hasStructure){
            // Do sth
            //System.out.println("It's got structure!");
        }else{
            // Do sth else
        }
        
        if(isSeparable){
            // Do sth
            //System.out.println("It's separable!");
        }else{
            // Do sth else
        }
    }
    
	public void run()
	{
		// Run your algorithm here
        int evals    = 0;
        int pop_size = 100;
        
        // init population
        Individual[] parents  = new Individual[pop_size];
        Individual[] children = new Individual[pop_size];
        for (int i = 0; i < pop_size; i++) {
            parents[i]  = new Individual(i, rnd_);
            children[i] = new Individual(i, rnd_);
        }
        
        // calculate fitness
        for (int i = 0; i < pop_size; i++) {
            parents[i].fitness = (double) evaluation_.evaluate(parents[i].x);
            evals++;
        }
        
        // sort by fitness
        Arrays.sort(parents);
        
        while(evals<evaluations_limit_){
            // Select parents
            // Apply crossover / mutation operators
            //for (int i = 0; i < pop_size; i++) {
            //    children[i].crossover(parents[0], parents[1]);
            //    children[i].mutate();
            //}
            for (int i = 0; i < 10; i++) {
                for (int j = 0; j < 10; j++) {
                    children[i*10+j].crossover(parents[i], parents[j]);
                    children[i*10+j].mutate();
                }
            }
            
            // Check fitness of unknown fuction
            for (int i = 0; i < pop_size; i++) {
                children[i].fitness = (double) evaluation_.evaluate(children[i].x);
                evals++;
            }
            // Select survivors
            Arrays.sort(children);
            
            Individual[] tmp = parents;
            parents          = children;
            children         = tmp;
        }
    }
}
