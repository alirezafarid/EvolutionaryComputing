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
        double best   = pop[0][0].fitness;
        int    best_n = 0;
        for (int n = 1; n < pop.length; n++) {
            double f = pop[n][0].fitness;
            if (best < f) {
                best   = f;
                best_n = n;
            }
        }
        System.out.print(" ");
        System.out.print(best);
        for (int i = 0; i < 10; i++) {
            System.out.print(" ");
            System.out.print(pop[best_n][0].x[i]);
        }
        System.out.println();
    }
    
	public void run()
	{
		// default parameters
        int islands       = 10;  // number of population islands
        int pop_size      = 5;   // number of individuals per island
        int exchange      = 2;   // number of individuals to exchange
        double sigma_init = 1.0; // initial sigma
        double pressure   = 0.5; // selection pressure
        int children      = 1;   // children per pair
        
        // bespoke parameters
        if (function[0]) {
            // SphereEvaluation
            // default
            
        } else if (function[1]) {
            // BentCigarFunction
            // default
            
        } else if (function[2]) {
            // SchaffersEvaluation
            islands    = 80;
            pop_size   = 2;
            exchange   = 1;
            sigma_init = 0.5;
            pressure   = 2.0;
        
        } else if (function[3]) {
            // KatsuuraEvaluation
            islands    = 1;
            pop_size   = 5000;
            exchange   = 50;
            sigma_init = 0.01;
            pressure   = 0.25;
            children   = 100;
        
        }
        
        // initialize population
        Individual[][] pop = new Individual[islands][pop_size];
        for (int i = 0; i < islands; i++) {
            for (int j = 0; j < pop_size; j++) {
                pop[i][j]  = new Individual(i*pop_size+j, sigma_init, rnd_);
            }
        }
        
        // initialize king of the hill (best known solution) (used only for the online contest)
        for (int n = 0; n < 0; n++) {
            Individual king = pop[n][0];
            if (function[0]) {
                // SphereEvaluation
                for (int i = 0; i < 10; i++) {
                    king.x[i] = 0.0;
                }
                
            } else if (function[1]) {
                // BentCigarFunction
                king.x[0] = -0.8905925387573322;
                king.x[1] =  3.990939725568024;
                king.x[2] =  0.1750075270168977;
                king.x[3] = -3.802956450666941;
                king.x[4] = -0.47512647426215243;
                king.x[5] = -2.091568842970519;
                king.x[6] =  1.385563099525254;
                king.x[7] = -0.7344883846697737;
                king.x[8] =  1.1424466100299642;
                king.x[9] = -0.30334927927877997;
                
            } else if (function[2]) {
                // SchaffersEvaluation
                king.x[0] =  3.65599992829552;
                king.x[1] =  2.5495999819725057;
                king.x[2] = -1.5296000114776622;
                king.x[3] =  1.4696000198400503;
                king.x[4] =  1.3959999270541783;
                king.x[5] = -1.9079999874954796;
                king.x[6] =  3.501600076536934;
                king.x[7] = -2.3504000192512433;
                king.x[8] = -0.38400000191915734;
                king.x[9] = -2.0359999988926867;
                
            } else if (function[3]) {
                // KatsuuraEvaluation
                king.x[0] =   1.1037196250834258;
                king.x[1] =   1.0554968388090176;
                king.x[2] =   0.8901799189934814;
                king.x[3] =  -1.9621457782643699;
                king.x[4] =  -0.1486274750313281;
                king.x[5] =  -1.2090665653828514;
                king.x[6] =   1.0939294722088457;
                king.x[7] =   0.4365273299171618;
                king.x[8] =  -2.211612145175111;
                king.x[9] =   0.02011991541794217;
                
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
            double r = rnd_.nextDouble();
            for (int i = 0; i < pop_size; i++) {
                if (r < 1 - Math.exp(-(i+1)*pressure)) {
                    // individual i selected as mother
                    Individual mother = pop[n][i];
                    
                    // search for father
                    for (int j = 0; j < pop_size; j++) {
                        if (i != j && rnd_.nextDouble() < mother.smf0(pop[n][j])) {
                            // individual j selected as father
                            Individual father = pop[n][j];
                            
                            // print distance between parents
                            System.out.print(mother.d(father));
                            
                            // print fitness and coordinates of best current individual
                            printBest(pop);
                            
                            // create children of mother and father
                            for (int k = 0; k < children; k++) {
                                Individual child = pop[n][pop_size-k-1];
                                child.crossover(mother, father);
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
