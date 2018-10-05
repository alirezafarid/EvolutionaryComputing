import org.vu.contest.ContestSubmission;
import org.vu.contest.ContestEvaluation;

import java.util.Random;
import java.util.Properties;
import java.util.Arrays;
import java.util.Collections;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import java.util.*;

public class player99 implements ContestSubmission
{
    Random rnd_;
    ContestEvaluation evaluation_;
    private int evaluations_limit_;

    double range;
    int evals;
    
    ArrayList<individual> population;
    individual indv; 
    
	public player99()
	{
         
	 rnd_ = new Random();
         population = new ArrayList<individual>();

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
        boolean isSeparable = Boolean.parseBoolean(props.getProperty("Separable"));

		// Do sth with property values, e.g. specify relevant settings of your algorithm
        if(isMultimodal){
            // Do sth
        }else{
            // Do sth else
        }
    }
    
    public void initialize_population(int size)
    {
        int counter =0;
        double  upper = 5.0;
        double lower = -5.0;
        for(int i = 0; i < size; i++){
            
            indv = new individual();
            double vector[] = new  double[10];
            for (int j=0; j< vector.length;j++)
            {
                vector[j]= this.rnd_.nextDouble() * (upper - lower) + lower;
            } 
            
            indv.set_vector(vector);
            indv.set_age(1);
            indv.set_id(counter++);
            population.add(indv);

        }
        
    }
    
    public void caclcuate_fitness()
    {
        for(individual indv: this.population)
        {
            indv.set_fitness((double)evaluation_.evaluate(indv.get_vector()));
            this.evals++;
        }
        
        Collections.sort(this.population,new CustomComparator(-1));
     }
    
    public double[] crossover(double[] vector1,double[] vector2)
    {
        double child_vector[] = new double[10];
        int crossPoint = 4+rnd_.nextInt(2);
        for (int i=0; i<10;i++)
        {
            if (i<crossPoint)
            { child_vector[i] = vector1[i];}
            else
            {child_vector[i] = vector2[i];}
        
        }
        return child_vector;
    }
    
    public double[] bitwise_crossover(double[] vector1,double[] vector2)
    {
        double child_vector[] = new double[10];

        for (int i=0; i<10;i++)
        {   if (i % 2 == 0 )
            { child_vector[i] = vector1[i];}
            else
            {child_vector[i] = vector2[i];}
        }
        return child_vector;
    }
    

    public double[] mutate_v(double[] child_vector)
    {
        for (int i=0; i <10; i++)
        {
            if (rnd_.nextDouble() > 0.9)
            {
                double v = child_vector[i] + (2 * rnd_.nextDouble() - 1);
                if (v >= -5 && v <= 5)
                {
                    child_vector[i] = child_vector[i] + v;
                }
            }
            
        }
        
        return child_vector;
    }
    
	public void run()
	{
        int population_size = 1500;
        int num_children = 900;
        evals = 0;
        
        
        initialize_population(population_size);

        int eval_limit = evaluations_limit_ - population_size - num_children;
         while(evals< eval_limit)
         {
        
             //Fitness Calculation
             this.caclcuate_fitness();
            
            //Parent Selection And Vartion
            
            for(int i=0; i<num_children*2; i=i+2)
            {
                double[] parent1_vector = population.get(i).get_vector();
                double[] parent2_vector = population.get(i+1).get_vector();
            
                double[] child_vector =mutate_v(
                 crossover(parent1_vector, parent2_vector));
            
                individual indv = new individual();
                indv.set_vector(child_vector);
                indv.set_fitness((double)evaluation_.evaluate(indv.get_vector()));
                this.evals++;
                indv.set_age(1);
                indv.set_id(population_size+i);
                this.population.add(indv);
            

            }

            //Survial Selection, birth rate = death rate
             
            Collections.sort(this.population,new CustomComparator(1));
            for (int i=0; i<num_children;i++)
            {
                population.remove(i);
            }
           
        
             
        }
        

    
	}
}
