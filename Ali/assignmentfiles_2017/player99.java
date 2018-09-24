import org.vu.contest.ContestSubmission;
import org.vu.contest.ContestEvaluation;

import java.util.Random;
import java.util.Properties;
import java.util.Arrays;

public class player99 implements ContestSubmission
{
	Random rnd_;
	ContestEvaluation evaluation_;
    private int evaluations_limit_;
    double range;
    
    double[][] population;
    double[][] parents;
    
	public player99()
	{
		rnd_ = new Random();
        // the range of each dimension is between - 5 to 5
        range = 5;
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
        this.population = new double[size][10];
        
        for (int i=0; i< this.population.length;i++)
        {
            for (int j=0; j< this.population[i].length;j++)
            {
     
                this.population[i][j] = range * rnd_.nextInt(2) * rnd_.nextDouble();
               
            }
        }
    }
    
    public void initialize_parents(int size)
    {
        parents = new double[size][10];
        
        for (int i=0; i< parents.length;i++)
        {
            for (int j=0; j< parents[i].length;j++)
            {

                parents[i][j] = 0;
            }
        }
    }
    
    
    
    public double get_fitness(double[] individual)
    {
        return (double) evaluation_.evaluate(individual);
    }
    
    
    public double[] crossover(double[] parent1,double[] parent2)
    {
        double child[] = new double[10];
        int crossPoint = 4+rnd_.nextInt(2);
        for (int i=0; i<10;i++)
        {
            if (i<crossPoint)
            { child[i] = parent1[i];}
            else
            {child[i] = parent2[i];}
        
        }
        return child;
    }
    
    public double[] mutate(double[] child)
    {   //crossover point
        int point1 = rnd_.nextInt(10);
        
        int point2 = rnd_.nextInt(10);
        while ( point2 == point1)
        {
            point2 = rnd_.nextInt(10);
        }
        double temp = child[point1];
        child[point1] = child[point2];
        child[point1] = temp;
        return child;
    }
    
    
    
	public void run()
	{
        int population_size = 20;
        int parents_size = 100;
        initialize_population(population_size);
        
        int evals = 0;
        int temp = 0;
        double individual;
        
        double fitness[][] = new double[population_size][2];
        
        double  parent[] = new double[parents_size];
        
        System.out.println( evaluations_limit_ );
        
        evaluations_limit_ = 2000;
        
//        while(evals<evaluations_limit_){
        
            //select the fittest parents for variation
            for (int i=0; i <population_size;i++)
            {
                fitness[i][0] = i;
                fitness[i][1] = get_fitness(population[i]);
                
                evals++;
            }
        //sort based on fitness
        for(int i=1; i < population_size;i++)
        {
            int index = i;
            double indexValue = fitness[i][1];
            
            for (int j=i; j>0; j--)
            {
                if  (fitness[j-1][1] > indexValue)
                {
                    fitness[j][1] = fitness[j-1][1];
                }else
                {
                    fitness[j][1] = indexValue;
                }
            }
        }
            
            
            for (int i=0;i <population_size; i++)
            {
                System.out.println(fitness[i][1]);
            }
//            System.out.println("\n\n\n");
            // sort population
            //pick the fittest as parent
            
            
           
    
      //  }
        //select fittests as the parents
        //crossover parents
        //mutate child
        //select a survival mechanism

	}
}
