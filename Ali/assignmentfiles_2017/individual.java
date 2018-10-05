import java.util.Random;
import java.util.Properties;
import java.util.Arrays;

import java.util.*;

class individual
{
    int individual_id;
    int age;
    double fitness_value;
    double[] feature_vector;
    
    public individual()
    {
        this.individual_id = 0;
        this.fitness_value = 0.0;
        feature_vector = new double [10];
        this.age = 0;
        
    }
    
    public void print_individual()
    {
        for (int i =0; i < 10; i++)
            System.out.println(feature_vector[i]);
         System.out.println("\n\n");
    }
    
    public int get_id()
    {
        return individual_id;
    }
    public void set_id(int ID)
    {
        this.individual_id = ID;
    }
    public double get_fitness()
    {
        return fitness_value;
    }
    public void set_fitness(double fitness)
    {
        this.fitness_value = fitness ;
    }
    public double[] get_vector()
    {
        return feature_vector;
    }
    public void set_vector(double[] vector)
    {
        this.feature_vector = vector ;
    }
    
    public int get_age()
    {
        return this.age;
    }
    public void set_age(int age)
    {
        this.age = age ;
    }
   
   // @Override
   // public int compareTo(Object comparestu) {
   //     double compareFitness= ((individual)comparestu).get_fitness();
   //     /* For Ascending order*/   
   //     
   //    /* For Descending order do like this */
   //     
   // }
}

