import java.util.*; 
import java.lang.*; 
import java.io.*; 

public class CustomComparator implements Comparator<individual> {
    
    int order;
    public CustomComparator (int order)
    {
     super();
     this.order = order;
     
    } 
    @Override
    public int compare(individual o1, individual o2) {
      
	    if (order == 1)		
             return Double.valueOf(o1.get_fitness()).compareTo(Double.valueOf(o2.get_fitness()));
            else
	     return Double.valueOf(o2.get_fitness()).compareTo(Double.valueOf(o1.get_fitness()));		
    
    }
}
