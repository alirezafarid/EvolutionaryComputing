import java.util.Properties;
import org.vu.contest.ContestEvaluation;
















public class SphereEvaluation
  implements ContestEvaluation
{
  private static final int EVALS_LIMIT_ = 10000;
  private static final double BASE_ = 11.5356D;
  private static final double ftarget_ = 0.0D;
  private double best_;
  private int evaluations_;
  private String multimodal_ = "false";
  private String regular_ = "true";
  private String separable_ = "true";
  private String evals_ = Integer.toString(10000);
  
  public SphereEvaluation()
  {
    best_ = 0.0D;
    evaluations_ = 0;
  }
  

  private double function(double[] paramArrayOfDouble)
  {
    double d = 0.0D;
    for (int i = 0; i < 10; i++) d += paramArrayOfDouble[i] * paramArrayOfDouble[i];
    return d;
  }
  


  public Object evaluate(Object paramObject)
  {
    if (!(paramObject instanceof double[])) throw new IllegalArgumentException();
    double[] arrayOfDouble = (double[])paramObject;
    if (arrayOfDouble.length != 10) { throw new IllegalArgumentException();
    }
    if (evaluations_ > 10000) { return null;
    }
    

    double d = 10.0D - 10.0D * ((function(arrayOfDouble) - 0.0D) / 11.5356D);
    if (d > best_) best_ = d;
    evaluations_ += 1;
    
    return new Double(d);
  }
  

  public Object getData(Object paramObject)
  {
    return null;
  }
  

  public double getFinalResult()
  {
    return best_;
  }
  

  public Properties getProperties()
  {
    Properties localProperties = new Properties();
    localProperties.put("Multimodal", multimodal_);
    localProperties.put("Regular", regular_);
    localProperties.put("Separable", separable_);
    localProperties.put("Evaluations", evals_);
    return localProperties;
  }
}
