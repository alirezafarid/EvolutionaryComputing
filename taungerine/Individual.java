import java.lang.Math;
import java.util.Random;

class Individual implements Comparable<Individual>
{
    static final int N = 10;
    static final int K = N * (N - 1) / 2;
    
    double[] x = new double[N];
    double fitness;
    int age;
    int id;
    int rank;
    Random rnd;
    int is_parent;
    
    // CMA-ES variables
    double tau_prime = 0.22;
    double tau       = 0.39;
    double beta      = 5.00;
    double epsilon   = 0.001;
    double[] sigma   = new double[N];
    double[] alpha   = new double[K];
    
    public Individual(int id, double sigma, Random rnd)
    {
        this.id  = id;
        this.rnd = rnd;
        
        // initialize sigma
        for (int i = 0; i < N; i++) {
            this.sigma[i] = sigma;
        }
        
        // initialize alpha
        for (int i = 0; i < K; i++) {
            this.alpha[i] = 0.0;
        }
        
        // initialize x
        for (int i = 0; i < N; i++) {
            this.x[i] = rnd.nextDouble() * 10.0 - 5.0;
            //this.x[i] = rnd.nextDouble()*0.01;
        }
    }
    
    @Override
    public int compareTo(Individual other) {
        return (int) Math.signum(other.fitness - this.fitness);
    }
    
    
    // Add all setter and getter functuons
    
    
    public void crossover(Individual a, Individual b)
    {
        // set this.x as a function of a.x and b.x
        // crossover sigma
        for (int i = 0; i < N; i++) {
            double r = rnd.nextDouble();
            sigma[i] = (r * a.sigma[i] + (1.0 - r) * b.sigma[i] ) / 2;
        }
        
        // crossover alpha
        for (int i = 0; i < K; i++) {
            double r = rnd.nextDouble();
            alpha[i] = (r * a.alpha[i] + (1.0 - r) * b.alpha[i] ) / 2;
        }
        
        // crossover x
        for (int i = 0; i < N; i++) {
            double r = rnd.nextDouble();
            x[i]     = (r * a.x[i] + (1.0 - r) * b.x[i] ) / 2;
        }
    }
    
    public void crossover2(Individual a, Individual b)
    {
        boolean[] source = new boolean[N];
        for (int i = 0; i < N; i++) {
            source[i] = rnd.nextBoolean();
        }
        
        for (int i = 0; i < N; i++) {
            if (source[i]) {
                sigma[i] = a.sigma[i];
                x[i]     = a.x[i];
            } else {
                sigma[i] = b.sigma[i];
                x[i]     = b.x[i];
            }
        }
        
        for (int i = 1; i < N; i++) {
            for (int j = 0; j < i; j++) {
                int k = i * (i - 1) / 2 + j;
                if (source[i] && source[j]) {
                    alpha[k] = a.alpha[k];
                } else if (!source[i] && !source[j]) {
                    alpha[k] = b.alpha[k];
                } else {
                    alpha[k] = 0.0;
                }
            }
        }
    }
    
    public void crossover3(Individual a, Individual b)
    {
        double[] source = new double[N];
        for (int i = 0; i < N; i++) {
            source[i] = rnd.nextDouble();
        }
        
        for (int i = 0; i < N; i++) {
            sigma[i] = source[i]*a.sigma[i] + (1 - source[i])*b.sigma[i];
            x[i]     = source[i]*a.x[i]     + (1 - source[i])*b.x[i];
        }
        
        for (int i = 1; i < N; i++) {
            for (int j = 0; j < i; j++) {
                int k = i * (i - 1) / 2 + j;
                alpha[k] = source[i]*source[j]*a.alpha[k] + (1 - source[i])*(1 - source[j])*b.alpha[k];
            }
        }
    }
    
    public void mutationA()
    {
        // mutate sigma
        double r = rnd.nextGaussian();
        for (int i = 0; i < N; i++) {
            sigma[i] *= Math.exp(tau_prime * r + tau * rnd.nextGaussian());
            if (sigma[i] < epsilon) {
                sigma[i] = epsilon;
            }
        }
        
        // mutate alpha
        for (int i = 0; i < K; i++) {
            alpha[i] += beta * rnd.nextGaussian();
            if (Math.PI < Math.abs(alpha[i])) {
                alpha[i] -= 2 * Math.PI * Math.signum(alpha[i]);
            }
        }
    }
    
    public void mutationB()
    {
        // mutate x
        double[] v = new double[N];
        for (int i = 0; i < N; i++) {
            v[i] = rnd.nextGaussian();
        }
        double[][] C = covarianceMatrix();
        double[][] L = cholesky(C);
        //check(C, L);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (! Double.isNaN(L[i][j])) {
                    x[i] += L[i][j] * v[j];
                    if (x[i] < -5.0) {
                        x[i] = -5.0;
                    } else if (5.0 < x[i]) {
                        x[i] = 5.0;
                    }
                }
            }
        }
    }
    
    public void mutation()
    {
        mutationA();
        mutationB();
    }
    
    private double[][] covarianceMatrix() {
        double[][] C = new double[N][N];
        
        // fill diagonals
        for (int i = 0; i < N; i++) {
            C[i][i] = sigma[i] * sigma[i];
        }
        
        // fill off-diagonals
        for (int i = 1; i < N; i++) {
            for (int j = 0; j < i; j++) {
                int k = i * (i - 1) / 2 + j;
                C[i][j] = (sigma[i] * sigma[i] - sigma[j] * sigma[j]) / 2 * Math.tan(2*alpha[k]);
                C[j][i] = C[i][j];
            }
        }
        
        return C;
    }
    
    private double[][] diag(double[] v) {
        double[][] A = new double[N][N];
        
        for (int i = 0; i < N; i++) {
            A[i][i] = v[i];
        }
        
        return A;
    }
    
    // CMA-ES functions //////////////////////////////////////////////
    // return Cholesky factor L of psd matrix A = L L^T
    public double[][] cholesky(double[][] C) {
        double[][] L = new double[N][N];
        
        for (int i = 0; i < N; i++)  {
            for (int j = 0; j <= i; j++) {
                double sum = 0.0;
                for (int k = 0; k < j; k++) {
                    sum += L[i][k] * L[j][k];
                }
                //if (i == j) L[i][i] = Math.sqrt(C[i][i] - sum);
                //else        L[i][j] = 1.0 / L[j][j] * (C[i][j] - sum);
                if (i == j) {
                    if (sum < C[i][i]) {
                        L[i][i] = Math.sqrt(C[i][i] - sum);
                    } else {
                        for (int k = 0; k < K; k++) {
                            alpha[k] = 0.0;
                        }
                        return diag(sigma);
                    }
                } else {
                    L[i][j] = 1.0 / L[j][j] * (C[i][j] - sum);
                }
            }
            if (L[i][i] <= 0) {
                //if (Double.isNaN(L[i][i])) {
                throw new RuntimeException("Matrix not positive definite");
            }
        }
        return L;
    }
    
    private void check(double[][] C, double[][]L) {
        double[][] A = new double[N][N];
        
        // multiply L L^T
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                for (int k = 0; k < N; k++) {
                    A[i][j] += L[i][k] * L[j][k];
                }
            }
        }
        
        double diff = 0.0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                diff += Math.abs(C[i][j] - A[i][j]);
            }
        }
        System.out.print("Difference: ");
        System.out.print(diff);
        System.out.println();
    }
    
}
