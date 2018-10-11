import java.lang.Math;
import java.util.Random;

class Individual implements Comparable<Individual>
{
    static final int N = 10;
    static final int K = N * (N - 1) / 2;

    int      id;
    Random   rnd;
    
    double[] x       = new double[N];
    double   fitness = Double.NEGATIVE_INFINITY;
    
    // CMA-ES variables
    double tau_prime = 1.0 / Math.sqrt(2 * N);
    double tau       = 1.0 / Math.sqrt(2 * Math.sqrt(N));
    double beta      = 5.00;
    double epsilon   = 0.01;
    double[] sigma   = new double[N];
    double[] alpha   = new double[K];
    
    public Individual(int id, double sigma_init, Random rnd)
    {
        this.id  = id;
        this.rnd = rnd;
        
        // initialize sigma
        for (int i = 0; i < N; i++) {
            this.sigma[i] = sigma_init;
        }
        
        // initialize alpha
        for (int i = 0; i < K; i++) {
            this.alpha[i] = 0.0;
        }
        
        // initialize x
        for (int i = 0; i < N; i++) {
            this.x[i] = rnd.nextDouble() * 10.0 - 5.0;
        }
    }
    
    @Override
    public int compareTo(Individual other) {
        return (int) Math.signum(other.fitness - this.fitness);
    }
    
    public void crossover(Individual a, Individual b)
    {
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
            x[i]     = (r * a.x[i]     + (1.0 - r) * b.x[i]     ) / 2;
        }
    }
    
    public void crossover2(Individual a, Individual b)
    {
        int[] source = new int[N];
        for (int i = 0; i < N; i++) {
            source[i] = rnd.nextInt(2);
        }
        
        // recombine sigma and x
        for (int i = 0; i < N; i++) {
            sigma[i] = source[i]*a.sigma[i] + (1 - source[i])*b.sigma[i];
            x[i]     = source[i]*a.x[i]     + (1 - source[i])*b.x[i];
        }
        
        // recombine alpha
        for (int i = 1; i < N; i++) {
            for (int j = 0; j < i; j++) {
                int k = i * (i - 1) / 2 + j;
                alpha[k] = source[i]*source[j]*a.alpha[k] + (1 - source[i])*(1 - source[j])*b.alpha[k];
            }
        }
    }
    
    public void crossover3(Individual a, Individual b)
    {
        double[] source = new double[N];
        for (int i = 0; i < N; i++) {
            source[i] = rnd.nextDouble();
        }
        
        // recombine sigma and x
        for (int i = 0; i < N; i++) {
            sigma[i] = source[i]*a.sigma[i] + (1 - source[i])*b.sigma[i];
            x[i]     = source[i]*a.x[i]     + (1 - source[i])*b.x[i];
        }
        
        // recombine alpha
        for (int i = 1; i < N; i++) {
            for (int j = 0; j < i; j++) {
                int k = i * (i - 1) / 2 + j;
                alpha[k] = source[i]*source[j]*a.alpha[k] + (1 - source[i])*(1 - source[j])*b.alpha[k];
            }
        }
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
                //C[j][i] = C[i][j];
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
    
    public double[][] cholesky(double[][] C) {
        double[][] L = new double[N][N];
        
        for (int i = 0; i < N; i++)  {
            for (int j = 0; j <= i; j++) {
                double sum = 0.0;
                for (int k = 0; k < j; k++) {
                    sum += L[i][k] * L[j][k];
                }
                if (i == j) {
                    if (sum < C[i][i]) {
                        L[i][i] = Math.sqrt(C[i][i] - sum);
                    } else {
                        // matrix is not positive definite
                        for (int k = 0; k < K; k++) {
                            alpha[k] = 0.0;
                        }
                        return diag(sigma);
                    }
                } else {
                    L[i][j] = 1.0 / L[j][j] * (C[i][j] - sum);
                }
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
        System.out.println(diff);
    }
    
    public void mutation()
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
                x[i] += L[i][j] * v[j];
                if (x[i] < -5.0) {
                    x[i] = -5.0;
                } else if (5.0 < x[i]) {
                    x[i] = 5.0;
                }
            }
        }
    }
    
    private double norm(Individual a, Individual b, int L) {
        double sum = 0.0;
        for (int i = 0; i < 10; i++) {
            double y = 1.0;
            for (int j = 0; j < L; j++) {
                y *= Math.abs(b.x[i] - a.x[i]);
            }
            sum += y;
        }
        
        return Math.pow(sum, 1.0 / L);
    }
    
    public double step(Individual b, int L, double r, double p_0, double p_1) {
        double p;
        
        double norm = norm(this, b, L);
        
        if (norm < r) {
            p = p_0;
        } else {
            p = p_1;
        }
        
        return p;
    }
    
    public double normal(Individual b, int L, double sigma) {
        double norm = norm(this, b, L);
        
        return Math.exp(-norm*norm / (2 * sigma*sigma));
    }
    
    public double lognormal(Individual b, int L, double mu, double sigma) {
        double norm = norm(this, b, L);
        double arg  = Math.log(norm) - mu;
        
        return Math.exp(-arg*arg / (2 * sigma*sigma));
    }
    
    // example standard distance function
    public double sdf0(Individual b) {
        int L        = 2;
        double sigma = 10.0;
        
        return this.normal(b, L, sigma);
    }
}
