import java.lang.Math;
import java.util.Random;

class Individual implements Comparable<Individual>
{
    static final int N = 10;
    static final int K = N * (N - 1) / 2;

    int    id;
    Random rnd;
    
    double[] x       = new double[N];
    double   fitness = Double.NEGATIVE_INFINITY;
    
    static final double tau_prime = 1.0 / Math.sqrt(2 * N);
    static final double tau       = 1.0 / Math.sqrt(2 * Math.sqrt(N));
    static final double beta      = 5.0;
    static final double epsilon   = 0.0;
    
    double[] sigma = new double[N];
    double[] alpha = new double[K];
    
    static final double d_max = Math.sqrt(1000.0);
    
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
                    L[i][j] = (C[i][j] - sum) / L[j][j];
                }
            }
        }
        return L;
    }
    
    private void check(double[][] C, double[][] L) {
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
    
    private double d(Individual a, Individual b) {
        double sum = 0.0;
        for (int i = 0; i < 10; i++) {
            double d_i  = b.x[i] - a.x[i];
            sum += d_i*d_i;
        }
        
        return Math.sqrt(sum);
    }
    
    public double step(Individual b, double p_0, double r_0, double p_1, double r_1, double p_2) {
        double d = d(this, b);
        
        if (r_1 < d) {
            return p_2;
        } else if (r_0 < d) {
            return p_1;
        } else {
            return p_0;
        }
    }
    
    public double normal(Individual b, double mu, double sigma) {
        double d   = d(this, b);
        double arg = d - mu;
        
        return 0.9 * (Math.exp(-arg*arg / (2 * sigma*sigma)) + 1.0/9);
    }
    
    public double lognormal(Individual b, double mu, double sigma) {
        double d   = d(this, b);
        double arg = Math.log(d) - mu;
        
        return 0.9 * (Math.exp(-arg*arg / (2 * sigma*sigma)) + 1.0/9);
    }
    
    public double smf0(Individual b) {
        return 0.5;
    }
    
    public double smf1(Individual b) {
        double p_0 = 0.5;
        double r_0 = 1.0 / 3.0 * d_max;
        double p_1 = 0.1;
        double r_1 = 2.0 / 3.0 * d_max;
        double p_2 = 0.1;
        
        return this.step(b, p_0, r_0, p_1, r_1, p_2);
    }
    
    public double smf2(Individual b) {
        double p_0 = 0.1;
        double r_0 = 1.0 / 3.0 * d_max;
        double p_1 = 0.5;
        double r_1 = 2.0 / 3.0 * d_max;
        double p_2 = 0.1;
        
        return this.step(b, p_0, r_0, p_1, r_1, p_2);
    }
    
    public double smf3(Individual b) {
        double p_0 = 0.1;
        double r_0 = 1.0 / 3.0 * d_max;
        double p_1 = 0.1;
        double r_1 = 2.0 / 3.0 * d_max;
        double p_2 = 0.5;
        
        return this.step(b, p_0, r_0, p_1, r_1, p_2);
    }
    
    public double smf4(Individual b) {
        double mu    = 0.0;
        double sigma = d_max / 8.0;
        
        return this.normal(b, mu, sigma);
    }
    
    public double smf5(Individual b) {
        double mu    = 0.0;
        double sigma = d_max / 4.0;
        
        return this.normal(b, mu, sigma);
    }
    
    public double smf6(Individual b) {
        double mu    = 0.0;
        double sigma = d_max / 2.0;
        
        return this.normal(b, mu, sigma);
    }
    
    public double smf7(Individual b) {
        double mu    = d_max / 2.0;
        double sigma = d_max / 8.0;
        
        return this.normal(b, mu, sigma);
    }
    
    public double smf8(Individual b) {
        double mu    = d_max / 2.0;
        double sigma = d_max / 4.0;
        
        return this.normal(b, mu, sigma);
    }
    
    public double smf9(Individual b) {
        double mu    = d_max / 2.0;
        double sigma = d_max / 2.0;
        
        return this.normal(b, mu, sigma);
    }
    
    public double smf10(Individual b) {
        double mu    = d_max;
        double sigma = d_max / 4.0;
        
        return this.normal(b, mu, sigma);
    }
    
    public double smf11(Individual b) {
        double mu    = d_max;
        double sigma = d_max / 2.0;
        
        return this.normal(b, mu, sigma);
    }
    
    public double smf12(Individual b) {
        double mu    = d_max;
        double sigma = d_max;
        
        return this.normal(b, mu, sigma);
    }
    
    public double smf13(Individual b) {
        double mu    = 1.0;
        double sigma = 1.0;
        
        return this.lognormal(b, mu, sigma);
    }
    
    public double smf14(Individual b) {
        double mu    = 2.0;
        double sigma = 1.0;
        
        return this.lognormal(b, mu, sigma);
    }
    
    public double smf15(Individual b) {
        double mu    = 3.0;
        double sigma = 1.0;
        
        return this.lognormal(b, mu, sigma);
    }
}
