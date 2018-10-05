import java.lang.Math;
import java.util.Random;

public class Individual implements Comparable<Individual> {
    static final int N = 10;
    static final int K = N * (N - 1) / 2;
    
    public double fitness  = Double.NEGATIVE_INFINITY;
    public double[] x      = new double[N];
    private double[] sigma = new double[N];
    private double[] alpha = new double[K];
    private Random rnd_;
    private int id;
    
    // hyperparameters
    double tau_prime = 0.22;//0.22;
    double tau       = 0.39;//0.39;
    double beta      = 5.00;//5.00;
    double epsilon   = 0.01;//0.01;
    
    public Individual(int id, Random rnd_) {
        this.id   = id;
        this.rnd_ = rnd_;
        
        // initialize sigma
        for (int i = 0; i < N; i++) {
            sigma[i] = 1.0;
        }
        
        // initialize alpha
        for (int i = 0; i < K; i++) {
            alpha[i] = 0.0;
        }
        
        // initialize x
        for (int i = 0; i < N; i++) {
            x[i] = rnd_.nextDouble() * 10.0 - 5.0;
        }
    }
    
    @Override
    public int compareTo(Individual other) {
        return (int) Math.signum(other.fitness - this.fitness);
    }
    
    private void printMatrix(double[][] A) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                System.out.print(A[i][j]);
                System.out.print(" ");
            }
            System.out.println();
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
        //System.out.print("Difference: ");
        //System.out.print(diff);
        //System.out.println();
    }
    
    public void mutate() {
        // mutate sigma
        for (int i = 0; i < N; i++) {
            sigma[i] *= Math.exp(tau_prime * rnd_.nextGaussian() + tau * rnd_.nextGaussian());
            if (sigma[i] < epsilon) {
                sigma[i] = epsilon;
            }
        }

        // mutate alpha
        for (int i = 0; i < K; i++) {
            //alpha[i] += beta * rnd_.nextGaussian();
            if (Math.PI < Math.abs(alpha[i])) {
                alpha[i] -= 2 * Math.PI * Math.signum(alpha[i]);
            }
        }
        
        // mutate x
        double[] v = new double[N];
        for (int i = 0; i < N; i++) {
            v[i] = rnd_.nextGaussian();
        }
        double[][] C = covarianceMatrix();
        double[][] L = cholesky(C);
        check(C, L);
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
    
    public void crossover(Individual a, Individual b) {
        // crossover sigma
        for (int i = 0; i < N; i++) {
            double r = rnd_.nextDouble();
            sigma[i] = (r * a.sigma[i] + (1.0 - r) * b.sigma[i] ) / 2;
        }
        
        // crossover alpha
        for (int i = 0; i < K; i++) {
            double r = rnd_.nextDouble();
            alpha[i] = (r * a.alpha[i] + (1.0 - r) * b.alpha[i] ) / 2;
        }
        
        // crossover x
        for (int i = 0; i < N; i++) {
            double r = rnd_.nextDouble();
            x[i]     = (r * a.x[i] + (1.0 - r) * b.x[i] ) / 2;
        }
    }
}
