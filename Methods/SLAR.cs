using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Methods
{
    public class SLAR
    {
        public double[] gaus(double[,] matrix, double[] vector)
        {
            double temp;
            double max;
            double m;
            int maxp;
            int n = vector.Length;
            double[] x = new double[n];
            for (int k = 0; k < n - 1; ++k)
            {
                max = matrix[k, k];
                maxp = k;
                for (int i = k + 1; i < n; ++i)
                {
                    if (matrix[i, k] > max)
                    {
                        max = matrix[i, k];
                        maxp = i;
                    }
                }
                if (maxp != k)
                {
                    for (int c = k; c < n; c++)
                    {
                        temp = matrix[maxp, c]; 
                        matrix[maxp, c] = matrix[k, c];
                        matrix[k, c] = temp;
                    }
                    temp = vector[maxp];
                    vector[maxp] = vector[k];
                    vector[k] = temp;
                }
                for (int i = k + 1; i < n; i++)
                {
                    if (matrix[k, k] == 0) matrix[k, k] = 0.0000000000001;
                    m = -(double)(matrix[i, k] / matrix[k, k]);
                    vector[i] = vector[i] + m * vector[k];
                    for (int j = k; j < n; j++)
                    {
                        matrix[i, j] = matrix[i, j] + m * matrix[k, j];
                    }
                }
            }
            double suma = 0;
            if (matrix[n-1, n-1] == 0) matrix[n-1, n-1] = 0.0000000000001;
            x[n - 1] = (double)(vector[n - 1] / matrix[n - 1, n - 1]);
            for (int k = n - 2; k >= 0; k--)
            {
                suma = 0;
                for (int j = k + 1; j < n; j++)
                {
                    suma = suma + matrix[k, j] * x[j];
                }
                if (matrix[k, k] == 0) matrix[k, k] = 0.0000000000001;
                x[k] = (double)((vector[k] - suma) / matrix[k, k]);
            }
            return x;
        }
        public void LU_matrix(double[,] matrix, out double[,] l, out double[,] u)
        {
            int n = matrix.GetLength(0);
            l = new double[n, n];
            u = new double[n, n];
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    u[i, j] = 0;
                    l[i, j] = 0;
                    if (i == j)
                    {
                        l[i, j] = 1;
                    }
                }
            }
            double suma;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (j >= i)
                    {
                        suma = 0;
                        for (int k = 0; k < i; k++)
                        {
                            suma = suma + l[i, k] * u[k, j];
                        }

                        u[i, j] = (double)(matrix[i, j] - suma);
                    }
                    if (j < i)
                    {
                        suma = 0;
                        for (int k = 0; k < j; k++)
                        {
                            suma = suma + l[i, k] * u[k, j];
                        }

                        l[i, j] = (double)((matrix[i, j] - suma) / u[j, j]);
                    }
                }
            }
        }
        public double[] LU(double[,] matrix, double[] vector)
        {
            int n = vector.Length;
            double[,] l = new double[n, n];
            double[,] u = new double[n, n];
            LU_matrix(matrix, out l, out u);
            double suma;
            double[] y = new double[n];
            double[] x = new double[n];
            for (int i = 0; i < n; i++)
            {
                y[i] = 0;
            }
            for (int i = 0; i < n; ++i)
            {
                suma = 0;
                for (int k = 0; k < i; k++)
                {
                    suma = suma + l[i, k] * y[k];
                }
                y[i] = vector[i] - suma;

            }
            for (int i = n - 1; i >= 0; --i)
            {
                suma = 0;
                for (int k = i + 1; k < n; k++)
                {
                    suma = suma + u[i, k] * x[k];
                }
                x[i] = (y[i] - suma) / u[i, i];

            }
            return x;
        }
       
    }
    public class prohonka
    {
        int n;
        double[] a;
        double[] b;
        double[] c;
        double[] alpha;
        double[] beta;
        double[] ksi;
        double[] eta;
        double[] y;

        public prohonka(double[,] a, double[] b)
        {
            n = b.Length;
            this.a = new double[n];
            this.b = new double[n];
            c = new double[n];
            alpha = new double[n];
            beta = new double[n];
            ksi = new double[n];
            eta = new double[n];
            y = new double[n];
            for (int i = 0; i < n - 1; i++)
            {
                this.b[i] = a[i, i + 1];
            }
            for (int i = 0; i < n; i++)
            {
                c[i] = a[i, i];
                y[i] = b[i];
            }
            for (int i = 1; i < n; i++)
            {
                this.a[i] = a[i, i - 1];
            }
        }

        public double[] getX()
        {
            double[] x = new double[n];
            int p;
            p = (int)Math.Truncate((double)n / 2);
            alpha[1] = -b[0] / c[0];
            beta[1] = y[0] / c[0];
            for (int i = 1; i < p; i++)
            {
                alpha[i + 1] = -b[i] / (a[i] * alpha[i] + c[i]);
                beta[i + 1] = (y[i] - a[i] * beta[i]) / (a[i] * alpha[i] + c[i]);
            }
            ksi[n - 1] = -a[n - 1] / c[n - 1];
            eta[n - 1] = y[n - 1] / c[n - 1];
            for (int i = n - 2; i >= p - 1; i--)
            {
                ksi[i] = -a[i] / (c[i] + b[i] * ksi[i + 1]);
                eta[i] = (y[i] - b[i] * eta[i + 1]) / (c[i] + b[i] * ksi[i + 1]);
            }
            x[p - 1] = (alpha[p] * eta[p] + beta[p]) / (1 - alpha[p] * ksi[p]);
            for (int i = p - 2; i >= 0; i--)
            {
                x[i] = alpha[i + 1] * x[i + 1] + beta[i + 1];
            }
            x[p] = (y[p] - b[p] * eta[p + 1]) / (b[p] * ksi[p + 1] + c[p]);
            for (int i = p - 1; i < n - 1; i++)
            {
                x[i + 1] = ksi[i + 1] * x[i] + eta[i + 1];
            }
            return x;
        }
    }
}
