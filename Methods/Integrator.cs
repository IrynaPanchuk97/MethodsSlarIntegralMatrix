using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
namespace Methods
{
   
    public class Integrator
    {
        public delegate double func(double x);
        public delegate double func1(double x, double z);
        public double square;
        private int No;
        private int N;
        private double rect_integral(double a, double b, int n, func f)
        {
            square = 0;
            double h = (double)((b - a) / n);
            for (int i = 0; i < n; ++i)
            {
                square += h * f(a + (1 / 2 + i) * h);
            }
            return square;
        }

        public double rectangle(double a, double b, double eps, func f)
        {
            double I=0;
            double Io=0;
            No=100;
            Io=rect_integral(a, b, No, f);
            N = 2 * No;
            I = rect_integral(a, b, N, f);
            while ((double)(Math.Abs(I - Io) / Math.Abs(I)) > eps)
            {
                No = N;
                N = 2 * No;
                Io = rect_integral(a, b, No, f);
                I = rect_integral(a, b, N, f);
            }
            return square;
        }
        public double TrapIntegral(double a, double b, int n, func f)
        {
            double h = (b - a) / n, tmp = 0, val = a;
            double[] x = new double[n + 1];
            for (int i = 0; i < x.Length; i++)
            {
                x[i] = val;
                if ((x[i] == a) || (x[i] == b))
                {
                    tmp = tmp + 1 / 2.0 * f(x[i]);
                }
                else
                {
                    tmp = tmp + f(x[i]);
                }
                val = val + h;
            }
            return (h * tmp);
        }
        public double trapeze(double a, double b, double eps, func f)
        {

            double I = 0;
            double Io = 0;
            No = 100;
            Io = trapeze_integral(a, b, No, f);
            N = 2 * No;
            I = trapeze_integral(a, b, N, f);
            while ((double)(Math.Abs(I - Io) / Math.Abs(I)) > eps)
            {
                No = N;
                N = 2 * No;
                Io = trapeze_integral(a, b, No, f);
                I = trapeze_integral(a, b, N, f);
            } 
            return square;
        }
        private double trapeze_integral(double a, double b, int n, func f)
        {
            double h = (double)( (b - a) / ( n) );
            double h1 =(double ) (h / 2);
            square = h * f(a);
            for (int i = 1; i < n-1; ++i)
            {
                square += h1 * 2 * f(a +  i * h);
            }
            square += h * f(b);
            return square;
        }

        public double simpson(double a, double b, double eps, func f)
        {

            double I = 0;
            double Io = 0;
            No = 100;
            Io = simpson_integral(a, b, No, f);
            N = 2 * No;
            I = simpson_integral(a, b, N, f);
            while ((double)(Math.Abs(I - Io) / Math.Abs(I)) > eps)
            {
                No = N;
                N = 2 * No;
                Io = simpson_integral(a, b, No, f);
                I = simpson_integral(a, b, N, f);
            }
            return square;
        }
        
        private double simpson_integral(double a, double b, int n, func f)
        {
            int m=n/2;
            int k = 0;
            double h = (double)((b - a) / n);
            double h1 = (double)((b - a) / (6*m));
            square = h * f(a);
            for (int i = 1; i < 2*m-1; ++i)
            {
                if (k%2==0)
                {
                    square += h1 * 4 * f(a + i * h);
                    
                }
                else
                {
                    square += h1 * 2 * f(a + i * h);
                   
                }
                k++;
            }
            square += h * f(b);
            return square;
        }
    }
}

