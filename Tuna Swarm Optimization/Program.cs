
using System.Runtime.CompilerServices;

namespace Tuna_Swarm_Optimization
{
    class Program
    {
        public static double Rastrigin_function(params double[] x) // minimum (0,0, ..., 0) 0, df: [-5.12; 5.12]
        {
            int n = x.Length;
            double sum = 0;
            for (int i = 0; i < n; i++)
                sum += x[i] * x[i] - 10 * Math.Cos(2 * Math.PI * x[i]);
            return 10 * n + sum;
        }

        public static double Rosenbrock_function(params double[] x) //minimum (1, 1, ..., 1) 0 df: caly przedzial R
        {
            int dimension = x.Length;
            double sum = 0;
            if (dimension != 1)
            {
                for (int i = 0; i < dimension - 1; i++)
                {
                    sum += 100 * Math.Pow((x[i + 1] - (x[i] * x[i])), 2) + Math.Pow((1 - x[i]), 2);
                }
            }
            return sum;
        }

        public static double Sphere_function(params double[] x)// minimum (0,0, ..., 0) 0, df: caly przedzial R
        {
            int dimension = x.Length;
            double sum = 0;
            for (int i=0; i<dimension; i++)
            {
                sum += x[i] * x[i];
            }
            return sum;
        }

        public static double Beale_function(params double[] x) //minimum (3, 0.5) 0, df:[-4,5; 4,5] tylko dwa wymiary!
        {
            int dimension = x.Length;
            double sum = 0;
            if (dimension == 2)
            {
                sum = Math.Pow((1.5 - x[0] + (x[0] * x[1])), 2) + Math.Pow((2.25 - x[0] + (x[0] * x[1] * x[1])), 2)
                    + Math.Pow((2.625 - x[0] + (x[0] * x[1] * x[1] * x[1])), 2);

                return sum;
            }
            else
            {
                return 10000000;
            }
        }



    }

}
