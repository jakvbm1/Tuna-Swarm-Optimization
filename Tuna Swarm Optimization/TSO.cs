using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Transactions;

namespace Tuna_Swarm_Optimization
{
    internal class TSO : IOptimizationAlgorithm
    {
        //parametry algorytmu
       private int number_of_iterations;
       private int numb_of_population;
       private int dimension;
       private double const_z;
       private double const_a;
       private double[] upper_limit;
       private double[] lower_limit;
       private double[] probabilities;
       private double[] index;

        //rozwiazania
       private double[][] arguments;
       private double[] results;
       private double[][] temp_arguments;
       private double[] temp_results;

        //najlepsze rozwiazanie
        public double[] best_arguments;
        public double best_result;

        //zmienne pomocnicze
        public int number_of_calls = 0;
        public int current_iteration = 0;
        string file_name = "Tuna_swarm_optimization.txt";

        public delegate double tested_function(params double[] arg);
        private tested_function f;

        public TSO(int number_of_iterations, int numb_of_population, int dimension, tested_function f, double const_z, double const_a,)
        {
            this.number_of_iterations = number_of_iterations;
            this.numb_of_population = numb_of_population;
            this.dimension = dimension;
            this.const_z = const_z;
            this.const_a = const_a;
            this.upper_limit = new double[dimension];
            this.lower_limit = new double[dimension];
            this.arguments = new double[numb_of_population][];
            this.results = new double[numb_of_population];
            this.temp_arguments = new double[numb_of_population][];
            this.temp_results = new double[numb_of_population];
            best_arguments = new double[dimension];
            probabilities = new double[numb_of_population];
            index = new double[numb_of_population];
            best_result = 0;
            this.f = f;

            for(int i=0; i<numb_of_population; i++)
            {
                arguments[i] = new double[dimension];
                temp_arguments[i] = new double[dimension];
            }
        }
        void creating_initial_population()
        {

            for (int i=0; i<numb_of_population; i++)
            {
                equation_1(i);
            }
            current_iteration++;
            Sort_population();
            SaveToFileStateOfAlghoritm();
        }
        void Sort_population()
        {
            for(int i=1; i<numb_of_population; i++)
            {
                for (int j=0; j<numb_of_population-i; j++)
                {
                    if (results[j] > results[j+1])
                    {
                        double result_placeholder;
                        double[] arguments_placeholder = new double[dimension];
                        result_placeholder = results[j];
                        for (int k=0; k<dimension; k++)
                        {
                            arguments_placeholder[k] = arguments[j][k];
                        }

                        results[j] = results[j+1];
                        results[j+1] = result_placeholder;

                        for(int k=0; k<dimension; k++)
                        {
                            arguments[j][k] = arguments[j + 1][k];
                            arguments[j + 1][k] = arguments_placeholder[k];
                        }
                    }
                }
            }

            best_result = results[0];
            for(int i=0; i<dimension; i++)
            {
                best_arguments[i] = arguments[0][i];
            }
        }
        void calculate_parameters(double &alpha_1, double &alpha_2,double &beta, double &l, double &p)
        {
            Random random= new Random();
            double b = random.NextDouble();
            alpha_1 = const_a + (1 - cons_a) * ((double)current_iteration / number_of_iterations);
            alpha_2 = (1-const_a) - (1-const_a)*((double)current_iteration / number_of_iterations);
            l = Math.Exp(3 * Math.Cos((number_of_iterations + (1 / current_iteration)) + 1) * Math.PI);
            beta = Math.Exp(b * l) * Math.Cos(2 * Math.PI * b);

            double arg_for_p = current_iteration/ number_of_iterations;
            p = Math.Pow((1 - arg_for_p), arg_for_p);

        }
        void update_best(double[] args, double res)
        {
            if (res < best_result)
            {
                best_result = res;
                for(int i=0; i<dimension; i++)
                {best_arguments[i] = args[i];}
            }
        }

        void equation_1(int i)
        {
            Random random = new Random();
            for(int j=0; j<dimension; j++)
            {
                arguments[i][j] = random.NextDouble() * (upper_limit[j] - lower_limit[j]) + lower_limit[j];
            }
            results[i] = f(arguments[i]);
            number_of_calls++;
        }



        void equation_2_transformation(double alpha_1, double alpha_2, double beta, int i)
        {
            double[] new_args = new double[dimension];
            for (int j=0; j<dimension; j++)
            {
                if (i == 0)
                {
                    new_args[j] = alpha_1* (best_arguments[j] + beta*Math.Abs(best_arguments[j] - arguments[i][j])) + alpha_2*arguments[i][j];
                }
                else
                {
                    new_args[j] = alpha_1 * (best_arguments[j] + beta * Math.Abs(best_arguments[j] - arguments[i][j])) + alpha_2 * arguments[i-1][j];
                }
            }
            
            for(int j=0; j<dimension; j++)
            {
                arguments[i][j] = new_args[j];
            }
            results[i] = f(arguments[i]);
            number_of_calls++;
            update_best(arguments[i], results[i]);
        }
        void equation_7_transformation(double alpha_1, double alpha_2, double beta, int i)
        {
            double[] random_point = new double[dimension];
            for(int i=0; i<dimension; i++)
            {
                Random random = new Random();
                random_point[i] = random.NextDouble() * (upper_limit[j] - lower_limit[j]) + lower_limit[j];
            }

            double[] new_args = new double[dimension];
            for (int j = 0; j < dimension; j++)
            {
                if (i == 0)
                {
                    new_args[j] = alpha_1 * (random_point[j] + beta * Math.Abs(random_point[j] - arguments[i][j])) + alpha_2 * arguments[i][j];
                }
                else
                {
                    new_args[j] = alpha_1 * (best_arguments[j] + beta * Math.Abs(random_point[j] - arguments[i][j])) + alpha_2 * arguments[i - 1][j];
                }
            }

            for (int j = 0; j < dimension; j++)
            {
                arguments[i][j] = new_args[j];
            }
            results[i] = f(arguments[i]);
            number_of_calls++;
            update_best(arguments[i], results[i]);
        }




        public int NumberOfEvaluationFitnessFunction => throw new NotImplementedException();

        double[] IOptimizationAlgorithm.XBest => throw new NotImplementedException();

        double IOptimizationAlgorithm.FBest => throw new NotImplementedException();

        public void LoadFromFileStateOfAlghoritm()
        {
            throw new NotImplementedException();
        }
        public void SaveResult()
        {
        }
        public void SaveToFileStateOfAlghoritm()
        {
            StreamWriter sw = new StreamWriter(file_name);
            sw.WriteLine(current_iteration);
            sw.WriteLine(number_of_calls);

            for (int i = 0; i < numb_of_population; i++)
            {
                for (int j = 0; j < dimension; j++)
                {
                    sw.Write(arguments[i][j]);
                }
                sw.Write(results[i] + '\n');
            }
        }
        public double Solve()
        {
            LoadFromFileStateOfAlghoritm();
            if (current_iteration== 0) { creating_initial_population();}

            for(current_iteration; current_iteration<number_of_iterations; current_iteration++)
            {
                for(int i=0; i<numb_of_population; i++)
                {
                    double alpha_1, alpha_2, beta, l, p;
                    calculate_parameters(alpha_1, alpha_2, beta, l, p);

                }
            }


        }
    }
}
