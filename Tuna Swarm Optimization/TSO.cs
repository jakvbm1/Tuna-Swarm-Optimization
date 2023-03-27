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

        public TSO(int number_of_iterations, int numb_of_population, int dimension, tested_function f, double const_z, double const_a)
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
                equation_1(arguments, i);
                results[i] = f(arguments[i]);
                number_of_calls++;
            }
            current_iteration++;

            SaveToFileStateOfAlghoritm();
        }
        void calculate_parameters(ref double alpha_1, ref double alpha_2,ref double beta, ref double l, ref double p)
        {
            Random random= new Random();
            double b = random.NextDouble();
            alpha_1 = const_a + (1 - const_a) * ((double)current_iteration / number_of_iterations);
            alpha_2 = (1-const_a) - (1-const_a)*((double)current_iteration / number_of_iterations);
            l = Math.Exp(3 * Math.Cos((number_of_iterations + (1 / current_iteration)) + 1) * Math.PI);
            beta = Math.Exp(b * l) * Math.Cos(2 * Math.PI * b);

            double arg_for_p = current_iteration/ number_of_iterations;
            p = Math.Pow((1 - arg_for_p), arg_for_p);

        }
        void update_best()
        {

            best_result = results[0];
            for (int j = 0; j < dimension; j++)
            {
                best_arguments[j] = arguments[0][j];
            }
            for (int i=0; i<numb_of_population; i++)
            {
                if (results[i]<best_result)
                {
                    best_result = results[i];
                    for(int j=0; j<dimension; j++)
                    {
                        best_arguments[j] = arguments[i][j];
                    }
                }
            }

        }
        void equation_1(double[][]args, int i) //to równanie jako jedyne przyjmuje też tablice argumentów, bo będę go używać zarówno do wypełniania arguments[][] jak i temp_arguments[][]
        {
            Random random = new Random();
            for(int j=0; j<dimension; j++)
            {
                args[i][j] = random.NextDouble() * (upper_limit[j] - lower_limit[j]) + lower_limit[j];
            }
            //results[i] = f(args[i]);
            //number_of_calls++;
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
                temp_arguments[i][j] = new_args[j];
            }
            //results[i] = f(arguments[i]);
            //number_of_calls++;
            //update_best(arguments[i], results[i]);
        }
        void equation_7_transformation(double alpha_1, double alpha_2, double beta, int i)
        {
            double[] random_point = new double[dimension];
            for(int j=0; j<dimension; j++)
            {
                Random random = new Random();
                random_point[j] = random.NextDouble() * (upper_limit[j] - lower_limit[j]) + lower_limit[j];
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
                temp_arguments[i][j] = new_args[j];
            }
            //results[i] = f(arguments[i]);
            //number_of_calls++;
            //update_best(arguments[i], results[i]);
        }
        void equation_9_transformation(double rand, double p, int i)
        {
            Random random = new Random();
            //int neg_or_pos = random.Next(1, 2);
            //if (neg_or_pos == 2) { neg_or_pos = -1;}

            //double TF = neg_or_pos * random.NextDouble();
            double[] new_args = new double[dimension];
            for (int j = 0; j < dimension; j++)
            {
                //nie jestem pewny czy TF ma być generowane dla każdego wymiaru czy jedno dla wszystkich więc na razie ustawilem osobno dla kazdego wymiaru
                int neg_or_pos = random.Next(1, 2);
                if (neg_or_pos == 2) { neg_or_pos = -1; }

                double TF = neg_or_pos * random.NextDouble();

                if (rand < 0.5)
                {
                    new_args[j] = best_arguments[j] + rand * (best_arguments[j] - arguments[i][j]) + TF * p * p * (best_arguments[j] - arguments[i][j]);
                }

                else
                {
                    new_args[j] = TF * p * p * arguments[i][j];
                }
            }

                for (int j = 0; j < dimension; j++)
            {
                temp_arguments[i][j] = new_args[j];
            }
            //results[i] = f(arguments[i]);
            //number_of_calls++;
            //update_best(arguments[i], results[i]);
        }

        public int NumberOfEvaluationFitnessFunction => number_of_calls;
        double[] IOptimizationAlgorithm.XBest => best_arguments;
        double IOptimizationAlgorithm.FBest => best_result;
        public void LoadFromFileStateOfAlghoritm()
        {
            StreamReader sr = new StreamReader(file_name);
            string line = "";

            if (File.Exists(file_name))
            {
                line = sr.ReadLine();
                current_iteration = int.Parse(line);
                line = sr.ReadLine();
                number_of_calls = int.Parse(line);

                for (int i = 0; i < numb_of_population; i++)
                {
                    line = sr.ReadLine();
                    string[] numbers = line.Split(" ");
                    results[i] = int.Parse(numbers[0]);

                    for(int j=0; j< dimension; j++)
                    {
                        arguments[i][j] = int.Parse(numbers[j+1]);
                    }
                }
            }
            sr.Close();
        }
        public void SaveResult()
        {
            StreamWriter sw = new StreamWriter(file_name + "_END_RESULT");
            sw.WriteLine("Ilosc wywolan funkcji: " + number_of_calls + '\n');
            sw.WriteLine("Rozmiar populacji: " + numb_of_population + '\n');
            sw.WriteLine("Ilosc iteracji: " + number_of_iterations + '\n');
            sw.WriteLine("Ilosc wymiarow funkcji: " + dimension + '\n');



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
            sw.Close();
        }
        public double Solve()
        {
            current_iteration = 0;
            Random random = new Random();
            LoadFromFileStateOfAlghoritm();
            if (current_iteration == 0) { creating_initial_population();}
            update_best();
            


            for(int w = current_iteration; w < number_of_iterations; w++)
            {
                for(int i=0; i<numb_of_population; i++)
                {
                    double alpha_1=0, alpha_2 = 0, beta = 0, l = 0, p = 0;
                    calculate_parameters(ref alpha_1, ref alpha_2, ref beta, ref l, ref p);
                    double rand = random.NextDouble();

                    if (rand > const_z)
                    {
                        equation_1(temp_arguments, i);
                    }

                    else if(rand >= 0.5)
                    {
                        equation_9_transformation(rand, p, i);
                    }

                    else if((current_iteration/number_of_iterations) < rand)
                    {
                        equation_7_transformation(alpha_1, alpha_2, beta, i);
                    }

                    else
                    {
                        equation_2_transformation(alpha_1, alpha_2, beta, i);
                    }
                }
                for (int i = 0; i < numb_of_population; i++)
               { 
                for (int j = 0; j<dimension; j++)
                {
                    arguments[i][j] = temp_arguments[i][j];
                }
                results[i] = f(arguments[i]);
                number_of_calls++;
               }

                update_best();
                SaveToFileStateOfAlghoritm();
                current_iteration++;
            }

            return best_result;
        }
    }
}
