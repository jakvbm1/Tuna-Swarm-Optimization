using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Globalization;
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

        //rozwiazania
       private double[][] arguments;
       private double[] results;
       private double[][] temp_arguments;

        //najlepsze rozwiazanie
        public double[] best_arguments;
        public double best_result;

        //zmienne pomocnicze
        public int number_of_calls = 0;
        public int current_iteration = 0;

        string file_name = "Tuna_swarm_optimization.txt";

        public delegate double tested_function(params double[] arg);
        private tested_function f;

        public TSO(int number_of_iterations, int numb_of_population, int dimension, tested_function f, double const_z=0.05, double const_a=0.7)
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
            best_arguments = new double[dimension];
            best_result = 1000000;
            this.f = f;

            for(int i=0; i<numb_of_population; i++)
            {
                arguments[i] = new double[dimension];
                temp_arguments[i] = new double[dimension];
            }
        }

       public void limit_setter(double[] tab_a, double[] tab_b)
        {
            for (int i=0; i<dimension; i++)
            {
                upper_limit[i] = tab_a[i];
                lower_limit[i] = tab_b[i];
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
        void calculate_parameters(ref double alpha_1, ref double alpha_2,ref double beta, ref double l, ref double p, int iteration)
        {
            Random random= new Random();
            double b = random.NextDouble();
            double arg_for_p = (double)iteration / number_of_iterations;
            alpha_1 = const_a + ((1 - const_a) * arg_for_p);
            alpha_2 = (1-const_a) - ((1-const_a) * arg_for_p);


            double ins_ins_l = 1.0 / iteration;
            double inside_l = Math.Cos(((number_of_iterations + ins_ins_l) - 1) * Math.PI);
            double insx3l = 3 * inside_l;
            l = Math.Exp(insx3l);
            double inside_beta = Math.Cos(2 * Math.PI * b);
            beta = inside_beta * Math.Exp((b * l)); 

            p = Math.Pow((1 - arg_for_p), arg_for_p);
           

        }
        void update_best()
        {

           // best_result = results[0];
            //for (int j = 0; j < dimension; j++)
            //{
            //    best_arguments[j] = arguments[0][j];
            //}
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
        }
        void equation_2_transformation(double alpha_1, double alpha_2, double beta, int i)
        {
            double[] new_args = new double[dimension];
            for (int j=0; j<dimension; j++)
            {
                double diff = Math.Abs(best_arguments[j] - arguments[i][j]);
                double multiply_by_a = best_arguments[j] + (beta * diff);

                if (i == 0)
                {
                    new_args[j] = (alpha_1* multiply_by_a) + (alpha_2*arguments[i][j]);
                }
                else
                {
                    new_args[j] = (alpha_1 * multiply_by_a) + (alpha_2 * arguments[i-1][j]);
                }
            }
            
            for(int j=0; j<dimension; j++)
            {
                temp_arguments[i][j] = new_args[j];
            }
        }
        void equation_7_transformation(double alpha_1, double alpha_2, double beta, int i)
        {
            double[] random_point = new double[dimension];
            Random random = new Random();
            for (int j=0; j<dimension; j++)
            {
                random_point[j] = random.NextDouble() * (upper_limit[j] - lower_limit[j]) + lower_limit[j];
            }

            double[] new_args = new double[dimension];
            for (int j = 0; j < dimension; j++)
            {
                double diff = Math.Abs(random_point[j] - arguments[i][j]);
                double multiply_by_a = random_point[j] + (beta * diff);
                if (i == 0)
                {
                    new_args[j] = (alpha_1 * multiply_by_a) + (alpha_2 * arguments[i][j]);
                }
                else
                {
                    new_args[j] = (alpha_1 * multiply_by_a) + (alpha_2 * arguments[i - 1][j]);
                }
            }

            for (int j = 0; j < dimension; j++)
            {
                temp_arguments[i][j] = new_args[j];
            }
        }
        void equation_9_transformation(double rand, double p, int i)
        {
            Random random = new Random();

            double[] new_args = new double[dimension];
            for (int j = 0; j < dimension; j++)
            {
                
                int neg_or_pos = random.Next(1, 3);
                if (neg_or_pos == 2) { neg_or_pos = -1; }

                int TF = neg_or_pos;
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
        }

        void displaying_in_console(int iteration) //funkcja do testowania by "podgladac" wyniki dzialania programu w konsoli
        {
            Console.WriteLine("Number of iteration: " + iteration);
            Console.WriteLine("Results:");
            for (int i=0; i<numb_of_population; i++)
            {
                Console.Write(results[i] + ", ");
                for (int j=0; j<dimension; j++)
                {
                    Console.Write(arguments[i][j] + " ");
                }
                Console.Write('\n');
            }

            Console.WriteLine("BEST RESULT: " + best_result);
            Console.Write("BEST ARGS: ");
            for (int j = 0; j < dimension; j++)
            {
                Console.Write(best_arguments[j] + " ");
            }
        }

        public int NumberOfEvaluationFitnessFunction => number_of_calls;
        double[] IOptimizationAlgorithm.XBest => best_arguments;
        double IOptimizationAlgorithm.FBest => best_result;
        public void LoadFromFileStateOfAlghoritm()
        {


            if (File.Exists(file_name))
            {
                StreamReader sr = new StreamReader(file_name);
                string line = "";

                line = sr.ReadLine();
                current_iteration = Convert.ToInt32(line);
                line = sr.ReadLine();
                number_of_calls = Convert.ToInt32(line);

                for (int i = 0; i < numb_of_population; i++)
                {
                    line = sr.ReadLine();
                    string[] numbers = line.Split(", ");
                                      

                    for(int j=0; j< dimension; j++)
                    {
                        arguments[i][j] = Convert.ToDouble(numbers[j]);
                    }
                    results[i] = Convert.ToDouble(numbers[dimension]);
                }
                sr.Close();
            }
 
        }
        public void SaveResult()
        {
            StreamWriter sw = File.CreateText("TSO_END.txt");
            sw.WriteLine(number_of_calls);
            sw.WriteLine(numb_of_population);
            sw.WriteLine(number_of_iterations);
            sw.WriteLine(const_a);
            sw.WriteLine(const_z);
            sw.WriteLine(dimension);

            sw.WriteLine(best_result);
            for(int i=0; i<dimension; i++)
            {
                sw.Write(best_arguments[i]+", ");
                
            }
            sw.Write('\n');
            sw.Close();

        }
        public void SaveToFileStateOfAlghoritm()
        {
            
              StreamWriter  sw = File.CreateText(file_name);

            sw.WriteLine(current_iteration);
            sw.WriteLine(number_of_calls);

            for (int i = 0; i < numb_of_population; i++)
            {
                
                for (int j = 0; j < dimension; j++)
                {
                    sw.Write(arguments[i][j]+", ");
                }
                sw.Write(results[i]);
                sw.Write('\n');
            }
            sw.Close();
            string name = "TSO_iteracja_"+current_iteration.ToString() + ".txt";
            Console.WriteLine(name);
            StreamWriter sw2 = File.CreateText(name);

            sw2.WriteLine(const_a);
            sw2.WriteLine(const_z);
            sw2.WriteLine(dimension);
            sw2.WriteLine(number_of_iterations);
            sw2.WriteLine(numb_of_population);
            sw2.WriteLine(number_of_calls);

            //sortowanie aktualnych wynikow by w pliku byly zapisane od najlepszego do najgorszego, sortuje na przekopiowanych do tymczasowych zmiennych
            //zamiast "oryginalow", bo sortowanie na oryginalach wplyneloby na dzialanie argumentu
            double[] res_to_save = new double[numb_of_population];
            double[][] args_to_save = new double[numb_of_population][];

            for (int i = 0; i < numb_of_population; i++)
            {
                args_to_save[i] = new double[dimension];
            }

            for (int i=0; i< numb_of_population;i++)
            {
                res_to_save[i] = results[i];
                for (int j=0; j < dimension; j++)
                {
                    args_to_save[i][j] = arguments[i][j];
                }
            }

            for(int i=0; i<numb_of_population; i++)
            {
                for (int j=0; j<numb_of_population-i-1; j++)
                {
                    if (res_to_save[j] > res_to_save[j+1])
                    {
                        double temp_r = res_to_save[j];

                        res_to_save[j] = res_to_save[j + 1];
                        res_to_save[j+1] = temp_r;

                        for (int k=0; k<dimension; k++)
                        {
                            double temp_a = args_to_save[j][k];
                            args_to_save[j][k] = args_to_save[j+1][k];
                            args_to_save[j+1][k] = temp_a;
                        }
                    }
                }
            }

            for (int i=0; i<numb_of_population; i++)
            {
                for (int j=0; j<dimension; j++)
                {
                    sw2.Write(args_to_save[i][j]+ ", ");
                }
                sw2.Write(res_to_save[i] + '\n');
            }
            sw2.Close();

        }
        public double Solve()
        {
            //current_iteration = 0;
            Random random = new Random();
            LoadFromFileStateOfAlghoritm();
            if (current_iteration == 0) { creating_initial_population();}
            update_best();
           
            //aby widziec efekty pracy "na biezaco" w konsoli nalezy odkomentowac komendy displaying_in_console

            for (int w = current_iteration; w < number_of_iterations; w++)
            {
                //displaying_in_console(w);
                for(int i=0; i<numb_of_population; i++)
                {
                    double alpha_1=0, alpha_2 = 0, beta = 0, l = 0, p = 0;
                    calculate_parameters(ref alpha_1, ref alpha_2, ref beta, ref l, ref p, w);
                    double rand = random.NextDouble();

                    if (rand < const_z)
                    {
                        
                        equation_1(temp_arguments, i);
                    }

                    else if(rand >= 0.5)
                    {
                        
                        equation_9_transformation(rand, p, i);
                    }

                    else if(((double)w/number_of_iterations) < rand)
                    {
                        equation_7_transformation(alpha_1, alpha_2, beta, i);
                    }

                    else
                    {

                       equation_2_transformation(alpha_1, alpha_2, beta, i);

                    }

                    current_iteration = w;
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
                current_iteration = w;
                update_best();
                SaveToFileStateOfAlghoritm();
            }
            //displaying_in_console(number_of_iterations);
            SaveResult();
            return best_result;
        }
    }
}
