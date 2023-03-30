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
        string file_name = "C:\\tso results\\Tuna_swarm_optimization";

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

        void limit_setter()
        {
            string user_input;
            double limit_d, limit_u;
            Console.WriteLine("Czy granice dla wszystkich wymiarow sa identyczne? [t/n]");
            char ans;
            ans = Convert.ToChar(Console.ReadLine());

            if (ans == 't' || ans == 'T')
            {
                
                Console.WriteLine("podaj ograniczenie dolne");
                user_input = Console.ReadLine();
                limit_d = Double.Parse(user_input);


                Console.WriteLine("podaj ograniczenie gorne");
                user_input = Console.ReadLine();
                limit_u = Double.Parse(user_input);

                for (int i = 0; i < dimension; i++)
                {
                    upper_limit[i] = limit_u;
                    lower_limit[i] = limit_d;
                }
            }

            else
            {
                for (int i = 0; i < dimension; i++)
                {
                    Console.WriteLine("Wymiar: " + (i + 1));
                    Console.WriteLine("podaj ograniczenie dolne");
                    user_input = Console.ReadLine();
                    limit_d = Double.Parse(user_input);

                    Console.WriteLine("podaj ograniczenie gorne");
                    user_input = Console.ReadLine();
                    limit_u = Double.Parse(user_input);
                    upper_limit[i] = limit_u;
                    lower_limit[i] = limit_d;
                }
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
            beta = inside_beta* Math.Exp(b * l)/100000000; //tu to jest tak chwilowo bo inaczej jakeis magiczne argumenty wyskakauja a tak sa chociaz w przedzialexD

            p = Math.Pow((1 - arg_for_p), arg_for_p);
            //Console.WriteLine("a1 = " + alpha_1 + " a2 = " + alpha_2 + " beta = " + beta + " l = " + l + " p = " + p + " arg = " + arg_for_p);

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
            //results[i] = f(args[i]);
            //number_of_calls++;
        }
        void equation_2_transformation(double alpha_1, double alpha_2, double beta, int i)
        {
            double[] new_args = new double[dimension];
            for (int j=0; j<dimension; j++)
            {
                double diff = Math.Abs(best_arguments[j] - arguments[i][j]);
                Console.WriteLine("ROZNICA: " + diff);
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
            //results[i] = f(arguments[i]);
            //number_of_calls++;
            //update_best(arguments[i], results[i]);
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
                Console.WriteLine("ROZNICA: " + diff);
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
            results[i] = f(arguments[i]);
            number_of_calls++;
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


            if (File.Exists(file_name + ".txt"))
            {
                StreamReader sr = new StreamReader(file_name + ".txt");
                string line = "";

                line = sr.ReadLine();
                current_iteration = Convert.ToInt32(line);
                line = sr.ReadLine();
                number_of_calls = Convert.ToInt32(line);

                for (int i = 0; i < numb_of_population; i++)
                {
                    line = sr.ReadLine();
                    string[] numbers = line.Split(", ");
                    results[i] = Double.Parse(numbers[0]);

                    for(int j=0; j< dimension; j++)
                    {
                        arguments[i][j] = Double.Parse(numbers[j+1]);
                    }
                }
                sr.Close();
            }
 
        }
        public void SaveResult()
        {
            StreamWriter sw = File.CreateText(file_name + "_END_RESULT.txt");
            sw.WriteLine("Ilosc wywolan funkcji: " + number_of_calls);
            sw.WriteLine("Rozmiar populacji: " + numb_of_population);
            sw.WriteLine("Ilosc iteracji: " + number_of_iterations);
            sw.WriteLine("parametr a: " + const_a);
            sw.WriteLine("parametr z: " + const_z);
            sw.WriteLine("Ilosc wymiarow funkcji: " + dimension + '\n' );

            sw.Write("Najlepszy wynik: " + best_result + " jego argumenty: ");
            for(int i=0; i<dimension; i++)
            {
                sw.Write(best_arguments[i]+", ");
                sw.Write('\n');
            }



        }
        public void SaveToFileStateOfAlghoritm()
        {
            
              StreamWriter  sw = File.CreateText(file_name+".txt");

            sw.WriteLine(current_iteration);
            sw.WriteLine(number_of_calls);

            for (int i = 0; i < numb_of_population; i++)
            {
                sw.Write(results[i] + ", ");
                for (int j = 0; j < dimension; j++)
                {
                    sw.Write(arguments[i][j]+", ");
                }
                sw.Write('\n');
            }
            sw.Close();
        }
        public double Solve()
        {
            limit_setter();
            //current_iteration = 0;
            Random random = new Random();
            LoadFromFileStateOfAlghoritm();
            if (current_iteration == 0) { creating_initial_population();}
            update_best();
           


            for (int w = current_iteration; w < number_of_iterations; w++)
            {
                displaying_in_console(w);
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
            displaying_in_console(number_of_iterations);
            SaveResult();
            return best_result;
        }
    }
}
