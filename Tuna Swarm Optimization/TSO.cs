using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

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

        public TSO(int number_of_iterations, int numb_of_population, int dimension, double const_z,double const_a, tested_function f)
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
            Random random = new Random();
            for (int i=0; i<numb_of_population; i++)
            {
                for (int j=0; j<dimension; j++)
                {
                    arguments[i][j] = random.NextDouble() * (upper_limit[j] - lower_limit[j]) + lower_limit[j];
                }
                results[i] = f(arguments[i]);
                number_of_calls++;
            }
            current_iteration++;
            SaveToFileStateOfAlghoritm();
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
            throw new NotImplementedException();
        }
    }
}
