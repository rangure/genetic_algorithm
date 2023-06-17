#include <stdio.h>
#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <array>
#include "CSimulator.h"
#include "CCircuit.h"
#include "Genetic_Algorithm.h"
#include <omp.h>

// #define OUT_RESULT_TO_FILE 1
#ifdef PARALLEL_ENABLE
#include <mpi.h>
#endif
// #define PERFORMANCE_ANALYSIS_MACRO 2

/**
 * @brief Constructor of the GA class.
 * This constructor initializes the class member variables.
 *
 * @param AP A structure containing the algorithm parameters.
 * @param SP A structure containing the simulation parameters.
 */
GA::GA(Algorithm_Parameters AP, SimulationParameters SP)
{
   parameters.max_iterations = AP.max_iterations;
   parameters.tol = AP.tol;
   parameters.mutation_rate = AP.mutation_rate;
   parameters.mutation_size = AP.mutation_size;
   parameters.parent_pool_size = AP.parent_pool_size;
   parameters.population_size = AP.population_size;
   parameters.circuit_size = AP.circuit_size;
   parameters.selection_scheme = AP.selection_scheme;
   parameters.tournament_size = AP.tournament_size;
   parameters.crossover_rate = AP.crossover_rate;
   parameters.parallel_mpi = AP.parallel_mpi;
   parameters.communicate_interval = AP.communicate_interval;
   parameters.tournament_size_parallel = AP.tournament_size_parallel;
   parameters.parent_comm_size_parallel = AP.parent_comm_size_parallel;
   parameters.max_iter_without_progress = AP.max_iter_without_progress;
   parameters.mutation_rate_increase_factor = AP.mutation_rate_increase_factor;
   parameters.max_iter_before_fail = AP.max_iter_before_fail;
   parameters.write_interval = AP.write_interval;
   parameters.mutation_scheme = AP.mutation_scheme;

   sim_params.F0_ger = SP.F0_ger;     // Input flow of gerardium
   sim_params.F0_waste = SP.F0_waste; // Input flow of waste
   sim_params.k_ger = SP.k_ger;       // Rate constant for gerardium
   sim_params.k_waste = SP.k_waste;   // Rate constant for waste
   sim_params.volume = SP.volume;     // Volume of each cell
   sim_params.phi = SP.phi;
   sim_params.rho = SP.rho;                   // Density
   sim_params.ger_reward = SP.ger_reward;     // Reward for gerardium
   sim_params.waste_reward = SP.waste_reward; // Negative reward for waste
   sim_params.min_flowrate = SP.min_flowrate; // Minimum flowrate for the simulator
   sim_params.max_iter = SP.max_iter;         // Maximum number of iterations for simulator
   setup();
}

/**
 * @brief Default constructor of the GA class.
 * This constructor initializes the class member variables to default values.
 */
GA::GA()
{
   parameters.max_iterations = 1000;
   parameters.tol = 0.001;
   parameters.mutation_rate = 0.5;
   parameters.mutation_size = 3;
   parameters.parent_pool_size = 20;
   parameters.population_size = 100;
   parameters.circuit_size = 20;
   parameters.selection_scheme = 1;
   parameters.tournament_size = 5;
   parameters.communicate_interval = 10;
   parameters.crossover_rate = 0.9;
   parameters.max_iter_without_progress = 100;
   parameters.mutation_rate_increase_factor = 1.001;
   parameters.max_iter_before_fail = 20000;
   parameters.write_interval = 10;
   parameters.mutation_scheme = 1;
   parameters.parallel_mpi = 0;
   parameters.tournament_size_parallel = 4;
   parameters.parent_comm_size_parallel = parameters.parent_pool_size;

   sim_params.F0_ger = 10;     // Input flow of gerardium
   sim_params.F0_waste = 100;  // Input flow of waste
   sim_params.k_ger = 5.e-3;   // Rate constant for gerardium
   sim_params.k_waste = 5.e-4; // Rate constant for waste
   sim_params.volume = 10;     // Volume of each cell
   sim_params.phi = 0.1;
   sim_params.rho = 3.e3;            // Density
   sim_params.ger_reward = 100.0;    // Reward for gerardium
   sim_params.waste_reward = -500.0; // Negative reward for waste
   sim_params.min_flowrate = 1.e-10; // Minimum flowrate for the simulator
   sim_params.max_iter = 1000;       // Maximum number of iterations for simulator
   setup();
}

/**
 * @brief This function returns an index into a one-dimensional array given two indices. The 1d array contains n elements, each one has size of size_per_cir.
 *
 * @param pos_in_array The position in the array.
 * @param pos_in_vec The position in the vector.
 * @return An index in a one-dimensional array.
 */
int GA::get_idx(int pos_in_array, int pos_in_vec)
{
   return pos_in_array * size_per_cir + pos_in_vec;
}
// return the idx in 1d array
// pos=1 node = 10 val = 0 means it will return the index of the
// concentration out-node of 10th node of first parent/population
/**
 * @brief This function returns an index into a one-dimensional array given three indices. 
 *
 * @param pos The position in the array.
 * @param node The node index.
 * @param val The value index.
 * @return An index in a one-dimensional array.
 */
int GA::get_idx_node(int pos, int node, int val)
{
   return pos * size_per_cir + node * 2 + 1 + val;
}

/**
 * @brief This function sets up the genetic algorithm, initializing various parameters and data structures. The setup is based on the user passed parameters.
 */
void GA::setup()
{

   size_per_cir = parameters.circuit_size * 2 + 1;
   vec_size_parent = parameters.parent_pool_size * size_per_cir;
   vec_size_population = parameters.population_size * size_per_cir;
   file_write_buffer = NULL;
#ifdef OUT_RESULT_TO_FILE
   std::cout << "Output result to file functionality is enabled."<< std::endl;
   file_write_buffer = new double[size_per_cir - 1];
#endif
#ifdef PARALLEL_ENABLE
   std::cout << "parallel enable is set " << std::endl;
   if (parameters.parallel_mpi == 1)
   {
      std::cout << "Setting up MPI" << std::endl;
      setup_parallel();
   }
#endif
   mate_count_all = 0;
   same_mate = 0;
   nottouched = true;
   parents = new int[vec_size_parent];
   population = new int[vec_size_population];
   best_circuit = new int[size_per_cir];
   // -1 is a place holder as it is possible to have negative fitness
   best_fitness = -1;
   fitness_population = new double[parameters.population_size];
   fitness_parents = new double[parameters.parent_pool_size];
   if (parameters.selection_scheme == 2)
   {
      roulette_sum = new double[parameters.population_size];
   }
   else
   {
      roulette_sum = nullptr;
   }
   if (parameters.selection_scheme == 1)
   {
      rank_sum = new double[parameters.population_size];
      rank_val = new double[parameters.population_size];
      population_idex_list = new int[parameters.population_size];
   }
   else
   {
      rank_sum = nullptr;
      rank_val = nullptr;
      population_idex_list = nullptr;
   }
   srand(time(NULL));
   int flag = 0;
// Generating first batch of parents pool in parallel
#pragma omp parallel
   {
      CCircuit checker(size_per_cir);
#pragma omp for schedule(dynamic, 1)
      for (int i = 0; i < parameters.parent_pool_size; i++)
      {
         int starting_idx = get_idx(i, 0);
         int cur_idx;
         int local_iter_count_1 = 0;
         // loop until a valid circuit is found
         while (local_iter_count_1 < parameters.max_iter_before_fail)
         {
            local_iter_count_1++;
            // randomly generating each entry and enploying some of the rules to make it efficient
            parents[starting_idx] = rand() % parameters.circuit_size;
            for (int j = 0; j < parameters.circuit_size; j++)
            {
               cur_idx = get_idx_node(i, j, 0);
               if (j == starting_idx)
               {
                  flag = 0;
               }
               else
               {
                  flag = 2;
               }
               int local_iter_count_2 = 0;
               while (local_iter_count_2 < parameters.max_iter_before_fail)
               {
                  local_iter_count_2++;
                  parents[cur_idx] = rand() % (parameters.circuit_size + flag);
                  if (parents[cur_idx] != parameters.circuit_size + 1 && parents[cur_idx] != j)
                  {
                     break;
                  }
               }
               local_iter_count_2 = 0;
               while (local_iter_count_2 < parameters.max_iter_before_fail)
               {
                  local_iter_count_2++;
                  parents[cur_idx + 1] = rand() % (parameters.circuit_size + flag);
                  if (parents[cur_idx + 1] != parameters.circuit_size && parents[cur_idx + 1] != j)
                  {
                     break;
                  }
               }
            }
            // add concentration and tail to the array incase there is none
            int outpos1 = rand() % parameters.circuit_size;
            while (outpos1 == parents[starting_idx])
            {
               outpos1 = rand() % parameters.circuit_size;
            }
            parents[get_idx_node(i, outpos1, 1)] = parameters.circuit_size + 1;
            int outpos2 = rand() % parameters.circuit_size;
            while (outpos2 == parents[starting_idx])
            {
               outpos2 = rand() % parameters.circuit_size;
            }
            parents[get_idx_node(i, outpos2, 0)] = parameters.circuit_size;
            // check circuit validity
            if (checker.check_validity(&parents[starting_idx]))
            {
               // evaluating fitness of the circuit and save the result to avoid re-evaluation
               bool is_valid = true;
               fitness_parents[i] = evaluate_circuit(size_per_cir, &parents[starting_idx], sim_params, is_valid);
               if (!is_valid)
                  continue;
               break;
            }
         }
         // std::cout << "finish generate parent " << i << std::endl;
      }
   }
   // update current best result outside the parallel region to avoid synchronization
   for (int i = 0; i < parameters.parent_pool_size; i++)
   {
      if (fitness_parents[i] > best_fitness || nottouched)
      {
         nottouched = false;
         best_fitness = fitness_parents[i];
         generalized_copy(0, i, best_circuit, parents, 1);
         old_best_fitness = best_fitness;
      }
   }
}

/**
 * @brief Destructor of the GA class.
 * This destructor releases dynamically allocated memory.
 */
GA::~GA()
{
   delete[] parents;
   delete[] best_circuit;
   delete[] population;
   delete[] fitness_parents;
   delete[] fitness_population;
   switch (parameters.selection_scheme)
   {
   case 1:
      delete[] rank_sum;
      delete[] rank_val;
      delete[] population_idex_list;
      break;
   case 2:
      delete[] roulette_sum;
      break;

   default:
      break;
   }

#ifdef OUT_RESULT_TO_FILE
   delete[] file_write_buffer;
#endif

#ifdef PARALLEL_ENABLE
   if (parameters.parallel_mpi == 1)
   {
      if (id == 0)
      {
         delete[] reduction_buffer_fitness;
         delete[] all_parents;
      }
      MPI_Finalize();
   }
#endif
}

/**
 * @brief This function checks if two arrays are the same.
 *
 * @param p1 Pointer to the first array.
 * @param p2 Pointer to the second array.
 * @return true if the two arrays are the same, false otherwise.
 */
bool GA::check_same(int *p1, int *p2)
{
   for (int i = 0; i < size_per_cir; i++)
   {
      if (p1[i] != p2[i])
      {
         return false;
      }
   }
   return true;
}

/**
 * @brief This function generates all the children for next iteration in the genetic algorithm.
 */
void GA::generate_child()
{
   int p1, p2;
   double mutation_rand;
   double crossover_rand;
   int ct1 = 0;
   int ct2 = 0;
   int ct3 = 0;
   // parallelizing the generation
#pragma omp parallel
   {
      CCircuit checker(size_per_cir);
#pragma omp for reduction(+ : ct1, ct2, ct3)
      for (int i = 0; i < parameters.population_size; i++)
      {
         int starting_idx = get_idx(i, 0);
         int cur_idx;

         int flag = 0;
         int local_iter_count_1 = 0;
         while (local_iter_count_1 < parameters.max_iter_before_fail)
         {
            local_iter_count_1++;
            p1 = rand() % parameters.parent_pool_size;
            // determin if crossover should be performed
            crossover_rand = rand() / double(RAND_MAX);
            if (crossover_rand <= parameters.crossover_rate)
            {

               flag = 1;
               p2 = rand() % parameters.parent_pool_size;
               while (p2 == p1)
               {
                  p2 = rand() % parameters.parent_pool_size;
               }
               // std::cout<<"mating "<<p1<<","<<p2<<std::endl;
               two_point_crossover(p1, p2, i);
            }
            else
            {
               flag = 0;
               copy_parent_to_child(p1, i);
            }
            // determin if and what kinds of mutation should be performed.
            switch (parameters.mutation_scheme)
            {
            // mutating individual gene independently
            case 1:
            {
               int starting_idx = get_idx(i, 0);

               for (int j = 1; j < size_per_cir; j++)
               {
                  mutation_rand = rand() / double(RAND_MAX);
                  if (mutation_rand <= parameters.mutation_rate)
                  {
                     population[starting_idx + j] = rand() % (parameters.circuit_size + 2);
                  }
               }
               break;
            }
            // determin the chromosome as a whole will have a mutation and where the mutation will be
            case 2:
            {
               mutation_rand = rand() / double(RAND_MAX);
               if (mutation_rand <= parameters.mutation_rate)
               {
                  for (int j = 0; j < parameters.mutation_size; j++)
                  {
                     mutation(i);
                  }
               }
               break;
            }

            default:
               break;
            }
            // check validity
            if (checker.check_validity(&population[starting_idx]))
            {
               // evaluating circuit
               bool is_valid = true;
               fitness_population[i] = evaluate_circuit(size_per_cir, &population[starting_idx], sim_params, is_valid);
               if (!is_valid)
                  continue;
#ifdef PERFORMANCE_ANALYSIS_MACRO
               if (flag)
               {
                  ct1++;
                  int start1 = get_idx(p1, 0);
                  int start2 = get_idx(p2, 0);
                  if (check_same(&parents[start1], best_circuit) || check_same(&parents[start2], best_circuit))
                     ct3++;
               }
               else
               {
                  ct2++;
               }
#endif

               break;
            }
         }
      }
   }
   // updating best result
   for (int i = 0; i < parameters.population_size; i++)
   {
      if (fitness_population[i] > best_fitness || nottouched)
      {
         best_fitness = fitness_population[i];
         generalized_copy(0, i, best_circuit, population, 1);
         nottouched = false;
      }
   }
   // write to file if the macro is defined
#ifdef OUT_RESULT_TO_FILE
   if (cur_iter % parameters.write_interval == 0)
      write_to_file();
#endif
   mate_count_all += ct1;
   same_mate += ct3;
}

/**
 * @brief This function writes the current best solution of the genetic algorithm to a file.
 */
void GA::write_to_file()
{
   int dummy_int;
   evaluate_circuit_write(size_per_cir, best_circuit, file_write_buffer, sim_params);
   // filename of format outputfile_10.dat where 10 is the current iteration count
   // outputfile_1_10.dat if MPI is used and 1 is the process id
   std::string filename = "out/outputfile_";
#ifdef PARALLEL_ENABLE
   filename += std::to_string(id)+"_";
#endif
   filename += std::to_string(cur_iter) + ".dat";
   std::stringstream fname;
   std::fstream f1;
   std::cout << filename << std::endl;
   fname << filename;
   f1.open(fname.str().c_str(), std::ios_base::out);
   for (int i = 0; i < size_per_cir; i++)
   {
      f1 << best_circuit[i] << " ";
   }
   f1 << std::endl;
   for (int i = 0; i < size_per_cir - 1; i++)
   {
      f1 << file_write_buffer[i] << " ";
   }
   f1.close();
}

/**
 * @brief Main function to optimize the genetic algorithm
 *
 * Depending on the defined macros, the function decides to execute parallel optimization
 * or non-parallel optimization.
 */
void GA::optimize()
{
#ifdef PARALLEL_ENABLE
   if (parameters.parallel_mpi == 1)
   {
      optimize_parallel();
   }
   else
   {
      optimize_without_mpi();
   }
#endif
#ifndef PARALLEL_ENABLE
   optimize_without_mpi();
#endif
}

/**
 * @brief Executes the genetic algorithm optimization without MPI parallelism
 *
 * The function iterates through a series of generations, each time generating a new population of possible solutions.
 * The solutions are selected and modified to form the next generation.
 */
void GA::optimize_without_mpi()
{
   int iter = 0;
   double diff = parameters.tol + 1;
   not_improving_count = 0;
   cur_iter = 0;
   std::cout << "Starting GA:" << std::endl;
   // iteration is stopped when the best result is not improving in certain amount of iterations
   while (iter < parameters.max_iterations && not_improving_count <= parameters.max_iter_without_progress)
   {
      if (iter % 20 == 0)
      std::cout << "Iteration: " << iter << std::endl;
#ifdef PERFORMANCE_ANALYSIS_MACRO
      std::cout << best_fitness << std::endl;
      std::cout << "mating " << mate_count_all << " times,mate with best vector" << same_mate << " times " << (double)same_mate / mate_count_all << std::endl;
#endif
      generate_child();
      selection();
      if (old_best_fitness >= best_fitness)
      {
         not_improving_count++;
         parameters.mutation_rate *= parameters.mutation_rate_increase_factor;
      }
      else
      {
         old_best_fitness = best_fitness;
         not_improving_count = 0;
      }
      iter++;
      cur_iter = iter;
   }
   std::cout << "Converged in " << iter << " iterations." << std::endl;
   std::cout << "Best circuit vector: " << std::endl;
   for (int i = 0; i < size_per_cir; i++)
   {
      std::cout << best_circuit[i] << " ";
   }
   std::cout << std::endl;
   std::cout << "Fitness score: " << best_fitness << std::endl;
#ifdef PERFORMANCE_ANALYSIS_MACRO
   std::cout << "mating " << mate_count_all << " times,mate with best vector" << same_mate << " times " << std::endl;
#endif

#ifdef OUT_RESULT_TO_FILE
   write_to_file();
#endif
}

/**
 * @brief Performs mutation operation on the population at a given index
 *
 * @param population_idx Index of the population member to be mutated
 *
 * The function randomly selects an index within the member of the population,
 * then replaces the value at that index with a new random value.
 */
void GA::mutation(int population_idx)
{
   // mutation index number
   int idx = rand() % (size_per_cir - 1) + 1;
   // mutation node number
   int node_num = get_node_num(idx);
   // first element of mutation array in 1d structure
   int starting_idx = get_idx(population_idx, 0);
   int valid_val = rand() % (parameters.circuit_size + 2);
   while (valid_val == node_num || (valid_val >= parameters.circuit_size && node_num == population[starting_idx]))
   {

      valid_val = rand() % (parameters.circuit_size + 2);
   }
   population[starting_idx + idx] = valid_val;
}

/**
 * @brief Sets the best solution to the parent pool
 *
 * This is typically done after a new generation is created,
 * and the best solution from the previous generation is kept for comparison.
 */
void GA::set_best_to_parent()
{
   fitness_parents[0] = best_fitness;
   generalized_copy(0, 0, parents, best_circuit, 1);
}

/**
 * @brief Executes tournament selection
 *
 * Selects the best solutions from a randomly chosen subset of the population.
 * The chosen solutions are then used to form the parent pool for the next generation.
 */
void GA::tournament_selection()
{
   int cur_idx;
   set_best_to_parent();
   for (int i = 1; i < parameters.parent_pool_size; i++)
   {
      double best = -1;
      int best_idx = -1;
      for (int k = 0; k < parameters.tournament_size; k++)
      {
         cur_idx = rand() % parameters.population_size;
         if (fitness_population[cur_idx] > best)
         {
            best = fitness_population[cur_idx];
            best_idx = cur_idx;
         }
      }
      generalized_copy(i, best_idx, parents, population, 1);
      fitness_parents[i] = fitness_population[best_idx];
   }
}

/**
 * @brief Performs the selection process on the current population
 *
 * Depending on the specified selection scheme (rank, roulette or tournament),
 * it selects the most suitable individuals from the population to be parents for the next generation.
 */
void GA::selection()
{
   switch (parameters.selection_scheme)
   {
   case 1:
      rank_selection();
      break;
   case 2:
      roulette_selection();
      break;
   case 3:
      tournament_selection();
   default:
      break;
   }
}

/**
 * @brief Perform a binary search
 * @param list pointer of array of interest
 * @param target the value wanted
 * @param size size of the array
 * 
 */
int binary_search(double *list, double target, int size)
{
   int idx_l = 0;
   int idx_r = size - 1;
   int idx_mid = size / 2;
   while (true)
   {
      if (target <= list[idx_mid])
      {
         if (idx_mid == 0 || target > list[idx_mid - 1])
         {
            return idx_mid;
         }
         idx_r = idx_mid;
         idx_mid = (idx_l + idx_r) / 2;
      }
      else
      {
         if (idx_mid >= size - 1)
         {
            return idx_mid;
         }
         if (target <= list[idx_mid + 1])
         {
            return idx_mid + 1;
         }
         idx_l = idx_mid;
         idx_mid = (idx_l + idx_r) / 2;
      }
   }
}
/**
 * @brief Performs the roulette selection process on the current population
 */
void GA::roulette_selection()
{
   double min_val = 0;
   int flag = 1;
   for (int i = 0; i < parameters.population_size; i++)
   {
      if (min_val > fitness_population[i] || flag == 1)
      {
         min_val = fitness_population[i];
         flag = 0;
      }
   }
   // make the roulette selection works for negative fitness
   roulette_sum[0] = fitness_population[0] - min_val + 1;
   for (int i = 1; i < parameters.population_size; i++)
   {
      roulette_sum[i] = roulette_sum[i - 1] + fitness_population[i] - min_val + 1;
   }
   double p_rand;
   int res;
   set_best_to_parent();
   for (int i = 1; i < parameters.parent_pool_size; i++)
   {
      p_rand = rand() / (double)RAND_MAX;
      res = binary_search(roulette_sum, p_rand * roulette_sum[parameters.population_size - 1], parameters.population_size);
      generalized_copy(i, res, parents, population, 1);
      fitness_parents[i] = fitness_population[res];
   }
}
/**
 * @brief calculating a value based on ranking. in general, higher ranking will result in higher value
 * @param ranking ranking
 * @return the corresponding value of this ranking
 */
double GA::rank_func(int ranking)
{
   return parameters.population_size - ranking;
}
/**
 * @brief Performs the rank selection process on the current population
 *
 */
void GA::rank_selection()
{
   reinitialize_pil();
   in_place_sort(population_idex_list, parameters.population_size, fitness_population);
   for (int i = 0; i < parameters.population_size; i++)
   {
      rank_val[i] = rank_func(parameters.population_size - i);
   }
   rank_sum[0] = rank_val[0];
   for (int i = 1; i < parameters.population_size; i++)
   {
      rank_sum[i] = rank_sum[i - 1] + rank_val[i];
   }
   double p_rand;
   int res;
   set_best_to_parent();
   for (int i = 1; i < parameters.parent_pool_size; i++)
   {
      p_rand = rand() / (double)RAND_MAX;
      res = binary_search(rank_sum, p_rand * rank_sum[parameters.population_size - 1], parameters.population_size);
      res = population_idex_list[res];
      generalized_copy(i, res, parents, population, 1);
      fitness_parents[i] = fitness_population[res];
   }
}
/**
 * @brief Performs the two point crossover
 * @param idx1 index of the first parent
 * @param idx2 index of the second parent
 * @param dest index of the population which will be overwritten
 */
void GA::two_point_crossover(int idx1, int idx2, int dest)
{
   int idx_p1 = get_idx(idx1, 0);
   int idx_p2 = get_idx(idx2, 0);
   int *parent1 = &parents[idx_p1];
   int *parent2 = &parents[idx_p2];
   int *temp_p = NULL;
   int pos1 = rand() % (size_per_cir);
   int pos2 = rand() % (size_per_cir);
   int pop_starting = get_idx(dest, 0);
   while (pos1 == pos2)
   {
      pos2 = rand() % (size_per_cir);
   }

   if (pos1 > pos2)
   {
      std::swap(pos1, pos2);
   }
   for (int i = 0; i < size_per_cir; i++)
   {
      if (i <= pos2 && i >= pos1)
      {
         temp_p = parent1;
      }
      else
      {
         temp_p = parent2;
      }
      population[pop_starting + i] = temp_p[i];
   }
}
/**
 * @brief Helper function for sorting
 * @param arr pointer to the array to sort
 * @param tailing_val pointer to the array that should sort according to arr
 * @param start starting index
 * @param end ending index
 * @return pivoting index
 */
int partition(double *arr, int *tailing_val, int start, int end)
{

   double pivot = arr[start];

   int count = 0;
   for (int i = start + 1; i <= end; i++)
   {
      if (arr[i] <= pivot)
         count++;
   }

   // Giving pivot element its correct position
   int pivotIndex = start + count;
   std::swap(arr[pivotIndex], arr[start]);
   std::swap(tailing_val[pivotIndex], tailing_val[start]);

   // Sorting left and right parts of the pivot element
   int i = start, j = end;

   while (i < pivotIndex && j > pivotIndex)
   {

      while (arr[i] <= pivot)
      {
         i++;
      }

      while (arr[j] > pivot)
      {
         j--;
      }

      if (i < pivotIndex && j > pivotIndex)
      {
         std::swap(arr[i], arr[j]);
         std::swap(tailing_val[i], tailing_val[j]);
         i++;
         j--;
      }
   }

   return pivotIndex;
}
/**
 * @brief Custom in-place sort, tailing_val will be sort according to arr.
 * @param arr pointer to the array to sort
 * @param tailing_val pointer to the array that should sort according to arr
 * @param start starting index
 * @param end ending index
 */
void quickSort(double *arr, int *tailing_val, int start, int end)
{

   // base case
   if (start >= end)
      return;

   // partitioning the array
   int p = partition(arr, tailing_val, start, end);

   // Sorting the left part
   quickSort(arr, tailing_val, start, p - 1);

   // Sorting the right part
   quickSort(arr, tailing_val, p + 1, end);
}
/**
 * @brief in-place sorting wrapper
 * @param tailing_val pointer to the array that should sort according to arr
 * @param size size of array
 * @param tosort pointer to the array to sort
 */
void GA::in_place_sort(int *tailing_val, int size, double *tosort)
{
   quickSort(tosort, tailing_val, 0, size - 1);
}
/**
 * @brief copy from parent to population
 * @param parentid index of parent to be copied
 * @param childid index of population to be copied to
 */
void GA::copy_parent_to_child(int parentid, int childid)
{
   int starting_p_idx = get_idx(parentid, 0);
   int starting_c_idx = get_idx(childid, 0);
   for (int i = 0; i < size_per_cir; i++)
   {
      population[starting_c_idx + i] = parents[starting_p_idx + i];
   }
}
// copy id2 th item to id2 + count -1 th item of arr2 to id1 th item to arr1
//  id1 , arr1 : dest
//  id2, arr2: org
/**
 * @brief generalized copy of two array of same structure
 * @param id1 idex in arr1 to copy to 
 * @param id2 idex in arr2 to copy
 * @param arr1 destination array
 * @param arr2 source array
 * @param count number of items to copy
 */
void GA::generalized_copy(int id1, int id2, int *arr1, int *arr2, int count)
{
   int starting_idx1 = get_idx(id1, 0);
   int starting_idx2 = get_idx(id2, 0);
   for (int i = 0; i < size_per_cir * count; i++)
   {
      arr1[starting_idx1 + i] = arr2[starting_idx2 + i];
   }
}
/**
 * @brief get the node number 
 * @param pos index in 1d array
 * @return the node number 
 */
int GA::get_node_num(int pos)
{
   return (pos - 1) / 2;
}
/**
 * @brief reset the population_idex_list
 */
void GA::reinitialize_pil()
{
   for (int i = 0; i < parameters.population_size; i++)
   {
      population_idex_list[i] = i;
   }
}
/**
 * @brief update best result
 */
void GA::update_best()
{
   for (int i = 0; i < parameters.parent_pool_size; i++)
   {
      if (fitness_parents[i] > best_fitness)
      {
         best_fitness = fitness_parents[i];
         generalized_copy(0, i, best_circuit, parents, 1);
      }
   }
}
#ifdef PARALLEL_ENABLE
/**
 * @brief Sets up the parallel computation environment
 *
 * If parallel computation is enabled, this function sets up the necessary data structures
 * and communication parameters to execute the genetic algorithm across multiple processors.
 */
void GA::setup_parallel()
{

   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   MPI_Comm_size(MPI_COMM_WORLD, &p_count);
   if (id == 0)
   {
      reduction_buffer_fitness = new double[parameters.parent_comm_size_parallel * p_count];
      all_parents = new int[parameters.parent_comm_size_parallel * p_count * size_per_cir];
   }
   else
   {
      reduction_buffer_fitness = nullptr;
      all_parents = nullptr;
      // reduction_buffer_index = NULL;
      // reduction_buffer_pid = NULL;
   }
}

/**
 * @brief Executes the genetic algorithm optimization in parallel
 *
 * This function executes the optimization process across multiple processors in parallel,
 * using MPI to handle inter-processor communication.
 * The function iterates through a series of generations, with each processor maintaining its own population.
 */
void GA::optimize_parallel()
{

   int iter = 0;
   double diff = parameters.tol + 1;
   cur_iter = 0;
   if(id==0)
   {
      std::cout << "GA starting: " << std::endl;

   }
   while (iter < parameters.max_iterations)
   {
      generate_child();
      selection();
      if (iter % parameters.communicate_interval == 0)
      {
         mpi_reduction_wrapper();
      }
      iter++;
      cur_iter=iter;
      parameters.mutation_rate *= parameters.mutation_rate_increase_factor;
   }

   comm_loc.fitness = best_fitness;
   comm_loc.pid = id;
   MPI_Allreduce(MPI_IN_PLACE, &comm_loc, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
   MPI_Bcast(best_circuit, size_per_cir, MPI_DOUBLE, comm_loc.pid, MPI_COMM_WORLD);
   if (id == 0)
   {
      std::cout << "Best circuit vector: " << std::endl;
      for (int i = 0; i < size_per_cir; i++)
      {
         std::cout << best_circuit[i] << " ";
      }
      std::cout << std::endl;
      std::cout << "Fitness score: " << best_fitness << std::endl;
   }
   #ifdef OUT_RESULT_TO_FILE
   write_to_file();
   #endif
   // return 0;
}

/**
 * @brief Wrapper function for MPI reduction operation
 *
 * This function organizes the gathered data from each processor, then carries out an MPI reduction operation
 * to combine the data in a meaningful way, typically to find the best solution amongst all processors.
 */
void GA::mpi_reduction_wrapper()
{
   MPI_Gather(fitness_parents, parameters.parent_pool_size, MPI_DOUBLE,
              reduction_buffer_fitness, parameters.parent_pool_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Gather(parents, parameters.parent_pool_size * size_per_cir, MPI_INT,
              all_parents, parameters.parent_comm_size_parallel * size_per_cir, MPI_INT, 0, MPI_COMM_WORLD);
   if (id == 0)
   {
      // for (int i = 0; i < parameters.parent_pool_size * parameters.parallel; i++)
      // {
      //    reduction_buffer_index[i] = i % parameters.parallel;
      //    reduction_buffer_pid[i] = i / parameters.parent_pool_size;
      // }
      tournament_selection_parallel();
   }
   MPI_Bcast(fitness_parents, parameters.parent_comm_size_parallel, MPI_DOUBLE,
             0, MPI_COMM_WORLD);
   MPI_Bcast(parents, parameters.parent_comm_size_parallel * size_per_cir, MPI_INT,
             0, MPI_COMM_WORLD);
   update_best();
}
// void GA::roulette_selection_parallel()
// {
//    roulette_sum[0] = fitness_population[0];
//    for (int i = 1; i < parameters.population_size; i++)
//    {
//       roulette_sum[i] = roulette_sum[i - 1] + fitness_population[i];
//    }
//    double p_rand;
//    int res;
//    // set_best_to_parent();
//    for (int i = 1; i < parameters.parent_pool_size; i++)
//    {
//       p_rand = rand() / (double)RAND_MAX;
//       res = binary_search(roulette_sum, p_rand * roulette_sum[parameters.population_size - 1], parameters.population_size);
//       for (int j = 0; j < parameters.circuit_size; j++)
//       {
//          parents[i]->units[j].conc_num = population[res]->units[j].conc_num;
//          parents[i]->units[j].tails_num = population[res]->units[j].tails_num;
//       }
//       fitness_parents[i] = fitness_population[res];
//    }
// }
/**
 * @brief perform tournament selection on the gathered information
 */
void GA::tournament_selection_parallel()
{

   int cur_idx;
   // set_best_to_parent();
   for (int i = 0; i < parameters.parent_pool_size; i++)
   {
      double best = -1;
      int best_idx = -1;
      for (int k = 0; k < parameters.tournament_size_parallel; k++)
      {
         cur_idx = rand() % (parameters.parent_comm_size_parallel * p_count);
         if (reduction_buffer_fitness[cur_idx] > best)
         {
            best = reduction_buffer_fitness[cur_idx];
            best_idx = cur_idx;
         }
      }
      generalized_copy(i, best_idx, parents, all_parents, 1);
      fitness_parents[i] = reduction_buffer_fitness[best_idx];
   }
}

#endif
