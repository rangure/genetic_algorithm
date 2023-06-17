#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>

#include "CUnit.h"
#include "CCircuit.h"
#include "CSimulator.h"
#include "Genetic_Algorithm.h"
#ifdef PARALLEL_ENABLE
#include <mpi.h>
#endif
using namespace std;

/** \brief Test GA performance on a circuit size of 5
 */
void test_circuit5()
{
    double rewards_ratio = 0.2;
    double flow_ratio = 0.1;

    Algorithm_Parameters AP;
    SimulationParameters SP;

    SP.F0_waste = 100;                    // Input flow of waste
    SP.F0_ger = SP.F0_waste * flow_ratio; // Input flow of gerardium

    SP.k_ger = 5.e-3;   // Rate constant for gerardium
    SP.k_waste = 5.e-4; // Rate constant for waste

    SP.volume = 10; // Volume of each cell
    SP.phi = 0.1;
    SP.rho = 3.e3; // Density

    SP.waste_reward = -500.0;                         // Negative reward for waste
    SP.ger_reward = -rewards_ratio * SP.waste_reward; // Reward for gerardium

    SP.min_flowrate = 1.e-10; // Minimum flowrate for the simulator
    SP.max_iter = 1000;

    AP.max_iterations = 1000;
    AP.tol = 0.001;
    AP.mutation_size = 3;
    AP.parent_pool_size = 10;
    AP.population_size = 150;
    AP.circuit_size = 5;
    AP.selection_scheme = 2;
    AP.mutation_rate = 0.04;
    AP.mutation_scheme = 1;
    AP.tournament_size = 5;
    AP.communicate_interval = 10;
    AP.crossover_rate = 0.7;
    AP.max_iter_without_progress = 100;
    AP.max_iter_before_fail = 10000;
    AP.mutation_rate_increase_factor = 1.001;

    AP.parallel_mpi = 0;
    AP.tournament_size_parallel = 4;
    AP.parent_comm_size_parallel = AP.parent_pool_size;
    AP.write_interval = 10000;
    GA ga = GA(AP, SP);
    ga.optimize();
    assert(ga.best_fitness>100);
}

 /** \brief Test GA performance on a circuit size of 10
 */
void test_circuit10()
{
    double rewards_ratio = 0.2;
    double flow_ratio = 0.1;

    Algorithm_Parameters AP;
    SimulationParameters SP;

    SP.F0_waste = 100;                    // Input flow of waste
    SP.F0_ger = SP.F0_waste * flow_ratio; // Input flow of gerardium

    SP.k_ger = 5.e-3;   // Rate constant for gerardium
    SP.k_waste = 5.e-4; // Rate constant for waste

    SP.volume = 10; // Volume of each cell
    SP.phi = 0.1;
    SP.rho = 3.e3; // Density

    SP.waste_reward = -500.0;                         // Negative reward for waste
    SP.ger_reward = -rewards_ratio * SP.waste_reward; // Reward for gerardium

    SP.min_flowrate = 1.e-10; // Minimum flowrate for the simulator
    SP.max_iter = 1000;

    AP.max_iterations = 1000;
    AP.tol = 0.001;
    AP.mutation_size = 3;
    AP.parent_pool_size = 10;
    AP.population_size = 150;
    AP.circuit_size = 10;
    AP.selection_scheme = 2;
    AP.mutation_rate = 0.04;
    AP.mutation_scheme = 1;
    AP.tournament_size = 10;
    AP.communicate_interval = 10;
    AP.crossover_rate = 0.7;
    AP.max_iter_without_progress = 100;
    AP.max_iter_before_fail = 10000;
    AP.mutation_rate_increase_factor = 1.001;

    AP.parallel_mpi = 0;
    AP.tournament_size_parallel = 4;
    AP.parent_comm_size_parallel = AP.parent_pool_size;
    AP.write_interval = 10000;
    GA ga = GA(AP, SP);
    ga.optimize();
    assert(ga.best_fitness>350);
}

/** \brief Test GA performance on a circuit size of 20
 */
void test_circuit20()
{
    double rewards_ratio = 0.2;
    double flow_ratio = 0.1;

    Algorithm_Parameters AP;
    SimulationParameters SP;

    SP.F0_waste = 100;                    // Input flow of waste
    SP.F0_ger = SP.F0_waste * flow_ratio; // Input flow of gerardium

    SP.k_ger = 5.e-3;   // Rate constant for gerardium
    SP.k_waste = 5.e-4; // Rate constant for waste

    SP.volume = 10; // Volume of each cell
    SP.phi = 0.1;
    SP.rho = 3.e3; // Density

    SP.waste_reward = -500.0;                         // Negative reward for waste
    SP.ger_reward = -rewards_ratio * SP.waste_reward; // Reward for gerardium

    SP.min_flowrate = 1.e-10; // Minimum flowrate for the simulator
    SP.max_iter = 1000;

    AP.max_iterations = 1000;
    AP.tol = 0.001;
    AP.mutation_size = 3;
    AP.parent_pool_size = 10;
    AP.population_size = 150;
    AP.circuit_size = 20;
    AP.selection_scheme = 2;
    AP.mutation_rate = 0.04;
    AP.mutation_scheme = 1;
    AP.tournament_size = 5;
    AP.communicate_interval = 10;
    AP.crossover_rate = 0.7;
    AP.max_iter_without_progress = 100;
    AP.max_iter_before_fail = 10000;
    AP.mutation_rate_increase_factor = 1.001;

    AP.parallel_mpi = 0;
    AP.tournament_size_parallel = 4;
    AP.parent_comm_size_parallel = AP.parent_pool_size;
    AP.write_interval = 10000;
    GA ga = GA(AP, SP);
    ga.optimize();
    assert(ga.best_fitness>600);
}
/** \brief Main function to run the tests
 *
 *  This main function runs all the defined test functions.
 */
// Temporary main
int main()
{
    test_circuit5();
    test_circuit10();
    test_circuit20();
    return 0;
}

// //------------------------------------------------------------------------------
