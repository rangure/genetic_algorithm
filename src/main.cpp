#include <iostream>

#include "CUnit.h"
#include "CCircuit.h"
#include "CSimulator.h"
#include "Genetic_Algorithm.h"
#ifdef PARALLEL_ENABLE
#include <mpi.h>
#endif

int main(int argc, char *argv[])
{

    double rewards_ratio = 0.2;
    double flow_ratio = 0.1;
    

    Algorithm_Parameters AP;
    SimulationParameters SP;

    SP.F0_waste = 100;  // Input flow of waste
    SP.F0_ger = SP.F0_waste * flow_ratio;     // Input flow of gerardium
    
    SP.k_ger = 5.e-3;   // Rate constant for gerardium
    SP.k_waste = 5.e-4; // Rate constant for waste
    
    SP.volume = 10;     // Volume of each cell
    SP.phi = 0.1;
    SP.rho = 3.e3;            // Density

    SP.waste_reward = -500.0; // Negative reward for waste
    SP.ger_reward = -rewards_ratio * SP.waste_reward;    // Reward for gerardium
    

    SP.min_flowrate = 1.e-10; // Minimum flowrate for the simulator
    SP.max_iter = 1000;

    AP.max_iterations = 1000;
    AP.tol = 0.001;
    AP.mutation_size = 3;
    AP.parent_pool_size = 10;
    AP.population_size = 150;
    AP.circuit_size = 15;
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
#ifdef PARALLEL_ENABLE
    if (AP.parallel_mpi)
        MPI_Init(&argc, &argv);
#endif
            GA ga = GA(AP,SP);
            ga.optimize();
}
