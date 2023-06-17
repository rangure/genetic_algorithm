/** Header for the Genetic Algorithm library
 *
 */
#pragma once
#include "CCircuit.h"
#include "CSimulator.h"
struct Algorithm_Parameters
    {
        // maximum iteration for GA
        // if MPI is used, each process will use this as their individual maximum iteration for GA
        int max_iterations;
        // tolerance for stopping, not used
        double tol;
        // individual mutation rate if mutation_scheme=1
        // chromosome mutation rate if mutation_scheme=2
        double mutation_rate;
        // number of points to mutate, only used when mutation_scheme=2
        int mutation_size;
        // parent pool size GA is keeping
        int parent_pool_size;
        // population size in each generation GA is keeping
        int population_size;
        // problem size
        int circuit_size;
        // selection scheme to use
        // rank selection: 1
        // roulette selection: 2
        // tournament selection: 3
        int selection_scheme;
        // tournament size, only used when selection_scheme=3
        int tournament_size;
        // crossover rate
        double crossover_rate;
        // if MPI should be enabled. only if parallel_mpi=1 and PARALLEL_ENABLE macro is defined the code will use MPI
        int parallel_mpi;
        // communication inverval, only used when MPI is used
        int communicate_interval;
        // tournament size, only used when MPI is used
        int tournament_size_parallel;
        // parent communication size, only used when MPI is used. currently this need to equal to parent_pool_size
        int parent_comm_size_parallel;
        // maximum continus iteration number without finding better result before stopping 
        int max_iter_without_progress;
        // mutation rate increase factor
        double mutation_rate_increase_factor;
        
        // maximum continus iteration number before break out of the loop
        // this is a fail safe mechanism for halting, should set to a large number
        int max_iter_before_fail;
        // iteration between output circuits to file, only used when OUT_RESULT_TO_FILE is defined
        int write_interval;
        // individual mutation if mutation_scheme=1
        // chromosome mutation if mutation_scheme=2
        int mutation_scheme;
    };

class GA
{
public:
    Algorithm_Parameters parameters;
    SimulationParameters sim_params;
    /** \brief Total size of parent-pool array*/
    int vec_size_parent;
    /** \brief Total size of population-pool array*/
    int vec_size_population;
    /** \brief Size per circuit array*/
    int size_per_cir;
    /** \brief Parameters for performance analysis*/
    int mate_count_all;
    /** \brief Parameters for performance analysis*/
    int same_mate;
    /** \brief Flag for updating the best fitness value*/
    bool nottouched;
    /** \brief Parent pool*/
    int *parents;
    /** \brief Population pool*/
    int *population;
    /** \brief Current best fitness*/
    double best_fitness;
    /** \brief Current best circuit*/
    int *best_circuit;
    /** \brief Fitness of corresponding population*/
    double *fitness_population;
    /** \brief Fitness of corresponding parent*/
    double *fitness_parents;
    /** \brief Rolling sum for roulette selection*/
    double *roulette_sum;
    /** \brief Rolling sum for rank selection*/
    double *rank_sum;
    /** \brief Ranking function values for each individual ranking*/
    double *rank_val;
    /** \brief Ranking of corresponding population*/
    int *population_idex_list;
    /** \brief Last best fitness*/
    double old_best_fitness;
    /** \brief Iteration count during which the best fitness is not improving*/
    int not_improving_count;
    /** \brief Current iteration number of Genetic algorithm*/
    int cur_iter;
    /** \brief Buffer to store the unit in/out flow information*/
    double * file_write_buffer;
    #ifdef PARALLEL_ENABLE
    /** \brief MPI process id*/
    int id;
    /** \brief MPI process count*/
    int p_count;
    /** \brief Buffer storing fitness for reduction*/
    double *reduction_buffer_fitness;
    /** \brief Buffer storing circuits for reduction*/
    int *all_parents;
    /** \brief Data structure for MPI_MAXLOC*/
    struct comm_loc
    {
        double fitness;
        int pid;
    } comm_loc;
    #endif
    GA();
    GA(Algorithm_Parameters AP, SimulationParameters SP);

    ~GA();
    void setup();
    void setup_parallel();
    void optimize();
    // int optimize(int, int *,
    //              double (&)(int, int *),
    //              struct Algorithm_Parameters);

    // Other overloaded functions for optimize

    void optimize_without_mpi();
    void two_point_crossover(int idx1, int idx2, int dest);
    void generate_child();
    void mutation(int population_idx);
    void selection();
    void rank_selection();
    void roulette_selection();
    void tournament_selection();
    double rank_func(int ranking);
    void set_best_to_parent();
    void in_place_sort(int *tailing_val, int size, double *tosort);
    void optimize_parallel();
    void mpi_reduction_wrapper();
    void tournament_selection_parallel();
    void roulette_selection_parallel();
    int get_idx(int pos_in_array, int pos_in_vec);
    int get_idx_node(int pos, int node, int val);
    int get_node_num(int pos);
    void copy_parent_to_child(int parentid, int childid);
    void generalized_copy(int id1, int id2, int*arr1, int*arr2, int count);
    void reinitialize_pil();
    void update_best();
    void select_mate(int& parent1, int& parent2);
    bool check_same(int *p1,int *p2);
    void write_to_file();
};

// Other functions and variables