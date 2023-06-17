/** header file for the circuit simulator
 *
 * This header file defines the function that will be used to evaluate the circuit
 */

#include <vector>
#include "CUnit.h"

#pragma once

// Define a struct to hold simulation parameters
struct SimulationParameters
{
    double F0_ger;
    double F0_waste;
    double k_ger;
    double k_waste;
    double volume;
    double phi;
    double rho;

    double ger_reward;
    double waste_reward;

    double min_flowrate;
    int max_iter;
};

class CSimulator
{
public:
    // Constructors
    CSimulator(){};
    CSimulator(int vector_size, int *circuit_vector, const SimulationParameters &params);

    // Public functions
    bool solve();
    double calculate_reward(bool is_converge);

    // Accessors
    double get_total_concentrate_ger() { return total_concentrate_ger; };
    double get_total_concentrate_waste() { return total_concentrate_waste; };
    double get_total_tail_ger() { return total_tail_ger; };
    double get_total_tail_waste() { return total_tail_waste; };
    std::vector<int> get_final_concentrate_list() const;
    std::vector<int> get_final_tail_list() const;
    std::vector<CUnit> get_units_vector() { return units_vector; }
    void write_vector_values_to_array(double* dest);
    
private:
    void instantiate_circuit_units();
    void update_feed_unit();
    void update_from_guess();
    void swap_feed_flows();
    void add_feed_flows();
    void calculate_total_outputs();

    // Member variables
    std::vector<CUnit> units_vector;
    int vector_size;
    int *circuit_vector;
    SimulationParameters params;

    double total_concentrate_ger, total_concentrate_waste;
    double total_tail_ger, total_tail_waste;
    std::vector<int> final_concentrate_list;
    std::vector<int> final_tail_list;
};

double evaluate_circuit(int vector_size, int *circuit_vector, const SimulationParameters &params, bool & is_valid);
double evaluate_circuit_write(int vector_size, int *circuit_vector,double * edge_vector, const SimulationParameters &params );

