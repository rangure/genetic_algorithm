#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>

#include "CUnit.h"
#include "CCircuit.h"
#include "CSimulator.h"

using namespace std;

/** \brief Test instantiation of circuit units
 *  
 *  This function sets up simulation parameters and creates a CSimulator instance 
 *  to verify correct instantiation and behavior of the circuit units.
 */
void test_instantiate_circuit_units() {

    // Set up the simulation parameters
    SimulationParameters sim_params;
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

    // Set up the circuit vector
    int circuit_vector[3] = {0, 1, 2};

    // Instantiate the simulator
    CSimulator simulator(3, circuit_vector, sim_params);

    // Check if the simulator was instantiated with correct number of units
    std::vector<CUnit> units_vector_test = simulator.get_units_vector();

    assert(units_vector_test.size() == 1);
}

/** \brief Check if two double values are close to each other within a certain epsilon.
 *
 *  \param a First double value
 *  \param b Second double value
 *  \param epsilon Maximum allowed difference between a and b
 *  \return True if |a - b| < epsilon, false otherwise
 */
bool isClose(double a, double b, double epsilon)
{
    return abs(a - b) < epsilon;
}

/** \brief Test the reward calculation of the simulator
 *  
 *  This function initializes simulation parameters and creates a CSimulator instance 
 *  to verify the correct calculation of the reward value.
 */
void test_calculate_reward() {
    // Initialize the simulation parameters
    SimulationParameters sim_params;
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

    // Set up the circuit vector
    int circuit_vector[11] = {4, 5, 1, 2, 4, 0, 1, 1, 6, 1, 3};

    // Instantiate the simulator
    CSimulator simulator(11, circuit_vector, sim_params);

    // Run the solve() function
    bool value = simulator.solve();
    double calc_reward = simulator.calculate_reward(value); //calculate the reward for that circuit
    double expected_reward = 107.204; //expected reward for that circuit

    assert(isClose(calc_reward, expected_reward, 1e-3)); //check if the calculated reward is close to the expected reward
}

/** \brief Test a small circuit simulation
 *  
 *  This function sets up simulation parameters for a small circuit, creates a CSimulator instance 
 *  and verifies the correct calculation of total concentrate and tail flows.
 */
void test_small_circuit()
{
    SimulationParameters sim_params;
    sim_params.F0_ger = 20;     // Input flow of gerardium
    sim_params.F0_waste = 80;  // Input flow of waste
    sim_params.k_ger = 5.e-3;   // Rate constant for gerardium
    sim_params.k_waste = 5.e-4; // Rate constant for waste
    sim_params.volume = 10;     // Volume of each cell
    sim_params.phi = 0.1;
    sim_params.rho = 3.e3;            // Density
    sim_params.ger_reward = 100.0;    // Reward for gerardium
    sim_params.waste_reward = -500.0; // Negative reward for waste
    sim_params.min_flowrate = 1.e-10; // Minimum flowrate for the simulator
    sim_params.max_iter = 1000;       // Maximum number of iterations for simulator

    // Set up the circuit vector
    int circuit_vector[3] = {0, 1, 2};

    // Instantiate the simulator
    CSimulator simulator(3, circuit_vector, sim_params);

    // Run the solve() function
    bool value = simulator.solve();

    // Check if the total concentrate and tail flows are correct
    double total_concentrate_ger = simulator.get_total_concentrate_ger();
    double total_concentrate_waste = simulator.get_total_concentrate_waste();
    double total_tail_ger = simulator.get_total_tail_ger();
    double total_tail_waste = simulator.get_total_tail_waste();

    // Check if the values are close to the expected values
    assert(isClose(total_concentrate_ger, 2.61, 1e-2));
    assert(isClose(total_concentrate_waste, 1.18, 1e-2));
    assert(isClose(total_tail_ger, 17.39, 1e-2));
    assert(isClose(total_tail_waste, 78.82, 1e-2));
}

/** \brief Main function to run the tests
 *  
 *  This main function runs all the defined test functions.
 */
//Temporary main
int main() {
    test_instantiate_circuit_units();
    test_calculate_reward();
    test_small_circuit();
    return 0;
}

// //------------------------------------------------------------------------------
