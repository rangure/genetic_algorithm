/**
 * @file
 * @brief Main code file for the CSimulator class.
 */

#include "CUnit.h"
#include "CCircuit.h"
#include "CSimulator.h"

// #define DEBUG

#include <vector>
#include <iostream>

#include <cmath>

/**
 * @brief CSimulator constructor.
 * @param vector_size Size of the circuit_vector.
 * @param circuit_vector Vector representing the circuit.
 * @param params Parameters for the simulation.
 */
CSimulator::CSimulator(int vector_size, int *circuit_vector, const SimulationParameters &params) : vector_size(vector_size), circuit_vector(circuit_vector), params(params)
{
  instantiate_circuit_units();
};

//------------------------------------------------------------------------------
/**
 * @brief Update the feed unit of the simulator.
 */
void CSimulator::update_feed_unit()
{
  units_vector[circuit_vector[0]].ger_flow_in = params.F0_ger;
  units_vector[circuit_vector[0]].waste_flow_in = params.F0_waste;
}
//------------------------------------------------------------------------------
/**
 * @brief Instantiate the circuit units.
 */
void CSimulator::instantiate_circuit_units()
{
  int n_units = vector_size / 2; // Rounds down, e.g. 5/2 = 2
  // std::vector<CUnit> units_vector;
  units_vector.resize(n_units);

  if (n_units != 0)
  {
    // Instantiate the feed unit
        units_vector[circuit_vector[0]] = CUnit(circuit_vector[0], true);
    update_feed_unit();
  }

  for (int i = 0; i < n_units; i++)
  {
     if (i != circuit_vector[0])
     {
       units_vector[i] = CUnit(i, false);
     }
  }

  //
  for (int i = 0; i < n_units; i++)
  {
    int current_concentrate = circuit_vector[2 * i + 1];
    int current_tails = circuit_vector[2 * i + 2];

    // For the current unit, add the concentrate and tail to the appropriate lists
    units_vector[i].concentrate_list.push_back(current_concentrate); // Concentrate
    units_vector[i].tails_list.push_back(current_tails);             // Tail

    // For the concentrate unit, add the current unit to the feed list
    if (current_concentrate >= n_units)
    {
      // Final concentrate unit
      final_concentrate_list.push_back(i);
    }
    else
    {
      // Add the current unit to the feed list of the concentrate unit
      units_vector[current_concentrate].concentrate_feed_list.push_back(i);
    }

    if (current_tails >= n_units)
    {
      // Final tails unit
      final_tail_list.push_back(i);
    }
    else
    {
      units_vector[current_tails].tails_feed_list.push_back(i);
    }
  }
}
//------------------------------------------------------------------------------
/**
 * @brief Update units from guess.
 */
void CSimulator::update_from_guess()
{
  for (int i = 0; i < units_vector.size(); i++)
  {

    // Here, we don't calculate but rather take the guess as the input
    double F_ger_sum = std::max(units_vector[i].ger_flow_in, params.min_flowrate);
    double F_waste_sum = std::max(units_vector[i].waste_flow_in, params.min_flowrate);

    // Calculate the residence time (temporary variable used to calculate recovery)
    double tau = params.volume * params.phi * params.rho / (F_ger_sum + F_waste_sum);

    // Calculate the recovery
    units_vector[i].ger_recovery = params.k_ger * tau / (1 + params.k_ger * tau);
    units_vector[i].waste_recovery = params.k_waste * tau / (1 + params.k_waste * tau);

    // Calculate the concentrate flow out
    units_vector[i].c_ger_flow_out = units_vector[i].ger_recovery * F_ger_sum;
    units_vector[i].c_waste_flow_out = units_vector[i].waste_recovery * F_waste_sum;

    // Calculate the tail flow out
    units_vector[i].t_ger_flow_out = (1 - units_vector[i].ger_recovery) * F_ger_sum;
    units_vector[i].t_waste_flow_out = (1 - units_vector[i].waste_recovery) * F_waste_sum;
  }
}

//------------------------------------------------------------------------------
/**
 * @brief Add feed flows.
 */
void CSimulator::add_feed_flows()
{
  // Add the tail and concentrate flows to the feed flows
  for (int i = 0; i < units_vector.size(); i++)
  {
    for (int j = 0; j < units_vector[i].concentrate_feed_list.size(); j++)
    {
      // ID of the unit that feeds concentrate into this unit
      int feed_unit_id = units_vector[i].concentrate_feed_list[j];
      // Add the flows from the concentrate and tail to the feed flows in the current unit
      units_vector[i].ger_flow_in += units_vector[feed_unit_id].c_ger_flow_out;
      units_vector[i].waste_flow_in += units_vector[feed_unit_id].c_waste_flow_out;
    }
    for (int j = 0; j < units_vector[i].tails_feed_list.size(); j++)
    {
      // ID of the unit that feeds tails into this unit
      int feed_unit_id = units_vector[i].tails_feed_list[j];
      // Add the flows from the concentrate and tail to the feed flows in the current unit
      units_vector[i].ger_flow_in += units_vector[feed_unit_id].t_ger_flow_out;
      units_vector[i].waste_flow_in += units_vector[feed_unit_id].t_waste_flow_out;
    }
  }
}
//------------------------------------------------------------------------------
/**
 * @brief Swap feed flows.
 */
void CSimulator::swap_feed_flows()
// Assing new feeds to old feeds and zero the new feeds
{
  for (int i = 0; i < units_vector.size(); i++)
  {
    units_vector[i].ger_flow_in_old = units_vector[i].ger_flow_in;
    units_vector[i].ger_flow_in = 0.0;

    units_vector[i].waste_flow_in_old = units_vector[i].waste_flow_in;
    units_vector[i].waste_flow_in = 0.0;
  }
}
//------------------------------------------------------------------------------
// Calculate the total outputs
/**
 * @brief Calculate the total outputs.
 */
void CSimulator::calculate_total_outputs()
{
  total_concentrate_ger = 0.0;
  total_concentrate_waste = 0.0;
  total_tail_ger = 0.0;
  total_tail_waste = 0.0;

  // Go through all the final concentrate units
  for (int i = 0; i < final_concentrate_list.size(); i++)
  {
    int current_unit_id = final_concentrate_list[i];
    total_concentrate_ger += units_vector[current_unit_id].c_ger_flow_out;
    total_concentrate_waste += units_vector[current_unit_id].c_waste_flow_out;
  }
  for (int i = 0; i < final_tail_list.size(); i++)
  {
    int current_unit_id = final_tail_list[i];
    total_tail_ger += units_vector[current_unit_id].t_ger_flow_out;
    total_tail_waste += units_vector[current_unit_id].t_waste_flow_out;
  }
}
//------------------------------------------------------------------------------
/**
 * @brief Calculate the reward.
 * @param is_converge Boolean parameter determining if the simulator has converged.
 * @return Reward for the simulation.
 */
double CSimulator::calculate_reward(bool is_converge)
{
  if(is_converge)
    return params.ger_reward * total_concentrate_ger + params.waste_reward * total_concentrate_waste;
  else
    return -1;
}
//------------------------------------------------------------------------------
/**
 * @brief Solve the simulation.
 * @return True if the simulation was solved successfully, false otherwise.
 */
bool CSimulator::solve()
{
  // Make an initial guess for feed flow rates to every component
  for (int i = 0; i < units_vector.size(); i++)
  {
    units_vector[i].ger_flow_in = params.F0_ger / units_vector.size();
    units_vector[i].waste_flow_in = params.F0_waste / units_vector.size();
  }

  double max_residual = 1.0; // Must be bigger than 1.e-6 so that we enter the while loop
  // Iterate until convergence
  int it = -1;

  while (max_residual > 1.e-6 && it < 1000)
  {
    it++;
    update_from_guess(); // Step 2
    swap_feed_flows();   // Step 3
    update_feed_unit();  // Step 4
    add_feed_flows();    // Step 5

    // Step 6: Check for convergence
    max_residual = 0.0;

    for (int i = 0; i < units_vector.size(); i++)
    {
      double residual = std::abs(units_vector[i].ger_flow_in - units_vector[i].ger_flow_in_old) + std::abs(units_vector[i].waste_flow_in - units_vector[i].waste_flow_in_old);
#ifdef DEBUG
      std::cout << "Residual for units " << i << " is " << residual << std::endl;
#endif
      if (residual > max_residual)
      {
        max_residual = residual;
      }
    }
#ifdef DEBUG
    std::cout << "Max residual: " << max_residual << std::endl;
#endif
  }

#ifdef DEBUG
  std::cout << "Iteration: " << it << std::endl;
#endif
if (it == 1000)
  return false;

  calculate_total_outputs(); // Step 7

  // Print the results

#ifdef DEBUG
  for (int i = 0; i < units_vector.size(); i++)
  {
    std::cout << "Unit: " << i << std::endl;
    std::cout << "ger_flow_in: " << units_vector[i].ger_flow_in << std::endl;
    std::cout << "waste_flow_in: " << units_vector[i].waste_flow_in << std::endl;
    std::cout << "ger_flow_in_old: " << units_vector[i].ger_flow_in_old << std::endl;
    std::cout << "waste_flow_in_old: " << units_vector[i].waste_flow_in_old << std::endl;
    std::cout << "ger_recovery: " << units_vector[i].ger_recovery << std::endl;
    std::cout << "waste_recovery: " << units_vector[i].waste_recovery << std::endl;
    std::cout << "c_ger_flow_out: " << units_vector[i].c_ger_flow_out << std::endl;
    std::cout << "c_waste_flow_out: " << units_vector[i].c_waste_flow_out << std::endl;
    std::cout << "t_ger_flow_out: " << units_vector[i].t_ger_flow_out << std::endl;
    std::cout << "t_waste_flow_out: " << units_vector[i].t_waste_flow_out << std::endl;
  }

  std::cout << "Total concentrate ger: " << total_concentrate_ger << std::endl;
  std::cout << "Total concentrate waste: " << total_concentrate_waste << std::endl;
  std::cout << "Total tail ger: " << total_tail_ger << std::endl;
  std::cout << "Total tail waste: " << total_tail_waste << std::endl;
#endif

return true;
}
/**
 * @brief writing the values in CUnits to an array
 * @return destination write buffer
 */
void CSimulator::write_vector_values_to_array(double * dest)
{
  for(int i =0 ; i<units_vector.size();i++)
  {
    double flow_con, flow_tail;
    flow_con = units_vector[i].c_ger_flow_out+units_vector[i].c_waste_flow_out;
    flow_tail = units_vector[i].t_ger_flow_out+units_vector[i].t_waste_flow_out;
    dest[units_vector[i].get_unit_id()*2] = flow_con;
    dest[units_vector[i].get_unit_id()*2+1] = flow_tail;
  }
}

//------------------------------------------------------------------------------
// This function wraps around the simulator and returns the reward directly
/**
 * @brief Evaluate the circuit.
 * @param vector_size Size of the circuit_vector.
 * @param circuit_vector Vector representing the circuit.
 * @param params Parameters for the evaluation.
 * @return Reward for the evaluation.
 */
double evaluate_circuit(int vector_size, int *circuit_vector, const SimulationParameters &params, bool &is_valid)
{
  // Instantiate the simulator
  CSimulator simulator(vector_size, circuit_vector, params);

  // Solve the circuit
  is_valid = simulator.solve();
  // Return the reward
  return simulator.calculate_reward(is_valid);
}
//------------------------------------------------------------------------------
// This function wraps around the simulator and returns the value of each unit
/**
 * @brief Evaluate the circuit and save the circuit information in the given buffer 
 * @param vector_size Size of the circuit_vector.
 * @param edge_vector Save buffer
 * @param circuit_vector Vector representing the circuit.
 * @param params Parameters for the evaluation.
 * @return Reward for the evaluation.
 */
double evaluate_circuit_write(int vector_size, int *circuit_vector,double * edge_vector, const SimulationParameters &params )
{
  // Instantiate the simulator
  CSimulator simulator(vector_size, circuit_vector, params);
  bool is_valid;
  // Solve the circuit
  is_valid = simulator.solve();
  // Return the reward
  simulator.write_vector_values_to_array(edge_vector);
  return simulator.calculate_reward(is_valid);
}