/** Header for the unit class
 *
 *
 */

#pragma once
#include <vector>

// Class representing a single unit in the circuit
class CUnit
{
public:
  // Constructor should be based on the circuit vector
  CUnit(){}; // Default constructor
  CUnit(int unit_id, bool is_feed = false) : unit_id(unit_id), is_feed(is_feed){};
  ~CUnit(){};

  std::vector<int> concentrate_list;
  std::vector<int> tails_list;
  std::vector<int> concentrate_feed_list;
  std::vector<int> tails_feed_list;

  // List of units that this unit can feed into/accept feed from
  // This is instantiated based on the circuit vector

  double ger_recovery, waste_recovery;       // R in the paper
  double ger_flow_in, waste_flow_in;         // F_i, current
  double ger_flow_in_old, waste_flow_in_old; // F_i, old

  double c_ger_flow_out, c_waste_flow_out; // C_i, concentrate (contains both ger and waste)
  double t_ger_flow_out, t_waste_flow_out; // T_i, tail

  // Accessors
  bool get_is_feed() const { return is_feed; }
  int get_unit_id() const { return unit_id; }

private:
  // Unit ID
  int unit_id;
  // Boolean to denote the '0' node in the circuit
  bool is_feed;
};

