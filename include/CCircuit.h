/** Header for the circuit class
 *
 * This header defines the circuit class and its associated functions
 *
*/

#pragma once

#include "CUnit.h"

#include <vector>

class CCircuit {
  public:
    CCircuit(int num_units);
    ~CCircuit();

    /**
     * @brief Checks the validity of the circuit.
     *
     * This function checks the validity of the circuit based on the following conditions:
     * 1. Every unit must be accessible from the feed.
     * 2. Every unit must have a route forward to both of the outlet streams.
     * 3. There should be no self-recycle.
     * 4. The destination for both products from a unit should not be the same unit.
     * 5. The destination id of each unit must be between 0 to num_units + 1.
     *    The first destination must not be tailings, the second destination must not be concentrate.
     *
     * @return True if the circuit is valid, false otherwise.
     */
    bool check_validity(int *circuit_vector);

  private:
    int num_units;
    int* circuit_vector;
    bool* marks;
    bool* have_conc;
    bool* have_tails;

    /**
     * @brief Marks the units in the circuit as accessible from the unit_num and
     *        marks whether the units have a route forward to both of the outlet streams.
     *
     * This recursive function marks the units in the circuit as accessible from the unit_num
     * and marks the uniis have a route forward to both of the outlet streams
     * by traversing the circuit graph starting from the given unit number.
     *
     * @param unit_num The unit number to start the traversal from.
     */
    void mark_accessible(int unit_num);
};

