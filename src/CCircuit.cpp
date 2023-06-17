#include <vector>
#include <iostream>
#include <stdio.h>

#include "CUnit.h"
#include "CCircuit.h"

/** 
 * \brief CCircuit Constructor.
 *
 * This function initializes the number of units in the circuit and 
 * dynamically allocates memory for marks array.
 * 
 * \param vector_size The size of the circuit vector.
 */
CCircuit::CCircuit(int vector_size){
//   this->units.resize(num_units);
  this->num_units = (vector_size - 1) / 2;
  this->marks = new bool[this->num_units];
  this->have_conc = new bool[this->num_units]();
  this->have_tails = new bool[this->num_units]();
}

/** 
 * \brief CCircuit Destructor.
 *
 * This function deletes the dynamically allocated memory for marks array.
 */
CCircuit::~CCircuit()
{
    delete[] this->marks;
    delete[] this->have_conc;
    delete[] this->have_tails;
}

/** 
 * \brief Check if the circuit is valid.
 *
 * This function checks for various conditions for circuit validity, 
 * including unit accessibility from the feed, route forward to both of the outlet streams,
 * no self-recycle, different destinations, and appropriate destination ID.
 *
 * \param circuit_vector A pointer to an array representing the circuit.
 * 
 * \return true if the circuit is valid, false otherwise.
 */
bool CCircuit::check_validity(int* circuit_vector)
{
    this->circuit_vector = circuit_vector;

    for (int i = 0; i < this->num_units; i++) {
        this->marks[i] = false;
        this->have_conc[i] = false;
        this->have_tails[i] = false;
    }

    // 1. Check worong cases
    for (int i = 0; i < this->num_units; i++) {
        int conc_num = circuit_vector[2 * i + 1];
        int tails_num = circuit_vector[2 * i + 2];

        // 3. Check no self-recycle
        if (conc_num == i || tails_num == i)
            return false;

        // 4. Check destinations not the same unit
        if (conc_num == tails_num)
            return false;

        // 5. Check destination id between 0 and num_units + 1
        if (conc_num < 0 || conc_num > this->num_units + 1)
            return false;
        if (conc_num == this->num_units + 1)  // check first destination is not tailings
            return false;

        if (tails_num < 0 || tails_num > this->num_units + 1)
            return false;
        if (tails_num == this->num_units) // check second destination is not concentrate
            return false;
    }

    // 2. Check being accessible from the feed
    mark_accessible(circuit_vector[0]);
    for (int i = 0; i < this->num_units; i++) {
        if (this->marks[i] == false)
            return false;
    }

    // 3. Check having a route forward to both of the outlet streams
    for (int i = 0; i < this->num_units; i++) {
        if (this->have_conc[i] == true && this->have_tails[i] == true)
            continue;

        int conc_num = circuit_vector[2 * i + 1];
        int tails_num = circuit_vector[2 * i + 2];

        for (int j = 0; j < this->num_units; j++)
            this->marks[j] = false;
        this->marks[i] = true;

        while (this->have_conc[i] != true) {
            if (this->marks[conc_num] == true)
                return false;

            if (this->have_conc[conc_num] == true)
                this->have_conc[i] = true;
            else {
                this->marks[conc_num] = true;
                conc_num = circuit_vector[2 * conc_num + 1];
            }
        }

        for (int j = 0; j < this->num_units; j++)
            this->marks[j] = false;
        this->marks[i] = true;

        while (this->have_tails[i] != true) {
            if (this->marks[tails_num] == true)
                return false;

            if (this->have_tails[tails_num] == true)
                this->have_tails[i] = true;
            else {
                this->marks[tails_num] = true;
                tails_num = circuit_vector[2 * tails_num + 2];
            }
        }
    }

    return true;
}

/** 
 * \brief Marks units that are accessible from a given unit.
 *
 * This function is used to determine which units in the circuit can be reached 
 * from a given unit. It uses recursion to visit all connected units.
 *
 * \param unit_num The number of the unit from which we want to determine accessibility.
 */
void CCircuit::mark_accessible(int unit_num)
{
    if (this->marks[unit_num]) return;  // Exit if we have seen this unit already

    this->marks[unit_num] = true;   // Mark that we have now seen the unit

    int conc_num = this->circuit_vector[2 * unit_num + 1];
    int tails_num = this->circuit_vector[2 * unit_num + 2];

    //If conc_num does not point at a circuit outlet recursively call the function
    if (conc_num < this->num_units) {
        mark_accessible(conc_num);
    } else {
        // ...Potentially do something to indicate that you have seen an exit
        if (conc_num == num_units)
            this->have_conc[unit_num] = true;
    }

    //If tails_num does not point at a circuit outlet recursively call the function 
    if (tails_num < this->num_units) {
        mark_accessible(tails_num);
    } else {
        // ...Potentially do something to indicate that you have seen an exit
        if (tails_num == num_units + 1)
            this->have_tails[unit_num] = true;
    }
}

