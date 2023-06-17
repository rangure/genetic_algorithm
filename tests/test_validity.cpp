#include <iostream>
#include <assert.h>

#include "CUnit.h"
#include "CCircuit.h"
#include <cstddef>

/** 
 * \brief Function to test circuit validity 
 *
 * \param vec Vector of integers representing a circuit.
 * \param expected_validity Expected result of circuit validity check.
 * 
 * This function creates a CCircuit object, checks the circuit validity and 
 * validates the result using assert.
 */
template<std::size_t N>
void test_circuit(int (&vec)[N], bool expected_validity) {
    CCircuit circuit(N);
    assert(circuit.check_validity(vec) == expected_validity);
}

/** 
 * \brief Main function
 *
 * This function includes test cases to check the validity of a circuit. 
 * Tests are conducted on both valid and invalid circuit configurations.
 */
int main(int argc, char * argv[]){

    // Test 0: valid circuit
    int valid[3] = {0, 1, 2};
    test_circuit(valid, true);

    // Test 1
    int vec1[] = {0, 1, 2, 3, 4, 4, 5, 6, 7, 7, 8, 8, 9,
				10, 11, 10, 11, 10, 11, 10, 11};
    test_circuit(vec1, true);

    // Test 2
    int vec2[] = {0, 1, 11, 2, 11, 3, 11, 4, 11, 5, 11, 6, 11,
				7, 11, 8, 11, 9, 11, 10, 11};
    test_circuit(vec2, true);


    // Test 3: valid, but need to traverse twice to find all outlets
    int vec3[] = {0, 1, 2, 3, 0, 0, 4};
    test_circuit(vec3, true);

    // Test 4:
    int vec4[] = {4, 5, 1, 2, 4, 0, 1, 1, 6, 1, 3};
    test_circuit(vec4, true);

    // Invalid circuits
    // Test 1: Every unit must be accessible from the feed.
    int invalid_accessibility[11] = {0, 1, 2, 3, 5, 1, 6, 5, 0, 1, 2};
    test_circuit(invalid_accessibility, false);

    // Test 2: Every unit must have a route forward to both of the outlet streams.
    // no one of outlet
    int invalid_final_outlet[7] = {0, 1, 2, 0, 4, 0, 4};
    test_circuit(invalid_final_outlet, false);

    // Test 3: outlet not reachable
    int invalid_final_outlet2[11] = {0, 1, 4, 2, 3, 1, 3, 1, 2, 5, 6};
    test_circuit(invalid_final_outlet2, false);

    int invalid_final_outlet3[] = {0, 1, 2, 5, 0, 3, 4, 2, 6, 2, 6};
    test_circuit(invalid_final_outlet3, false);

    // Test 4: There should be no self-recycle.
    int invalid_recycle[3] = {0, 1, 0};
    test_circuit(invalid_recycle, false);

    // test 5: The destination for both products from a unit should not be the same unit.
    int invalid_same_destination[3] = {0, 1, 1};
    test_circuit(invalid_same_destination, false);

    // Test 6: The destination id of each unit must be between 0 to num_units + 1.
    int invalid_destination[3] = {0, 1, 3};
    test_circuit(invalid_destination, false);

    // test 7: The first destination must not be tailings, the second destination must not be concentrate.
    int invalid_destination2[3] = {0, 2, 1};
    test_circuit(invalid_destination2, false);

    return 0;
}