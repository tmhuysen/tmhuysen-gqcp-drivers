/**
 *  Test driver
 */

#include <gqcp.hpp>
#include <hdf5.h>

int main() {
    GQCP::Atom atom1 (GQCP::elements::elementToAtomicNumber("N"), -1/2, 0, 0);
    std::cout << "Hello world";
    return 0;
}
