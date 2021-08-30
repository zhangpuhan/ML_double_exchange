#include <iostream>
#include "lattice_base.h"
#include "analysis.h"

int main(int argc, char* argv[]) {

    std::string status = argc > 1 ? argv[1] : "finalpass";
    int lattice_length = argc > 2 ? atoi(argv[2]) : 30;
    double h = argc > 3 ? atof(argv[3]) : 0.01;
    double alpha = argc > 4 ? atof(argv[4]) : 1.0;
    int round = argc > 5 ? atoi(argv[5]) : 0;

    std::vector<Square_Site> test_lattice;
    test_lattice.reserve(lattice_length * lattice_length);
    for (int i = 0; i < lattice_length * lattice_length; ++i) {
        test_lattice.emplace_back(i, lattice_length);
    }

    if (status == "initialize") {
        initialize(test_lattice);
        return 1;
    }
    if (status == "intermediate") {
        intermediate_passing(test_lattice, h, alpha);
        return 2;
    }
    if (status == "finalpass") {
        final_passing(test_lattice, h, alpha, round);
        return 3;
    }


    return 0;
}