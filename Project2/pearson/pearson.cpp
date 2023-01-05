/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "analysis.hpp"
#include "dataset.hpp"
#include <iostream>
#include <cstdlib>

int main(int argc, char const* argv[])
{
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " [dataset] [outfile] [threads]" << std::endl;
        std::exit(1);
    }

    unsigned threads = static_cast<unsigned>(std::stoul(argv[3]));

    std::vector<Vector> datasets = Dataset::read(argv[1]);

#ifndef DMAKE_BASERINO
    std::vector<double> result_corrs = Analysis::paracorrelation_coefficients(datasets, threads);
#else
    std::vector<double> result_corrs = Analysis::correlation_coefficients(datasets);
#endif
    
    Dataset::write(result_corrs, argv[2]);

    return 0;
}
