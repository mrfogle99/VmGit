/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "analysis.hpp"
#include "dataset.hpp"
#include <iostream>
#include <cstdlib>

int main(int argc, char const* argv[])
{
    #ifndef MAKE_BASERINO
        if (argc != 4) {
            std::cerr << "Usage: " << argv[0] << " [dataset] [outfile] [threads]" << std::endl;
            std::exit(1);
        }
    
        unsigned threads = static_cast<unsigned>(std::stoul(argv[3]));

    #else
        if (argc != 3) {
            std::cerr << "Usage: " << argv[0] << " [dataset] [outfile]" << std::endl;
            std::exit(1);
        }
    #endif

    std::vector<Vector> datasets = Dataset::read(argv[1]);


#ifndef MAKE_BASERINO
    std::vector<double> result_corrs = Analysis::paracorrelation_coefficients(datasets, threads);
#else
    std::vector<double> result_corrs = Analysis::correlation_coefficients(datasets);
#endif
    

    printf("rows GOTTEM: %llu\n", result_corrs.size());

    Dataset::write(result_corrs, argv[2]);

    return 0;
}
