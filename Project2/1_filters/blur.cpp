/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "matrix.hpp"
#include "ppm.hpp"
#include "filters.hpp"
#include <cstdlib>
#include <iostream>

int main(int argc, char const* argv[])
{
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " [radius] [infile] [outfile] [threads]" << std::endl;
        std::exit(1);
    }

    PPM::Reader reader {};
    PPM::Writer writer {};

    Matrix m { reader(argv[2]) };
    unsigned radius { static_cast<unsigned>(std::stoul(argv[1])) };
    unsigned threads { static_cast<unsigned>(std::stoul(argv[4])) };
    

    Matrix blurred { Filter::parablur(m, radius, threads) };
    writer(blurred, argv[3]);

    return 0;
}