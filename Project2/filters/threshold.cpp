/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <cstdlib>
#include <iostream>

int main(int argc, char const* argv[])
{
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " [infile] [outfile]" << std::endl;
        std::exit(1);
    }

    unsigned threads { static_cast<unsigned>(std::stoul(argv[3])) };


    PPM::Reader reader {};
    PPM::Writer writer {};

    auto m { reader(argv[1]) };
    auto thresholded { Filter::parathreshold(m, threads) };

    writer(thresholded, argv[2]);

    return 0;
}
