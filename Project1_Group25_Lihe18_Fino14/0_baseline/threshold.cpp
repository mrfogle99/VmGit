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
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " [infile] [outfile]" << std::endl;
        std::exit(1);
    }
for(int i = 0; i < 20; i++) {

    PPM::Reader reader {};
    PPM::Writer writer {};

    auto m { reader(argv[1]) };
    auto thresholded { Filter::threshold(m) };

    writer(thresholded, argv[2]);
}

    return 0;
}
