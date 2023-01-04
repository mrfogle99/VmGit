/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "matrix.hpp"

#if !defined(FILTERS_HPP)
#define FILTERS_HPP

namespace Filter {

namespace Gauss {
    constexpr unsigned max_radius { 1000 };
    constexpr float max_x { 1.33 };
    constexpr float pi { 3.14159 };
    
    void get_weights(int n, double* weights_out);
}


Matrix parablur (Matrix m, const int radius, const int maxThreads);
void* calc_blur(void*);

Matrix parathreshold (Matrix m, int maxThreads);
Matrix threshold(Matrix m);
void* calcthreshold (void* args);

}

#endif