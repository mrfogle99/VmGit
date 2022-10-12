/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <cmath>

namespace Filter {

namespace Gauss {
    void get_weights(int n, double* weights_out)
    {
        for (auto i { 0 }; i <= n; i++) {
            double x { static_cast<double>(i) * max_x / n };
            weights_out[i] = exp(-x * x * pi);
        }
    }
}

Matrix blur(Matrix m, const int radius)
{
    Matrix scratch { PPM::max_dimension };
    auto dst { m };

    // Moved the gauss weight calculations from the double loop to just doing it a single time.
    double w[Gauss::max_radius] {};
    Gauss::get_weights(radius, w);

    for (auto x { 0 }; x < dst.get_x_size(); x++) {
        for (auto y { 0 }; y < dst.get_y_size(); y++) {


            auto r { w[0] * dst.r(x, y) }, g { w[0] * dst.g(x, y) }, b { w[0] * dst.b(x, y) }, n { w[0] };

            for (auto wi { 1 }; wi <= radius; wi++) {
                auto wc { w[wi] };
                auto x2 { x - wi };
                if (x2 >= 0) {
                    r += wc * dst.r(x2, y);
                    g += wc * dst.g(x2, y);
                    b += wc * dst.b(x2, y);
                    n += wc;
                }
                x2 = x + wi;
                if (x2 < dst.get_x_size()) {
                    r += wc * dst.r(x2, y);
                    g += wc * dst.g(x2, y);
                    b += wc * dst.b(x2, y);
                    n += wc;
                }
            }
            scratch.r(x, y) = r / n;
            scratch.g(x, y) = g / n;
            scratch.b(x, y) = b / n;
        }
    }

    for (auto x { 0 }; x < dst.get_x_size(); x++) {
        for (auto y { 0 }; y < dst.get_y_size(); y++) {

            auto r { w[0] * scratch.r(x, y) }, g { w[0] * scratch.g(x, y) }, b { w[0] * scratch.b(x, y) }, n { w[0] };

            for (auto wi { 1 }; wi <= radius; wi++) {
                auto wc { w[wi] };
                auto y2 { y - wi };
                if (y2 >= 0) {
                    r += wc * scratch.r(x, y2);
                    g += wc * scratch.g(x, y2);
                    b += wc * scratch.b(x, y2);
                    n += wc;
                }
                y2 = y + wi;
                if (y2 < dst.get_y_size()) {
                    r += wc * scratch.r(x, y2);
                    g += wc * scratch.g(x, y2);
                    b += wc * scratch.b(x, y2);
                    n += wc;
                }
            }
            dst.r(x, y) = r / n;
            dst.g(x, y) = g / n;
            dst.b(x, y) = b / n;
        }
    }

    return dst;
}

Matrix threshold(Matrix m)
{
    auto dst { m };
    unsigned nump { dst.get_x_size() * dst.get_y_size() };
    unsigned sum {};

    // Gets a pointer to the matrix data
    auto R = dst.get_R();
    auto G = dst.get_G();
    auto B = dst.get_B();

    for (auto i { 0 }; i < nump; i++)
    {
        sum += R[i] + G[i] + B[i];
    }

    sum /= nump;

    unsigned psum = 0;

    /*
     * Removed the potential branching by replacing the if statement with a simple boolean check
     * This will then be implicitly cast into a integer, this gives the threshold values, or 'd' an value of either 0 or 255
     * The matrix data is also manipulated directly without the need of calling it's set functions. A small performance increase.
     */
    for (auto i { 0 }; i < nump; i++) {
        psum = R[i] + G[i] + B[i];

        // removed if statement
        const unsigned char d = (unsigned char)((sum <= psum) * 255);
        R[i] = G[i] = B[i] = d;
    }

    return dst;
}

}
