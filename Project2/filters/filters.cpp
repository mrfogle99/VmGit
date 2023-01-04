/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <cmath>
#include <pthread.h>

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

struct blur_thread_data 
{
    unsigned index;
    unsigned length;
    unsigned radius;
    Matrix* mat;
};

Matrix parablur (Matrix m, const int radius, const int maxThreads)
{
    Matrix scratch { PPM::max_dimension };
    Matrix dst { m };

    unsigned totalLength = dst.get_x_size() * dst.get_y_size();
    unsigned length = totalLength / maxThreads;


    pthread_t threads [maxThreads];
    blur_thread_data ranges[maxThreads];

    for(int i = 0; i < maxThreads; i++)
    {
        unsigned l = length;
        unsigned end = (i * length) + length;
        if(end > totalLength)
            l = totalLength - (i * length);

        ranges[i].mat = &dst;
        ranges[i].length = l;
        ranges[i].index = i * length;
        ranges[i].radius = radius;
        pthread_create(&threads[i], NULL, calc_blur, (void*)&ranges[i]);
    }

    for(int i = 0; i < maxThreads; i++)
    {
        pthread_join(threads[i], NULL);
    }

    return dst;
}

void* calc_blur (void* args)
{
    blur_thread_data* range = static_cast<blur_thread_data*>(args);
    Matrix scratch { PPM::max_dimension };

    double w[Gauss::max_radius] {};
    Gauss::get_weights(range->radius, w);


    for(int i = 0; i < range->length; i++)
    {
        int j = (i + range->index);
        int x = j / range->mat->get_y_size();
        int y = j % range->mat->get_y_size();

        auto r { w[0] * range->mat->r(x, y) };
        auto g { w[0] * range->mat->g(x, y) };
        auto b { w[0] * range->mat->b(x, y) };
        auto n { w[0] };

        for (auto wi { 1 }; wi <= range->radius; wi++) {
            auto wc { w[wi] };
            auto x2 { x - wi };
            if (x2 >= 0) {
                r += wc * range->mat->r(x2, y);
                g += wc * range->mat->g(x2, y);
                b += wc * range->mat->b(x2, y);
                n += wc;
            }
            x2 = x + wi;
            if (x2 < range->mat->get_x_size()) {
                r += wc * range->mat->r(x2, y);
                g += wc * range->mat->g(x2, y);
                b += wc * range->mat->b(x2, y);
                n += wc;
            }
        }
        scratch.r(x, y) = r / n;
        scratch.g(x, y) = g / n;
        scratch.b(x, y) = b / n;
    }

    for(int i = 0; i < range->length; i++)
    {
        int j = (i + range->index);
        int x = j / range->mat->get_y_size();
        int y = j % range->mat->get_y_size();

        auto r { w[0] * scratch.r(x, y) };
        auto g { w[0] * scratch.g(x, y) };
        auto b { w[0] * scratch.b(x, y) };
        auto n { w[0] };

        for (auto wi { 1 }; wi <= range->radius; wi++) {
            auto wc { w[wi] };
            auto y2 { y - wi };
            if (y2 >= 0) {
                r += wc * scratch.r(x, y2);
                g += wc * scratch.g(x, y2);
                b += wc * scratch.b(x, y2);
                n += wc;
            }
            y2 = y + wi;
            if (y2 < range->mat->get_y_size()) {
                r += wc * scratch.r(x, y2);
                g += wc * scratch.g(x, y2);
                b += wc * scratch.b(x, y2);
                n += wc;
            }
        }

        range->mat->r(x, y) = r / n;
        range->mat->g(x, y) = g / n;
        range->mat->b(x, y) = b / n;
    }
}

struct threshold_thread_data 
{
    unsigned index;
    unsigned length;

    unsigned* sum;
    Matrix* mat;
    pthread_barrier_t* barr;
};

Matrix parathreshold (Matrix m, int maxThreads)
{
    Matrix scratch { PPM::max_dimension };
    Matrix dst { m };

    unsigned totalLength = dst.get_x_size() * dst.get_y_size();
    unsigned length = totalLength / maxThreads;

    pthread_t threads [maxThreads];
    threshold_thread_data data [maxThreads];

    pthread_barrier_t barr;
    pthread_barrier_init(&barr, NULL, maxThreads);

    unsigned sum = 0;

    for (auto i { 0 }; i < totalLength; i++) {
        sum += dst.r(i, 0) + dst.g(i, 0) + dst.b(i, 0);
    }

    sum /= totalLength;



    sum = 0;

    for(int i = 0; i < maxThreads; i++)
    {
        unsigned l = length;
        unsigned end = (i * length) + length;
        if(end > totalLength)
        {
            l = totalLength - (i * length);
            printf("new length %i\n", l);
        }

        data[i].mat = &dst;
        data[i].length = l;
        data[i].index = i * length;
        data[i].barr = &barr;
        data[i].sum = &sum;
        pthread_create(&threads[i], NULL, calcthreshold, (void*)&data[i]);
    }

    for(int i = 0; i < maxThreads; i++)
    {
        pthread_join(threads[i], NULL);
    }

    pthread_barrier_destroy(&barr);
    return dst;
}

void* calcthreshold (void* args)
{
    threshold_thread_data* data = static_cast<threshold_thread_data*>(args);
  
    for(int j = 0; j < data->length; j++)
    {
        int i = (j + data->index);  
        unsigned val = data->mat->r(i, 0) + data->mat->g(i, 0) + data->mat->b(i, 0); 
        *data->sum += val;
    }

    pthread_barrier_wait(data->barr);

    int sum = *data->sum;
    sum /= (data->mat->get_x_size() * data->mat->get_y_size());

    unsigned psum {};

    for (int j = 0; j < data->length; j++) 
    {
        int i = data->index + j;
        psum = data->mat->r(i, 0) + data->mat->g(i, 0) + data->mat->b(i, 0);

        if (sum > psum) {
            data->mat->r(i, 0) = data->mat->g(i, 0) = data->mat->b(i, 0) = 0;
        } else {
            data->mat->r(i, 0) = data->mat->g(i, 0) = data->mat->b(i, 0) = 255;
        }
    }

}

Matrix threshold(Matrix m)
{
    auto dst { m };

    unsigned sum = 0;
    unsigned area { dst.get_x_size() * dst.get_y_size() };

    for (auto i { 0 }; i < area; i++) {
        sum += dst.r(i, 0) + dst.g(i, 0) + dst.b(i, 0);
    }

    sum /= area;



    unsigned psum {};

    for (auto i { 0 }; i < area; i++) {
        psum = dst.r(i, 0) + dst.g(i, 0) + dst.b(i, 0);
        if (sum > psum) {
            dst.r(i, 0) = dst.g(i, 0) = dst.b(i, 0) = 0;
        } else {
            dst.r(i, 0) = dst.g(i, 0) = dst.b(i, 0) = 255;
        }
    }

    return dst;
}

}
