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

    pthread_barrier_t* barr;
    pthread_mutex_t* lock;
    Matrix* scratch;
};

Matrix parablur (Matrix m, const int radius, const int maxThreads)
{
    Matrix scratch { PPM::max_dimension };
    Matrix dst { m };

    unsigned totalLength = dst.get_x_size() * dst.get_y_size();
    unsigned length = totalLength / maxThreads;

    pthread_t threads [maxThreads];
    blur_thread_data ranges[maxThreads];

    pthread_mutex_t lock;
    pthread_mutex_init(&lock, NULL);

    pthread_barrier_t barr;
    pthread_barrier_init(&barr, NULL, maxThreads);

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
        ranges[i].lock = &lock;
        ranges[i].scratch = &scratch;
        ranges[i].barr = &barr;

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
    blur_thread_data* data = static_cast<blur_thread_data*>(args);
    double w[Gauss::max_radius] {};
    Gauss::get_weights(data->radius, w);

    for(int i = 0; i < data->length; i++)
    {
        int j = (i + data->index);
        int x = j % data->mat->get_x_size();
        int y = j / data->mat->get_x_size();

        auto r { w[0] * data->mat->r(x, y) };
        auto g { w[0] * data->mat->g(x, y) };
        auto b { w[0] * data->mat->b(x, y) };
        auto n { w[0] };

        for (auto wi { 1 }; wi <= data->radius; wi++) {
            auto wc { w[wi] };
            auto x2 { x - wi };
            if (x2 >= 0) {
                r += wc * data->mat->r(x2, y);
                g += wc * data->mat->g(x2, y);
                b += wc * data->mat->b(x2, y);
                n += wc;
            }
            x2 = x + wi;
            if (x2 < data->mat->get_x_size()) {
                r += wc * data->mat->r(x2, y);
                g += wc * data->mat->g(x2, y);
                b += wc * data->mat->b(x2, y);
                n += wc;
            }
        }

        pthread_mutex_lock(data->lock);
        data->scratch->r(x, y) = r / n;
        data->scratch->g(x, y) = g / n;
        data->scratch->b(x, y) = b / n;
        pthread_mutex_unlock(data->lock);
    }

    pthread_barrier_wait(data->barr);

    for(int i = 0; i < data->length; i++)
    {
        int j = (i + data->index);
        int x = j % data->mat->get_x_size();
        int y = j / data->mat->get_x_size();

        auto r { w[0] * data->scratch->r(x, y) };
        auto g { w[0] * data->scratch->g(x, y) };
        auto b { w[0] * data->scratch->b(x, y) };
        auto n { w[0] };

        for (auto wi { 1 }; wi <= data->radius; wi++) {
            auto wc { w[wi] };
            auto y2 { y - wi };
            if (y2 >= 0) {
                r += wc * data->scratch->r(x, y2);
                g += wc * data->scratch->g(x, y2);
                b += wc * data->scratch->b(x, y2);
                n += wc;
            }
            y2 = y + wi;
            if (y2 < data->mat->get_y_size()) {
                r += wc * data->scratch->r(x, y2);
                g += wc * data->scratch->g(x, y2);
                b += wc * data->scratch->b(x, y2);
                n += wc;
            }
        }

        data->mat->r(x, y) = r / n;
        data->mat->g(x, y) = g / n;
        data->mat->b(x, y) = b / n;
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
/*
    unsigned sum = 0;
    unsigned area { dst.get_x_size() * dst.get_y_size() };

    for (auto i { 0 }; i < area; i++) {
        sum += dst.r(i, 0) + dst.g(i, 0) + dst.b(i, 0);
    }

    printf("og sum: %u over area %u\n", sum, area);
    sum /= area;
    printf("og sum after: %u\n", sum);
*/

    // im4 ger 444, im1 ger 326 t.ex

    unsigned totalLength = dst.get_x_size() * dst.get_y_size();
    unsigned length = totalLength / (unsigned)maxThreads;

    pthread_t threads [maxThreads];
    threshold_thread_data data [maxThreads];

    pthread_barrier_t barr;
    pthread_barrier_init(&barr, NULL, maxThreads);

    unsigned sum = 0;

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

        printf("starting thread. index(%u) length(%u)\n", data[i].index, data[i].length);
        pthread_create(&threads[i], NULL, calcthreshold, &data[i]);
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
    Matrix* mat = data->mat;

    for(unsigned j = 0; j < data->length; j++)
    {
        const unsigned i = (j + data->index);  
        unsigned val = mat->r(i, 0) + mat->g(i, 0) + mat->b(i, 0); 
        __atomic_fetch_add(data->sum, val, __ATOMIC_SEQ_CST);
    }

    pthread_barrier_wait(data->barr);
    unsigned sum = *data->sum / (mat->get_x_size() * mat->get_y_size());
    unsigned psum = 0;

    for (int j = 0; j < data->length; j++) 
    {
        int i = data->index + j;
        psum = mat->r(i, 0) + mat->g(i, 0) + mat->b(i, 0);

        if (sum > psum) {
            mat->r(i, 0) = mat->g(i, 0) = mat->b(i, 0) = 0;
        } else {
            mat->r(i, 0) = mat->g(i, 0) = mat->b(i, 0) = 255;
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
