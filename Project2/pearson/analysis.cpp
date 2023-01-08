/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "analysis.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
#include <vector>
#include <pthread.h>

namespace Analysis {

struct pearson_thread_data 
{
    unsigned startRow;
    unsigned rowsCount;
    unsigned cap;

    unsigned* iter;
    std::vector<double>* result;
    
    pthread_mutex_t* lock;
    std::vector<Vector>* datasets;
};

std::vector<double> paracorrelation_coefficients(std::vector<Vector> datasets, unsigned threadCount)
{
    pthread_t threads[threadCount];
    pearson_thread_data data [threadCount];

    int rows = datasets.size();

    pthread_mutex_t lock;
    pthread_mutex_init(&lock, NULL);

    std::vector<double> result {};
    unsigned iter = 0;

    for(int i = 0; i < threadCount; i++)
    {
        data[i].result = &result;
        data[i].lock = &lock;

        data[i].rowsCount = (rows / threadCount);
        data[i].startRow = i * data[i].rowsCount;

        data[i].datasets = &datasets;
        data[i].iter = &iter;
        data[i].cap = rows;
        
        printf("starting %i range(%i-%i) rows(%i)\n", i, (i*data[i].rowsCount), ((i*rows)+rows-1), rows - data[i].startRow);
        pthread_create(&threads[i], NULL, calc_coefficients, &data[i]);
    }

    for(int i = 0; i < threadCount; i++)
        pthread_join(threads[i], NULL);

    //for(int i = 0; i < threadCount; i++)
     //   result.insert(result.end(), data[i].result.begin(), data[i].result.end());

    return result;
}

void* calc_coefficients (void* v)
{
    pearson_thread_data* data = static_cast<pearson_thread_data*>(v);
 
    const unsigned max = data->startRow + data->rowsCount;

    for(unsigned sample1 = data->startRow; sample1 < max; sample1++)
    {
        for(unsigned sample2 = (sample1 + 1); sample2 < data->cap; sample2++)
        {
            //__atomic_fetch_add(data->iter, 1, __ATOMIC_SEQ_CST);

            double corr = pearson(data->datasets->at(sample1), data->datasets->at(sample2));
            
            pthread_mutex_lock(data->lock);
            data->result->push_back(corr);
            pthread_mutex_unlock(data->lock);
        }
    }
    
}

std::vector<double> correlation_coefficients(std::vector<Vector> datasets)
{
    std::vector<double> result {};
    unsigned iter = 0;

    for (int sample1 = 0; sample1 < datasets.size() - 1; sample1++)
    {
        for (int sample2 = (sample1 + 1); sample2 < datasets.size(); sample2++) 
        {
            iter ++;
            double corr = pearson(datasets[sample1], datasets[sample2]);
            result.push_back(corr);
        }
    }

    return result;
}

double pearson(Vector vec1, Vector vec2)
{
    double x_mean = vec1.mean();
    double y_mean = vec2.mean();

    Vector x_mm = vec1 - x_mean;
    Vector y_mm = vec2 - y_mean;

    double x_mag = x_mm.magnitude();
    double y_mag = y_mm.magnitude();

    Vector x_mm_over_x_mag = x_mm / x_mag;
    Vector y_mm_over_y_mag = y_mm / y_mag;

    double r { x_mm_over_x_mag.dot(y_mm_over_y_mag) };

    return std::max(std::min(r, 1.0), -1.0);
}
};
