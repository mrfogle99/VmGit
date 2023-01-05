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
    std::vector<double> result;

    pthread_mutex_t* lock;
    std::vector<Vector>* datasets;
};

std::vector<double> paracorrelation_coefficients(std::vector<Vector> datasets, unsigned threadCount)
{
    pthread_t threads[threadCount];
    pearson_thread_data data [threadCount];

    int rows = datasets.size() / threadCount;

    pthread_mutex_t lock;
    pthread_mutex_init(&lock, NULL);

    std::vector<double> result {};

    for(int i = 0; i < threadCount; i++)
    {
        /*
        unsigned l = length;
        unsigned end = (i * length) + length;

        if(end > totalLength)
        {
            l = totalLength - (i * length);
        }
        */

        data[i].result = std::vector<double>();
        data[i].lock = &lock;
        data[i].startRow = i * rows;
        data[i].rowsCount = rows;
        data[i].datasets = &datasets;
        
        pthread_create(&threads[i], NULL, calc_coefficients, &data[i]);
    }

    for(int i = 0; i < threadCount; i++)
        pthread_join(threads[i], NULL);

    for(int i = 0; i < threadCount; i++)
        result.insert(result.end(), data[i].result.begin(), data[i].result.end());

    return result;
}

void* calc_coefficients (void* v)
{
    pearson_thread_data* data = static_cast<pearson_thread_data*>(v);
 
    int max = data->startRow + data->rowsCount;
    
    for(int sample1 = data->startRow; sample1 < max - 1; sample1++)
    {
        for(int sample2 = (sample1 + 1); sample2 < max; sample2++)
        {
            double corr = pearson((*data->datasets)[sample1], (*data->datasets)[sample2]);
            pthread_mutex_lock(data->lock);

            data->result.push_back(corr);
            
            pthread_mutex_unlock(data->lock);
        }
    }
}

std::vector<double> correlation_coefficients(std::vector<Vector> datasets)
{
    std::vector<double> result {};

    for (int sample1 = 0; sample1 < datasets.size() - 1; sample1++)
    {
        for (int sample2 = (sample1 + 1); sample2 < datasets.size(); sample2++) 
        {
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
