/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "matrix.hpp"
#include "ppm.hpp"
#include <fstream>
#include <stdexcept>
#include <cstring>

Matrix::Matrix(unsigned char* R, unsigned char* G, unsigned char* B, unsigned x_size, unsigned y_size, unsigned color_max)
    : R { R }
    , G { G }
    , B { B }
    , x_size { x_size }
    , y_size { y_size }
    , color_max { color_max }
{
}

Matrix::Matrix()
    : Matrix {
        nullptr,
        nullptr,
        nullptr,
        0,
        0,
        0,
    }
{
}

Matrix::Matrix(unsigned dimension)
    : R { new unsigned char[dimension * dimension] }
    , G { new unsigned char[dimension * dimension] }
    , B { new unsigned char[dimension * dimension] }
    , x_size { dimension }
    , y_size { dimension }
    , color_max { 0 }
{
}

Matrix::Matrix(const Matrix& other)
    : R { new unsigned char[other.x_size * other.y_size] }
    , G { new unsigned char[other.x_size * other.y_size] }
    , B { new unsigned char[other.x_size * other.y_size] }
    , x_size { other.x_size }
    , y_size { other.y_size }
    , color_max { other.color_max }
{
    const int size =other.x_size * other.y_size;
    std::memcpy(R, other.R, size);
    std::memcpy(G, other.G, size);
    std::memcpy(B, other.B, size);

    // changed to single loop
    // linear
    #pragma omp parallel for
    for(int i = 0; i< (x_size*y_size); i++)
    {
        R[i] = other.R[i];
        G[i] = other.G[i];
        B[i] = other.B[i];
    }
}

Matrix& Matrix::operator=(const Matrix other)
{
    if (this == &other) {
        return *this;
    }

    this->~Matrix();
    const int size =other.x_size * other.y_size;
    R = new unsigned char[size];
    G = new unsigned char[size];
    B = new unsigned char[size];

    x_size = other.x_size;
    y_size = other.y_size;
    color_max = other.color_max;

    /*
    std::memcpy(R, other.R, size);
    std::memcpy(G, other.G, size);
    std::memcpy(B, other.B, size);
    */


    // changed to single loop
    // linear
    #pragma omp parallel for
    for(int i = 0; i< (x_size*y_size); i++)
    {
        R[i] = other.R[i];
        G[i] = other.G[i];
        B[i] = other.B[i];
    }

    return *this;
}

Matrix::~Matrix()
{
    if (R) {
        delete[] R;
        R = nullptr;
    }
    if (G) {
        delete[] G;
        G = nullptr;
    }
    if (B) {
        delete[] B;
        B = nullptr;
    }

    x_size = y_size = color_max = 0;
}

unsigned Matrix::get_x_size() const
{
    return x_size;
}

unsigned Matrix::get_y_size() const
{
    return y_size;
}

unsigned Matrix::get_color_max() const
{
    return color_max;
}

unsigned char* Matrix::get_R()
{
    return R;
}

unsigned char * Matrix::get_G()
{
    return G;
}

unsigned char * Matrix::get_B()
{
    return B;
}

unsigned char Matrix::r(unsigned x, unsigned y) const
{
    return R[y * x_size + x];
}

unsigned char Matrix::g(unsigned x, unsigned y) const
{
    return G[y * x_size + x];
}

unsigned char Matrix::b(unsigned x, unsigned y) const
{
    return B[y * x_size + x];
}

unsigned char& Matrix::r(unsigned x, unsigned y)
{
    return R[y * x_size + x];
}

unsigned char& Matrix::g(unsigned x, unsigned y)
{
    return G[y * x_size + x];
}

unsigned char& Matrix::b(unsigned x, unsigned y)
{
    return B[y * x_size + x];
}
