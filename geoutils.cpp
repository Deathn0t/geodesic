#pragma once

#include "geoutils.h"


const double EPS = 1e-7;

Vector2d intersect(Vector2d u, Vector2d v)
{
    Vector2d inter;
    // x coord. of the intersection between l1 and l2
    inter(0) = (u(1) - v(1)) / (v(0) - u(0));
    // y coord. of the intersection between l1 and l2
    inter(1) = v(0) * inter(0) + v(1);
    return inter;
}

/**
 * Check if c is in range [a,b]
**/
bool point_in_range(Vector2d c, Vector2d a, Vector2d b)
{
    double x0, x1, y0, y1;
    if (a(0) > b(0))
    {
        x0 = b(0);
        x1 = a(0);
    }
    else
    {
        x0 = a(0);
        x1 = b(0);
    }

    if (a(1) > b(1))
    {
        y0 = b(1);
        y1 = a(1);
    }
    else
    {
        y0 = a(1);
        y1 = b(1);
    }

    if (x0 <= EPS)
    {
        x0 = 0.0;
    }

    if (x1 <= EPS)
    {
        x1 = 0.0;
    }

    if (y0 <= EPS)
    {
        y0 = 0.0;
    }

    if (y1 <= EPS)
    {
        y1 = 0.0;
    }

    return (x0 <= c(0) && c(0) <= x1 && y0 <= c(1) && c(1) <= y1);
}