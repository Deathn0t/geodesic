#include "window.h"

Window::Window() {}

Window::Window(double f_b0, double f_b1, double f_d0, double f_d1, double f_sigma, int f_dir, int f_edge_id, Vector3d &f_v0, Vector3d &f_v1, int f_v0id, int f_v1id)
{
    if (f_b0 < 1e-2)
    {
        b0 = 0.0;
    }
    else
    {
        b0 = f_b0;
    }
    b1 = f_b1;
    d0 = f_d0;
    d1 = f_d1;
    sigma = f_sigma;
    dir = f_dir;
    edge_id = f_edge_id;
    v0 = f_v0;
    v1 = f_v1;
    v0id = f_v0id;
    v1id = f_v1id;
    s = nullptr;
}

int Window::get_edge_id()
{
    return edge_id;
}

int Window::get_dir()
{
    return dir;
}

double Window::get_sigma()
{
    return sigma;
}

double Window::get_b0()
{
    return b0;
}

double Window::get_b1()
{
    return b1;
}

void Window::set_b0(double new_b0)
{
    b0 = new_b0;
}

void Window::set_b1(double new_b1)
{
    b1 = new_b1;
}

Vector3d Window::get_v0()
{
    return v0;
}

Vector3d Window::get_v1()
{
    return v1;
}

int Window::get_v0id()
{
    return v0id;
}

int Window::get_v1id()
{
    return v1id;
}

double Window::get_d0()
{
    return d0;
}

void Window::set_d0(double f_d0)
{
    d0 = f_d0;
}

void Window::set_d1(double f_d1)
{
    d1 = f_d1;
}

double Window::get_d1()
{
    return d1;
}

Vector2d Window::get_s()
{
    if (s == nullptr)
    {
        double x0, x1, c, delta, x, y_s1, y_s2;
        x0 = b0;
        x1 = b1;
        x = (d1 * d1 - d0 * d0 - x1 * x1 + x0 * x0) / (2 * (x0 - x1));
        c = x1 * x1 + x * x - 2 * x1 * x - d1 * d1;
        delta = sqrt(-4 * c);
        y_s1 = delta / 2;
        y_s2 = -delta / 2;

        // HERE VS SHOUDL BE CALLED S AND S SHOULD BE STORE IN WINDOW TO ALLOW FOR WINDOWS INTERSECTION.
        // we have the two possibles sources, we choosed the one with a positive y
        if (y_s1 > 0.)
        {
            s = new Vector2d(x, y_s1);
        }
        else
        {
            s = new Vector2d(x, y_s2);
        }
    }
    return *s;
}

void Window::print()
{
    cout
        << "W(b0=" << b0
        << ", b1=" << b1
        << ", d0=" << d0
        << ", d1=" << d1
        << ", sig=" << sigma
        << ", eid=" << edge_id
        << ", v0=[" << v0(0) << "," << v0(1) << "," << v0(2) << "]"
        << ", v1=[" << v1(0) << "," << v1(1) << "," << v1(2) << "]"
        << ")";
}

double Window::min_geodist()
{
    return 0.0;
}

double Window::max_geodist()
{
    return 0.0;
}