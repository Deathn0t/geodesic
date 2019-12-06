#include "window.h"

Window::Window() {}

Window::Window(double f_b0, double f_b1, double f_d0, double f_d1, double f_sigma, int f_dir, int f_edge_id, Vector3d &f_v0, Vector3d &f_v1, int f_v0id, int f_v1id)
{
    if (f_b0 < 1e-7)
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
    if (new_b0 < 0)
    {
        cout << "new_b0: " << new_b0 << endl;
        // exit(0);
    }
    b0 = new_b0;
}

void Window::set_b1(double new_b1)
{
    if (new_b1 > (v0 - v1).norm())
    {
        cout << "new_b1: " << new_b1 << ", norm: " << (v0 - v1).norm() << endl;
        // exit(0);
    }
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
        return Vector2d(x, y_s1);
    }
    else
    {
        return Vector2d(x, y_s2);
    }
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
    double min_dist_w_s = 0.0;

    Vector2d w_b0_2d = Vector2d(get_b0(), 0);
    Vector2d w_b1_2d = Vector2d(get_b1(), 0);

    Vector2d s_w = get_s();
    double dist_w_s_b0 = (s_w - w_b0_2d).norm() + get_sigma();
    double dist_w_s_b1 = (s_w - w_b1_2d).norm() + get_sigma();

    // w

    if (s_w[0] >= get_b0() && s_w[0] <= get_b1())
    {
        Vector2d point_on_segment = Vector2d(s_w[0], 0);

        min_dist_w_s = (s_w - point_on_segment).norm() + get_sigma();
    }
    else if (s_w[0] < get_b0())
    {
        Vector2d b0_2d = Vector2d(get_b0(), 0);
        min_dist_w_s = (s_w - b0_2d).norm() + get_sigma();
    }
    else
    {
        Vector2d b1_2d = Vector2d(get_b1(), 0);
        min_dist_w_s = (s_w - b1_2d).norm() + get_sigma();
    }
    return min_dist_w_s;
}

double Window::max_geodist()
{
    double max_dist_w_s = 0.0;

    Vector2d w_b0_2d = Vector2d(get_b0(), 0);
    Vector2d w_b1_2d = Vector2d(get_b1(), 0);

    Vector2d s_w = get_s();
    double dist_w_s_b0 = (s_w - w_b0_2d).norm() + get_sigma();
    double dist_w_s_b1 = (s_w - w_b1_2d).norm() + get_sigma();

    // w

    if (s_w[0] >= get_b0() && s_w[0] <= get_b1())
    {
        Vector2d point_on_segment = Vector2d(s_w[0], 0);

        max_dist_w_s = max(dist_w_s_b0, dist_w_s_b1);
    }
    else if (s_w[0] < get_b0())
    {
        Vector2d b0_2d = Vector2d(get_b0(), 0);
        max_dist_w_s = dist_w_s_b1;
    }
    else
    {
        Vector2d b1_2d = Vector2d(get_b1(), 0);
        max_dist_w_s = dist_w_s_b0;
    }
    return max_dist_w_s;
}