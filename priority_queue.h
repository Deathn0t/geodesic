#include <iostream>
#include <ostream>
#include <tuple>
#include <list>
#include <iterator>

#include "window.h"

using namespace std;

class PriorityQueue
{
    list<Window *> l;

public:
    void push(Window *v);

    void print();

    void remove(Window *v);

    void pop();

    Window *top();

    bool empty();
};