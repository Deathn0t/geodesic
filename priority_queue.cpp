#include "priority_queue.h"

void PriorityQueue::push(Window *v)
{
    list<Window *> tmp;
    while (!l.empty() && v->min_geodist() > l.front()->min_geodist())
    {
        tmp.push_front(l.front());
        l.pop_front();
    }
    l.push_front(v);
    while (!tmp.empty())
    {
        l.push_front(tmp.front());
        tmp.pop_front();
    }
}

void PriorityQueue::print()
{
    list<Window *>::iterator it;
    for (it = l.begin(); it != l.end(); ++it)
        cout << '\t' << *it;
    cout << '\n';
}

void PriorityQueue::remove(Window *v)
{
    l.remove(v);
}

void PriorityQueue::pop()
{
    l.pop_front();
};

Window *PriorityQueue::top()
{
    return l.front();
}

bool PriorityQueue::empty()
{
    return l.empty();
}
