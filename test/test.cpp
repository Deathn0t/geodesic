#include <iostream>
#include <ostream>
#include <tuple>
#include <list>
#include <iterator>

using namespace std;

class PriorityQueue
{
    list<int> l;

public:
    void push(int v)
    {
        list<int> tmp;
        int c;
        while (!l.empty() && v > l.front())
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

    void print()
    {
        list<int>::iterator it;
        for (it = l.begin(); it != l.end(); ++it)
            cout << '\t' << *it;
        cout << '\n';
    }

    void remove(int v)
    {
        l.remove(v);
    }

    void pop()
    {
        l.pop_front();
    };

    int front()
    {
        return l.front();
    }
};

int main(int argc, char *argv[])
{
    PriorityQueue Q;

    for (int i = 0; i < 10; i++)
    {
        Q.push(i);
    }

    Q.print();
    Q.remove(5);
    Q.print();

    cout << "front: " << Q.front() << endl;
    Q.pop();
    Q.print();

    return 0;
}