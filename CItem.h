#pragma once
#include <vector>
using namespace std;

class CItem
{
private:
    int num;
    int weight;
    bool loaded;

public:
    CItem();
    CItem(int _num, int _weight);
    int getNum();
    int getWeight();
    void setLoaded(bool b);
    bool getLoaded();
    void printItem();
};