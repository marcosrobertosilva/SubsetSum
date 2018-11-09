#include "CItem.h"
#include <iostream>
using namespace std;

CItem::CItem()
{

}

CItem::CItem(int _num, int _weight)
{
    num = _num;
    weight = _weight;
    loaded = false;
}

int CItem::getNum()
{
    return num;
}

int CItem::getWeight()
{
    return weight;
}

bool CItem::getLoaded()
{
    return loaded;
}

void CItem::setLoaded(bool b)
{
    loaded = b;
}

void CItem::printItem()
{
    cout << "Item no. = " << num << " peso = " << weight << " inserido em algum bin? " << loaded << endl;
}