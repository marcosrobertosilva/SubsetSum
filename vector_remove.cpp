#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;


void print(const vector<int> & v)
{
	for(vector<int>::const_iterator it = v.begin(); it != v.end(); it++)
		cout << *it << " ";
	cout << endl;
}

/* Remove a *primeira* ocorrencia de val no vector */
void remove_val(vector<int> & v, int val)
{
	std::vector<int>::iterator position = std::find(v.begin(), v.end(), val);
	if(position != v.end()) // == myVector.end() means the element was not found
    	v.erase(position);
}

int main()
{
	vector<int> v;
	v.push_back(0);
	v.push_back(1);
	v.push_back(2);
	v.push_back(3);
	v.push_back(3);
	v.push_back(4);
	v.push_back(5);
	v.push_back(3);
	v.push_back(6);


	print(v);

	cout << "Apos a remocao do valor 3...\n";
	remove_val(v, 3);
	print(v);

	return 0;
}