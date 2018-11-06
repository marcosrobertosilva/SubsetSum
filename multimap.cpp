#include <iostream> 
#include <map> 
#include <iterator> 
  
using namespace std; 
  
int main() 
{ 
    multimap <int, int> mymap;        // empty multimap container 
  
    // insert elements in random order 
    mymap.insert(pair <int, int> (1, 40)); 
    mymap.insert(pair <int, int> (2, 30)); 
    mymap.insert(pair <int, int> (2, 35)); 
    mymap.insert(pair <int, int> (2, 37)); 
    mymap.insert(pair <int, int> (3, 60)); 
    mymap.insert(pair <int, int> (4, 20)); 
    mymap.insert(pair <int, int> (5, 50)); 
    mymap.insert(pair <int, int> (6, 50)); 
    mymap.insert(pair <int, int> (6, 10)); 
  
    // printing multimap gquiz1 
    multimap <int, int> :: iterator itr; 
    cout << "\nThe multimap mymap is : \n"; 
    cout << "\tKEY\tELEMENT\n"; 
    for (itr = mymap.begin(); itr != mymap.end(); ++itr) 
    { 
        cout  <<  '\t' << itr->first 
              <<  '\t' << itr->second << '\n'; 
    } 
    cout << endl; 
  
    mymap.erase(mymap.lower_bound(2));
    
    cout << "\nThe multimap mymap after erase is : \n"; 
    cout << "\tKEY\tELEMENT\n"; 
    for (itr = mymap.begin(); itr != mymap.end(); ++itr) 
    { 
        cout  <<  '\t' << itr->first 
              <<  '\t' << itr->second << '\n'; 
    } 
    cout << endl; 
    
    return 0; 
  
} 
