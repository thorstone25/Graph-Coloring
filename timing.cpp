// Example of how to compute runtimes in C++

#include <time.h>
#include <iostream>
using namespace std;


int main()
{
   clock_t startTime, endTime;

   int x = 0;
   startTime = clock();
   
   for (int i = 0; i < 1000; i++)
   {
      x++;
      cout << x << endl;
   }

   cout << startTime << endl;
   cout << clock() << endl;

   int diff = clock()-startTime;
   
   cout << (float) diff / CLOCKS_PER_SEC << endl;
}
