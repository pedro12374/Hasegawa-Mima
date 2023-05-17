#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int exam1[4];// array that can hold 100 numbers for 1st column 
int exam2[4];// array that can hold 100 numbers for 2nd column 
int exam3[4];// array that can hold 100 numbers for 3rd column  
int t[4];
int main() // int main NOT void main
{ 
  ifstream infile;   

  int num = 0; // num must start at 0
  infile.open("test.dat");// file containing numbers in 3 columns 
     if(infile.fail()) // checks to see if file opended 
    { 
      cout << "error" << endl; 
      return 1; // no point continuing if the file didn't open...
    } 
       while(!infile.eof()) // reads file to end of *file*, not line
      {  
          infile >>t[num];
         infile >> exam1[num]; // read first column number
         infile >> exam2[num]; // read second column number
         infile >> exam3[num]; // read third column number

         ++num; // go to the next number

         // you can also do it on the same line like this:
         // infile >> exam1[num] >> exam2[num] >> exam3[num]; ++num;
      } 
  infile.close(); 
  for(int i=0;i<4;i++){
    cout<<t[i]<<"\t"<<exam1[i]<<"\t"<<exam2[i]<<"\t"<<exam3[i]<<endl;
  }
  return 0; // everything went right.
} 
