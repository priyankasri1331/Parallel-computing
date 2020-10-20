#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>

using namespace std;

class create_triu
{
    
    public:

    double **A;

    /*create_triu( int n)
    {
        A = new double*[n];
        for(int i = 0; i < n; i++)
            A[i] = new double[n];
    }*/

    double** create_random_matrix(int n)
    {
        A = new double*[n];
        for(int i = 0; i < n; i++)
            A[i] = new double[n];

        for(int i = 0; i < n; i++)
        {

            for(int j = 0; j < n; j++)
            {
                A[i][j] = (double) rand() / RAND_MAX;
            }
        }
        return (double**) A;
    }

    /*double** abs(double **A)
    {  
        
        return (double**) A;

    }
    double** sum(double **A)
    {

        return (double**) A;
    }

    double** diag(double **A)
    {


    }

    double** add(double **A)
    {


    }

    double** triu(double **A)
    {


    }*/
};

int main()
{
    int n = 2;
    double **A;
    int i;
    int j;
    
    create_triu obj1;
    
    A = obj1.create_random_matrix(n);

    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            cout << A[i][j] << endl;
        }
    }
    return 0;
}