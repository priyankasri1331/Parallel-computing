#include <chrono>
#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <cmath> 
#include <vector>
#include <math.h>

using namespace std;


void display( vector<vector<double> > &A) 
{ 
    for (int i=0; i<A.size(); i++) 
    { 
        for (int j=0; j<A[0].size(); j++)
        {
            //cout << "I'm here";
            cout << A[i][j] << " "; 

        } 

        cout << endl; 
    } 
} 
void display2( vector<vector<double> > &A)
{

    for (int i=0; i<A.size(); i++) 
    { 
        for (int j=0; j<A[0].size()*2; j++) 
            cout << A[i][j] << " "; 
        cout << endl; 
    } 
}
vector<vector<double> > mat_copy( vector<vector<double> > &A)
{
    vector<vector<double> > temp (A.size(), vector<double>(A[0].size(),0));
    for (int i=0; i<A.size(); i++) 
    { 
        for (int j=0; j<A[0].size(); j++) 
        {
            temp[i][j] = A[i][j];
        }
    }

    return temp;
}
void create_random_matrix(vector<vector<double> > &A)
{


    for(int i = 0; i < A.size(); i++)
    {
        for(int j = 0; j < A[0].size(); j++)
        {
            A[i][j] = (double) rand();
        }
    }

    return;
}

void absol(vector<vector<double> > &A)
{  
    double x;
    for(int i = 0; i < A.size(); i++)
    {
        for(int j = 0; j < A[0].size(); j++)
        {
            x = A[i][j];
            A[i][j] = abs(x);
        }
    }
    return;
}

void sum( vector<vector<double> > &A, vector<double> &A_arr)
{

    double temp_sum = 0.0;
    
    for(int i = 0; i < A.size(); i++)
    {
        temp_sum = 0.0;
        for(int j = 0; j < A[0].size(); j++)
        {
            temp_sum += A[j][i]; 
        }
        A_arr[i] = temp_sum;
    }

    return; 
}

void diag(vector<vector<double> > &A, vector<double> &B)
{

    for(int i = 0; i < A.size(); i++)
    {
        for(int j = 0; j < A[0].size(); j++)
        {
            if(i == j)
            {
                A[i][j] = B[i];
            }
            else
            {
                A[i][j] = 0.0;
            }
            
        }
    }
    return; 
}

void add(vector<vector<double> > &A, vector<vector<double> > &B)
{

    for(int i = 0; i < A.size(); i++)
    {
        for(int j = 0; j < A[0].size(); j++)
        {
            A[i][j] = A[i][j] + B[i][j];
        }
    }
    return; 
}

vector<vector<double> > triu(vector<vector<double> > &A)
{
    int n = A.size();
    vector<vector<double> > matrix (n, vector<double>(n,0));
    for(int i = 0; i < A.size(); i++)
    {
        for(int j = 0; j < A.size(); j++)
        {
            if((i == j) || (i < j))
            {
                matrix[i][j] = A[i][j];
            }
            else
            {
                matrix[i][j] = 0.0;
            }
            
        }
    }
    return matrix;
}

vector<vector<double> > inverse(vector<vector<double> > &A)
{
    int n = A.size();
    vector<vector<double> > matrix (A.size(), vector<double>(2 * A.size(),0));
    for(int i = 0; i < A.size(); i++)
    {
        for(int j = 0; j < A.size(); j++)
        {
            matrix[i][j] = A[i][j];
        }
    }

    double temp; 
    #pragma omp parallel for
    for(int i = 0; i < A.size(); i++)
    {
        for(int j = 0; j < A.size(); j++)
        {

            if (j == (i + n)) 
                matrix[i][j] = 1; 
            else if (j >= n)
                matrix[i][j] = 0;

        } 
    } 

    #pragma omp parallel for
    for (int i = n - 1; i > 0; i--) 
    { 

        if (matrix[i - 1][0] < matrix[i][0]) 
        { 
            vector<double> temp = matrix[i]; 
            matrix[i] = matrix[i - 1]; 
            matrix[i - 1] = temp; 
        } 
    } 
    
    #pragma omp parallel for
    for (int i = 0; i < n; i++) { 

        for (int j = 0; j < n; j++) { 

            if (j != i) { 

                temp = (float) matrix[j][i] / matrix[i][i]; 
                for (int k = 0; k < 2 * n; k++) { 

                    matrix[j][k] -= matrix[i][k] * temp; 
                } 
            } 
        } 
    } 
    #pragma omp parallel for

    for (int i = 0; i < n; i++) { 

        temp = matrix[i][i]; 
        for (int j = 0; j < 2 * n; j++) { 

            matrix[i][j] = (float) matrix[i][j] / temp; 
        } 
    } 

    #pragma omp parallel for
    for(int i = 0; i < n ; i++)
    {
        for(int j = n; j < n *2; j++)
        {
            A[i][j - n] = matrix[i][j];
        }
    }

    return A; 
}

vector<vector<double> > matrix_mul( vector<vector<double> > &A, vector<vector<double> > &B)
{
    int A_row = A.size();
    int A_col = A[0].size();
    int B_row = B.size();
    int B_col = B[0].size();

    vector<vector<double> > temp(A_row, std::vector<double>(B_col, 0));

    for(int i = 0; i < A_row; i++)
    {
        for(int j = 0; j < B_col; j++)
        {
            for(int k = 0; k < A_col; k++)
            {
                temp[i][j] += A[i][k]*B[k][j];
            }
        }
    }
    return temp;
}
void mat_negate(vector<vector<double> > &A)
{
    int row = A.size();
    int col = A[0].size();
    //if(row_dimension!=col_dimension){
    //    throw std::runtime_error(std::string("Calc Inverse Error: Size > 2\n"));
    //}
    //vector<vector<double> > temp(row, std::vector<double>(col, 0));

    for(int i=0; i<row; i++)
    {
        for(int j=0; j<col; j++)
        {
            A[i][j] = (-1.0) * A[i][j];
        }
    }
    return;

}

int main(int argc, char *argv[])
{
    int n = std::atoi(argv[1]);

    if(n == 0)
    {
        cout << "n should be greater than 0";
        exit(0);
    }
    vector<vector<double> > A(n, vector<double>(n,0));
    vector<double> B(n, 0);
    vector<vector<double> > C(n, vector<double>(n,0));
    vector<vector<double> > D(n, vector<double>(n,0));
    vector<vector<double> > E(n, vector<double>(n,0));
    vector<vector<double> > F(n, vector<double>(n,0));
    
    double G;

    create_random_matrix(A);
    absol(A);
    sum(A, B);
    diag(C, B);
    add(A,C);
    F = triu(A);
    display(F);
    //A = obj1.matrix_mul(A,A);
    //A = obj1.negate(A);
    //display(A);

    if(n == 0)
    {
        cout << "Dimension should be greater than 0" << endl;
        exit(0);
    }
    else if (n == 1)
    {
        A[0][0] = 1/A[0][0];
    }

    auto start_time = chrono::high_resolution_clock::now();

    //omp_set_nested(1);
    #pragma omp parallel
    #pragma omp single
    A = inverse(F);
    display(A);

    //std::cout<<"Inverse Matrix: \n";
    //print_matrix(A_inv);
    auto stop_time = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> dur_ms = stop_time - start_time;
    cout << "Time elapsed: " << dur_ms.count() << "ms" << std::endl;

    return 0;
}