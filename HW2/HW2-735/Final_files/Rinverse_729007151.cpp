#include <iostream>
#include <cstdlib>
#include <cmath> 
#include <vector>
#include <omp.h>
#include <chrono>

using namespace std;


void display( vector<vector<double> > &A) 
{
    //cout << "Display_func";  
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
vector<vector<double> > mat_copy( vector<vector<double> > &A)
{
    //cout << "mat_copy";
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

    //cout << "create_random_matrix";
    int n = A.size();
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            A[i][j] = (double) rand();
        }
    }

    return;
}

void absol(vector<vector<double> > &A)
{  
    //cout << "absol";
    int n = A.size();
    double x;
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
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
    int n = A.size();
    for(int i = 0; i < n; i++)
    {
        temp_sum = 0.0;
        for(int j = 0; j < n; j++)
        {
            temp_sum += A[j][i]; 
        }
        A_arr[i] = temp_sum;
    }

    return; 
}

void diag(vector<vector<double> > &A, vector<double> &B)
{
    int n = A.size();
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
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
    int n = A.size();
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
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
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
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
//Gauss Jordan inverse algorithm
vector<vector<double> > inverse(vector<vector<double> > &A)
{

    int n = A.size();
    vector<vector<double> > matrix (n, vector<double>(2 * n,0));
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            matrix[i][j] = A[i][j];
        }
    }
    //Augmented matrix is created
    double temp; 
    for (int i = 0; i < n; i++) { 

        for (int j = 0; j < 2 * n; j++) { 

            if (j == (i + n)) 
                matrix[i][j] = 1; 
            else if (j >= n)
                matrix[i][j] = 0;

        } 
    } 

    //Row exchange
    for (int i = n - 1; i > 0; i--) 
    { 

        if (matrix[i - 1][0] < matrix[i][0]) 
        { 
            vector<double> temp = matrix[i]; 
            matrix[i] = matrix[i - 1]; 
            matrix[i - 1] = temp; 
        } 
    } 
    
    //Subtract and multiply with a constant
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
    //multiply with a constant
    for (int i = 0; i < n; i++) { 

        temp = matrix[i][i]; 
        for (int j = 0; j < 2 * n; j++) { 

            matrix[i][j] = (float) matrix[i][j] / temp; 
        } 
    } 

    //copy from the augmented matrix to matrix A;
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

    for(int i=0; i<row; i++)
    {
        for(int j=0; j<col; j++)
        {
            A[i][j] = (-1.0) * A[i][j];
        }
    }
    return;

}

//recursive function for inverse computation

vector<vector<double> >  compute_inverse(vector<vector<double> > &A, int start_row, int end_row, int start_col, int end_col)
{
    int n = end_row-start_row + 1;


    vector<vector<double> > A_inverse(n , std::vector<double>(n , 0));

    if (n < 16) {

        vector<vector<double> > A_temp(n, std::vector<double>(n, 0));
        for (int i=0; i<n; i++){

            for (int j=0; j<n; j++){

                A_temp[i][j] = A[start_row+i][start_col+j];

            }
        }
        A_inverse = inverse(A_temp);

    }
    else{

        int mat1_dim = n/2;
        int mat2_dim = n/2;

        if (n % 2 !=0)
        {
            mat2_dim = (n/2)+1;
        }

        int mat3_row = mat1_dim;
        int mat3_col = mat2_dim;

        //create matrices for storage
        vector<vector<double> > Mat1(mat1_dim, std::vector<double>(mat1_dim, 0));
        vector<vector<double> > Mat2(mat2_dim, std::vector<double>(mat2_dim, 0));
        vector<vector<double> > Mat3_temp(mat3_row, std::vector<double>(mat3_col, 0));
        vector<vector<double> > Mat3(mat3_row, std::vector<double>(mat3_col, 0));
        //Recursively call for inverse computation
        #pragma omp task shared(A, Mat1)

        Mat1 = compute_inverse(A, start_row, start_row+(mat1_dim-1), start_col, start_col+(mat1_dim-1));
        
        #pragma omp task shared(A, Mat2)

        Mat2 = compute_inverse(A, start_row+mat1_dim, end_row, start_col+mat1_dim, end_col);
	#pragma omp taskwait

        for(int i=0; i<mat3_row; i++){
            for(int j=0; j<mat3_col; j++){
                Mat3_temp[i][j] = A[start_row+i][start_col+mat1_dim+j];
            }
        }
        //Computing inverse for matrix 3 using the given formula in matlab code given

        Mat3 = matrix_mul(Mat1, Mat3_temp);
        Mat3 = matrix_mul(Mat3, Mat2);
        mat_negate(Mat3);

        //Write them to corresponding locations
	#pragma omp parallel for
        for (int i=0; i<mat1_dim; i++){
            for (int j=0; j<mat1_dim; j++){
                A_inverse[i][j] = Mat1[i][j];
            }
        }
	#pragma omp parallel for
        for (int i=0; i<mat2_dim; i++){
            for (int j=0; j<mat2_dim; j++){
                A_inverse[mat1_dim+i][mat1_dim+j] = Mat2[i][j];
            }
        }

	#pragma omp parallel for
        for (int i=0; i<mat3_row; i++){
            for (int j=0; j<mat3_col; j++){
                A_inverse[i][mat1_dim+j] = Mat3[i][j];
            }
        }

    }
    return A_inverse;
}

int main(int argc, char *argv[])
{
    int n = atoi(argv[1]);

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
    //display(F);

    if(n == 0)
    {
        cout << "Dimension should be greater than 0" << endl;
        exit(0);
    }
    else if (n == 1)
    {
        A[0][0] = 1/A[0][0];
    }
    auto start_time = std::chrono::high_resolution_clock::now();

    #pragma omp parallel
    #pragma omp single
 
    A = compute_inverse(F, 0, n - 1, 0, n - 1);

    auto stop_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> dur_ms = stop_time - start_time;
    std::cout << "Time elapsed: " << dur_ms.count() << "ms" << std::endl;

    //display(A);
    return 0;
}
