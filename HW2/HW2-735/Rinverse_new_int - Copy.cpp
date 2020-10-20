#include <iostream>
#include <cstdlib>
#include <cmath> 
#include <vector>

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
        for (int j=0; j<A.size()*2; j++) 
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

    double temp; 
    for (int i = 0; i < n; i++) { 

        for (int j = 0; j < 2 * n; j++) { 

            if (j == (i + n)) 
                matrix[i][j] = 1; 
            else if (j >= n)
                matrix[i][j] = 0;

        } 
    } 


    for (int i = n - 1; i > 0; i--) 
    { 

        if (matrix[i - 1][0] < matrix[i][0]) 
        { 
            vector<double> temp = matrix[i]; 
            matrix[i] = matrix[i - 1]; 
            matrix[i - 1] = temp; 
        } 
    } 
    

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

    for (int i = 0; i < n; i++) { 

        temp = matrix[i][i]; 
        for (int j = 0; j < 2 * n; j++) { 

            matrix[i][j] = (float) matrix[i][j] / temp; 
        } 
    } 


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

vector<vector<double> >  compute_inverse(vector<vector<double> > &A, int start_row, int end_row, int start_col, int end_col)
{
    int dimension = end_row-start_row + 1;
    //cout << dimension;

    vector<vector<double> > A_inv(dimension , std::vector<double>(dimension , 0));
    //cout << A_inv.size();
    if (dimension <= 8) {
        //cout << "--------";
        vector<vector<double> > A_copy(dimension, std::vector<double>(dimension, 0));
        for (int i=0; i<dimension; i++){
            //cout << i << " ------------------\n";
            for (int j=0; j<dimension; j++){
                //cout << "----------------\n";
                A_copy[i][j] = A[start_row+i][start_col+j];
                //cout << j << "||||||||||||||||||\n";
            }
        }
        //cout<<"Test 1\n";
        A_inv = inverse(A_copy);
        //cout<<"Test 2\n";
    }
    else{
        //cout << "Else --" << endl;
        //dimension can be even or odd
        //For odd cases we need to do handling slight differently
        int top_dim;
        int btm_dim;

        if (dimension%2==0){
            top_dim = dimension/2;
            btm_dim = dimension/2;
        }
        else{
            top_dim = dimension/2;
            btm_dim = (dimension/2)+1;
        }
        int corner_row_dim = top_dim;
        int corner_col_dim = btm_dim;

        //Creating Empty Matrices for calculating inverse of sub-problems
        vector<vector<double> > A_inv_top(top_dim, std::vector<double>(top_dim, 0));
        vector<vector<double> > A_inv_btm(btm_dim, std::vector<double>(btm_dim, 0));
        vector<vector<double> > A_corner(corner_row_dim, std::vector<double>(corner_col_dim, 0));
        vector<vector<double> > A_inv_corner(corner_row_dim, std::vector<double>(corner_col_dim, 0));

        A_inv_top = compute_inverse(A, start_row, start_row+(top_dim-1), start_col, start_col+(top_dim-1));

        A_inv_btm = compute_inverse(A, start_row+top_dim, end_row, start_col+top_dim, end_col);


        for(int i=0; i<corner_row_dim; i++){
            for(int j=0; j<corner_col_dim; j++){
                A_corner[i][j] = A[start_row+i][start_col+top_dim+j];
            }
        }

        A_inv_corner = matrix_mul(A_inv_top, A_corner);
        A_inv_corner = matrix_mul(A_inv_corner, A_inv_btm);
        mat_negate(A_inv_corner);

        for (int i=0; i<top_dim; i++){
            for (int j=0; j<top_dim; j++){
                A_inv[i][j] = A_inv_top[i][j];
            }
        }

        for (int i=0; i<btm_dim; i++){
            for (int j=0; j<btm_dim; j++){
                A_inv[top_dim+i][top_dim+j] = A_inv_btm[i][j];
            }
        }

        for (int i=0; i<corner_row_dim; i++){
            for (int j=0; j<corner_col_dim; j++){
                A_inv[i][top_dim+j] = A_inv_corner[i][j];
            }
        }

        //std::cout<<"A_inv:\n";
        //display(A_inv);

    }
    return A_inv;
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
    //else if(n < 16)
    //{
        //A = inverse(F);
        //display(A);
    //}
    //else
    //{
    A = compute_inverse(F, 0, n - 1, 0, n - 1);
    //}
    display(A);
    return 0;
}