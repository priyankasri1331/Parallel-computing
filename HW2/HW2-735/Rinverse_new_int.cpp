#include <iostream>
#include <cstdlib>
#include <cmath> 
#include <vector>

#define N 6

using namespace std;


void display( vector<vector<double> > &A) 
{ 
    for (int i=0; i<N; i++) 
    { 
        for (int j=0; j<N; j++) 
            cout << A[i][j] << " "; 
        cout << endl; 
    } 
} 
void display2( vector<vector<double> > &A)
{

    for (int i=0; i<N; i++) 
    { 
        for (int j=0; j<N*2; j++) 
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


class create_triu
{
    public:
    int n;
    create_triu(int _n)
    {
        n = _n;
    }

    void create_random_matrix(vector<vector<double> > &A)
    {


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

        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                A[i][j] = A[i][j] + B[i][j];
            }
        }
        return; 
    }

    void triu(vector<vector<double> > &A)
    {

        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                if((i == j) || (i < j))
                {
                    A[i][j] = A[i][j];
                }
                else
                {
                    A[i][j] = 0.0;
                }
                
            }
        }
        return;
    }

    vector<vector<double> > inverse(vector<vector<double> > &A)
    {

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

    vector<vector<double> > negate( vector<vector<double> > &A)
    {
        int row = A.size();
        int col = A[0].size();
        //if(row_dimension!=col_dimension){
        //    throw std::runtime_error(std::string("Calc Inverse Error: Size > 2\n"));
        //}
        vector<vector<double> > temp(row, std::vector<double>(col, 0));

        for(int i=0; i<row; i++)
        {
            for(int j=0; j<col; j++)
            {
                temp[i][j] = (-1.0) * A[i][j];
            }
        }
        return temp;
    }
    vector<vector<double> >  compute_inverse( vector<vector<double> > &A, int start_row, int end_row, int start_col, int end_col)
    {
        int n = end_row - start_row + 1;
        vector<vector<double> > A_copy (n, vector<double>(n,0));
        vector<vector<double> > A_inv (n, vector<double>(n,0));

        if (n<=6)
        {
            for (int i=0; i<n; i++)
            {
                for (int j=0; j<n; j++)
                {
                    A_copy[i][j] = A[start_row+i][start_col+j];
                }
            }   
            cout<<"Test 1\n";
            inverse(A_copy);
            cout<<"Test 2\n";
        }
        else
        {
            //int n1;// = round(n/2);
            int mat1_dim = n / 2;
            int mat2_dim = n / 2;
            int mat3_row = n / 2;
            int mat3_col = n / 2;
            if(n % 2 != 0)
            {
                mat2_dim = n / 2 + 1;
                mat3_row = n / 2;
                mat3_col = n / 2 + 1;
            }

            //cout << "Reached here";

            //static int count = 0;
            vector<vector<double> > Mat_1 (mat1_dim, vector<double>(mat1_dim,0));
            vector<vector<double> > Mat_2 (mat2_dim, vector<double>(mat2_dim,0));
            vector<vector<double> > Mat_3 (mat3_row, vector<double>(mat3_col,0));
            vector<vector<double> > Mat_temp (mat3_row, vector<double>(mat3_col,0));
            //cout << count++ << endl;

            Mat_1 = compute_inverse(A, start_row , start_row + (mat1_dim-1), start_col, start_col + (mat1_dim-1));
            Mat_2 = compute_inverse(A, start_row + mat2_dim, end_row, start_col+mat1_dim, end_col);
            //display(Mat_1);
            cout << endl;
            //display(Mat_2);

            for(int i=0; i < mat3_row; i++)
            {
                for(int j=0; j < mat3_col; j++)
                {
                    Mat_temp[i][j] = A[start_row+i][start_col+mat1_dim+j];
                }
            }

            Mat_3 = matrix_mul(Mat_1, Mat_temp);
            Mat_3 = matrix_mul(Mat_3, Mat_2);
            Mat_3 = negate(Mat_3);

            for (int i=0; i<mat1_dim; i++){
                for (int j=0; j<mat1_dim; j++){
                    A_copy[i][j] = Mat_1[i][j];
                }
            }

            for (int i=0; i<mat2_dim; i++){
                for (int j=0; j<mat2_dim; j++){
                    A_copy[mat1_dim+i][mat1_dim+j] = Mat_2[i][j];
                }
            }

            for (int i=0; i<mat3_row; i++){
                for (int j=0; j<mat3_col; j++){
                    A_copy[i][mat1_dim+j] = Mat_3[i][j];
                }
            }

            //std::cout<<"A_inv:\n";
            //print_matrix(A_inv);

        }
        display(A_copy);
        return A_copy;
    }
    
};



int main()
{
    int n = N;
    int n1 = 2 * N;

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

    create_triu obj1(n);
    obj1.create_random_matrix(A);
    obj1.absol(A);
    obj1.sum(A, B);
    obj1.diag(C, B);
    obj1.add(A,C);
    obj1.triu(A);
    //display(A);
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
    //    obj1.inverse(A);
    //}
    //else
    //{
    A = obj1.compute_inverse(A, 0, n, 0, n);
    //}
    //display(A);
    return 0;
}