#include <iostream>
#include <cstdlib>
#include <cmath> 

#define N 3

using namespace std;


void display(double **A) 
{ 
    for (int i=0; i<N; i++) 
    { 
        for (int j=0; j<N; j++) 
            cout << A[i][j] << " "; 
        cout << endl; 
    } 
} 
void display2(double **A)
{
    for (int i=0; i<N; i++) 
    { 
        for (int j=0; j<N*2; j++) 
            cout << A[i][j] << " "; 
        cout << endl; 
    } 
}


class create_triu
{
    public:
    int n;
    //double **A;
    create_triu(int _n)
    {
        n = _n;
    }

    double** create_random_matrix()
    {
        double **A = new double*[n];
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

    double** absol(double **A)
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
        return (double**) A;
    }

    double* sum(double **A)
    {
        double *A_arr = new double[n];
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
        return (double*) A_arr;
    }

    double** diag(double *B)
    {
        double **A_diag = new double *[n];
        for(int i = 0; i < n; i++)
            A_diag[i] = new double[n];
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                if(i == j)
                {
                    A_diag[i][j] = B[i];
                }
                else
                {
                    A_diag[i][j] = 0.0;
                }
                
            }
        }
        return (double**) A_diag;
    }

    double** add(double **A, double **B)
    {
        double **A_add = new double *[n];
        for(int i = 0; i < n; i++)
            A_add[i] = new double[n];
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                A_add[i][j] = A[i][j] + B[i][j];
            }
        }
        return (double**) A_add;
    }

    double** triu(double **A)
    {
        double **A_triu = new double *[n];
        for(int i = 0; i < n; i++)
            A_triu[i] = new double[n];
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                if((i == j) || (i < j))
                {
                    A_triu[i][j] = A[i][j];
                }
                else
                {
                    A_triu[i][j] = 0.0;
                }
                
            }
        }
        return (double**) A_triu;
    }

    void inverse(double **A)
    {
        double **matrix = new double *[n];
        for(int i = 0; i < n * 2; i++)
            matrix[i] = new double[n];

        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                matrix[i][j] = A[i][j];
            }
        }

        double temp; 
        printf("=== Matrix ===\n");
        display(matrix);


    
        // Create the augmented matrix 
        // Add the identity matrix 
        // of order at the end of original matrix. 
        for (int i = 0; i < n; i++) { 
    
            for (int j = 0; j < 2 * n; j++) { 
    
                // Add '1' at the diagonal places of 
                // the matrix to create a identity matirx 
                if (j == (i + n)) 
                    matrix[i][j] = 1; 
            } 
        } 
        display2(matrix);

        double ratio = 0;

		 for(int i=0;i< n;i++)
		 {
			  if(matrix[i][i] == 0.0)
			  {
				   printf("Mathematical Error!");
				   exit(0);
			  }
			  for(int j=1;j<n;j++)
			  {
				   if(i!=j)
				   {
					    ratio = matrix[j][i]/matrix[i][i];
					    for(int k=0;k < 2*n;k++)
					    {
					     	matrix[j][k] = matrix[j][k] - ratio*matrix[i][k];
					    }
				   }
			  }
		 }
		 /* Row Operation to Make Principal Diagonal to 1 */
		 for(int i=0;i< n;i++)
		 {
			  for(int j=n;j<2*n;j++)
			  {
			   	matrix[i][j] = matrix[i][j]/matrix[i][i];
			  }
		 }

        /*// Interchange the row of matrix, 
        // interchanging of row will start from the last row 
        for (int i = n - 1; i > 0; i--) { 
    
            // Swapping each and every element of the two rows 
            // if (matrix[i - 1][0] < matrix[i][0]) 
            // for (int j = 0; j < 2 * order; j++) { 
            // 
            //        // Swapping of the row, if above 
            //        // condition satisfied. 
            // temp = matrix[i][j]; 
            // matrix[i][j] = matrix[i - 1][j]; 
            // matrix[i - 1][j] = temp; 
            //    } 
    
            // Directly swapping the rows using pointers saves time 
    
            if (matrix[i - 1][0] < matrix[i][0]) 
            { 
                double* temp = matrix[i]; 
                matrix[i] = matrix[i - 1]; 
                matrix[i - 1] = temp; 
            } 
        } 
        display2(matrix);

        // Replace a row by sum of itself and a 
        // constant multiple of another row of the matrix 
        for (int i = 0; i < n; i++) { 
    
            for (int j = 0; j < n; j++) { 
    
                if (j != i) { 
    
                    temp = matrix[j][i] / matrix[i][i]; 
                    for (int k = 0; k < 2 * n; k++) { 
    
                        matrix[j][k] -= matrix[i][k] * temp; 
                    } 
                } 
            } 
        } 

        display2(matrix);
    
        // Multiply each row by a nonzero integer. 
        // Divide row element by the diagonal element 
        for (int i = 0; i < n; i++) { 
    
            temp = matrix[i][i]; 
            for (int j = 0; j < 2 * n; j++) { 
    
                matrix[i][j] = matrix[i][j] / temp; 
            } 
        } 
        */
    
        // print the resultant Inverse matrix. 
        //printf("\n=== Inverse Matrix ===\n"); 
        //PrintInverse(matrix, order, 2 * order); 

        for(int i = 0; i < n ; i++)
        {
            for(int j = n; j < n *2; j++)
            {
                A[i][j] = matrix[i][j];
            }
        }

        display2(matrix);
    
        return; 
    } 

};

/*double** compute_inverse(double **R)
{
    int n = N;
    double **Ri = new double *[n];
    for(int i = 0; i < n; i++)
        Ri[i] = new double[n];

    if (n == 1)
    {
        Ri[0][0] = 1/R[0][0];
    }
    else if(n < 16)
    {
        Ri = inv(R);
    }
    {
        int n1 = round(n/2);
        if (n1 == 2)
        {

        }

        Ri = R;
        Ri();

    }
    return (double**) Ri;
}*/


int main()
{
    int n = N;

    if(n == 0)
    {
        cout << "n should be greater than 0";
        exit(0);
    }
    double **A;
    double *B;
    double **C;
    double **D;
    double **E;
    double **F;
    double G;

    create_triu obj1(n);
    A = obj1.create_random_matrix();
    A = obj1.absol(A);
    B = obj1.sum(A);
    C = obj1.diag(B);
    D = obj1.add(A,C);
    //E = obj1.triu(D);


    //E = obj1.getCofactor(A, 0, 2, 2);



    //display(A);
    obj1.inverse(A); 
    //display(A);
    //E has the triangular matrix

   
    cout << endl;
    //display(E);
    cout << endl;
    //cout << G << endl;
    

    //F = compute_inverse(E);
    return 0;
}
/*
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            cout << A[i][j] << endl;
        }
    }

    for(int i = 0; i < n; i++)
    {
        cout << B[i] << endl;
    }

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
            cout << C[i][j] << endl;
    }
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
            cout << D[i][j] << endl;
    }

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
            cout << E[i][j] << endl;
    }
*/