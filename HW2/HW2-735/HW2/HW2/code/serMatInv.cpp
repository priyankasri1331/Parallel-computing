//Computing Inverse Sequentially in Code:
//Using Recursion or Gaussian - Elimination Method for Computing the inverse of a matrix
//Creating a random matrix based on given dimension, provided it is leading diagonal are non-zeros.
#include <iomanip>
#include<iostream>
#include <vector>
#include <iterator>
#include <cstdlib>
#include <string>
#include <assert.h>
#include <omp.h>
#include <cmath>
#include <chrono>
using namespace std;
//task0 : Create a class matrix, that has both print and generate function.
//task1 : create a random upper triangular matrix where none of the leading diagonal not 0
//task2 : Create a per class inverse function, for computing the inverse.
class MATRIX
{
private:	
	vector<vector<float>> matrix;
	unsigned long rows;//Make sure that rows and other concepts start from 0, prevent this
	unsigned long columns;//Fucking index bug.
public:
	MATRIX(string,unsigned long,unsigned long);//Constructor, letskeep it dummy, until we get the algorithms and test matrix correct.
	void PrintMatrix();//Displays the output in terminal
	MATRIX Inverse(string) const;//Returns the inverse of matrix.
	MATRIX operator -();//-ve sign operator for R12 computation
	float Element(unsigned long,unsigned long) const;//Returns the element based on row,column index.
	void PushVals(unsigned long,unsigned long,float);//Push values onto matrix at specified location.
	inline unsigned long Rows() const { return this->rows;}
	inline unsigned long Columns() const { return this->columns;}
	//Helper Funcitons
	MATRIX SubMatrix(unsigned long,unsigned long,unsigned long,unsigned long)const;//Returns submatrix of given size
	void AssignSubMatrix(unsigned long,unsigned long,unsigned long,unsigned long,const MATRIX);//Assigns a matrix to given dimensions
};

void MATRIX::AssignSubMatrix(unsigned long rowStartIndex,unsigned long rowEndIndex,unsigned long colStartIndex,unsigned long colEndIndex,const MATRIX inpMatrix)
{
	//Check if the dimensions of Inputmatrix are in accordance with the inputIndices,and \
	check if those indices are in accordance with calling Matrix.
	assert(rowStartIndex<=this->rows && rowEndIndex<=this->rows && colStartIndex<=this->columns && colEndIndex<=this->columns);
	assert(inpMatrix.Rows()== (rowEndIndex-rowStartIndex) && inpMatrix.Columns()== (colEndIndex-colStartIndex));

	unsigned long inpRows = 0;
	while((inpRows<=inpMatrix.Rows()) && (rowStartIndex<=rowEndIndex))
	{
		unsigned long inpCols = 0;
		unsigned long _colStartIndex =colStartIndex;
		while((_colStartIndex<=colEndIndex) && (inpCols<=inpMatrix.Columns()))
		{
			this->matrix[rowStartIndex][_colStartIndex]=inpMatrix.Element(inpRows,inpCols);
			inpCols++;
			_colStartIndex++;
		}
		inpRows++;
		rowStartIndex++;
	}
}



MATRIX MATRIX::SubMatrix(unsigned long rowStartIndex,unsigned long rowEndIndex,unsigned long colStartIndex,unsigned long colEndIndex) const
{
	assert((rowStartIndex<=this->rows) && (rowEndIndex<=this->rows)&&(rowEndIndex>=rowStartIndex));//So that it lies within the range of Actual Matrix
	assert((colStartIndex<=this->columns) && (colEndIndex<=this->columns)&&(colEndIndex>=colStartIndex));
	//Comparison of Unsigned expression is bydefault true of >=0.
	MATRIX outputMatrix("zero",rowEndIndex-rowStartIndex+1,colEndIndex-colStartIndex+1);//Since it is actual number of rows not indices we add 1
	//MATRIX outputMatrix("zero",2,2);
	unsigned long outputRows = 0;
	unsigned long outputCols =0;
	while((rowStartIndex<=rowEndIndex) && (outputRows<=outputMatrix.Rows()))
	{	
		unsigned long _colStartIndex = colStartIndex;
		outputCols=0;
		while((_colStartIndex<=colEndIndex) && (outputCols<=outputMatrix.Columns()))
		{
			//cout<<this->matrix[rowStartIndex][_colStartIndex];
			outputMatrix.PushVals(outputRows,outputCols,this->matrix[rowStartIndex][_colStartIndex]);
			_colStartIndex++;
			outputCols++;
		}
		//cout<<endl;
		outputRows++;
		rowStartIndex++;
	}
return outputMatrix;
}


MATRIX MATRIX::operator - ()
{
	MATRIX outputMatrix("zero",this->rows+1,this->columns+1);
	assert(outputMatrix.Rows()==this->rows);
	assert(outputMatrix.Columns()==this->columns);
	
	for(unsigned long rows = 0;rows<=this->rows;rows++)
	{
		for(unsigned long columns = 0;columns<=this->columns;columns++)
		{
			outputMatrix.PushVals(rows,columns,-(this->matrix[rows][columns]));
		}
	}
//It must return a MATRIX that is negative of the calling matrix operator
return outputMatrix;
}

MATRIX MATRIX::Inverse(string execType) const
{	//
	MATRIX inverse("zero",this->Rows()+1,2*(this->Columns()+1));//Creates the Adjunct matrix keeping the rowConstant.
	assert((inverse.Rows()== this->Rows()) && (inverse.Columns()==(2*this->Columns())+1));
	
	// Augment the Identity Matrix.
	inverse.AssignSubMatrix(0,this->Rows(),0,this->Columns(),*this);
	inverse.AssignSubMatrix(0,inverse.Rows(),(inverse.Columns()+1)/2,inverse.Columns(),MATRIX("identity",this->Rows()+1,this->Columns()+1));
	
	//Check DiagElement ==0,
	unsigned long rank = this->Rows();
	for(unsigned long check=0;check<=rank;check++)
	{
		if(this->Element(check,check)==0) 
		{
			cout<<"No Inverse Exists"<<endl;
			exit(1);
		}
	}
	
	//########################## Actual Matrix Computation #####################################
	//Apply Gauss Jordan elimination.
	#pragma omp parallel if(execType == "parallel")
	#pragma omp for schedule(dynamic,rank/omp_get_num_threads()) 
	for(unsigned long n = 0; n<=rank;n++)//Splitting upper loops.
	{
		//First divide the Row[n] by A[n][n].(Rowise.)
		float norm = inverse.Element(n,n);
		for(unsigned long _cols = 0;_cols<=inverse.Columns();_cols++) 
		{
			inverse.PushVals(n,_cols,inverse.Element(n,_cols)/norm);
		}
		
		//Then make the column[n] == 0 except at when row == column for n.
		for(unsigned long _rows = 0;_rows<=this->Rows();_rows++)//j
		{	if(n!=_rows)
			{
				float multiplier = inverse.Element(_rows,n);
				for(unsigned long _k=0;_k<=inverse.Columns();_k++)
				{
					float temp = inverse.Element(_rows,_k)-(multiplier*(inverse.Element(n,_k)));
					inverse.PushVals(_rows,_k,temp);
				}
			}
			
		}
		//cout<<" After the rank = "<<n<<endl;
		//inverse.PrintMatrix();
		//Repeat the Procedure for all the Rows and Columns(Upto the Entire Value)
	}
	//cout<<endl<<endl;
	//cout<<"Step Wise Inverse Checks"<<endl;
	//inverse.PrintMatrix();

	return inverse.SubMatrix(0,inverse.Rows(),(inverse.Columns()+1)/2,inverse.Columns());
}

float MATRIX::Element(unsigned long rows,unsigned long columns) const
{	
	assert(rows<=this->rows);
	assert(columns<=this->columns);
	return matrix[rows][columns];
}

MATRIX::MATRIX(string matrixType="upperTriangular",unsigned long rows=0,unsigned long columns=0)
{	
	this->rows = rows-1;//Assigned rows and columns with starting at 0 index,great.
	this->columns = columns-1;
	srand(1);
    for(unsigned long rowIndex =0;rowIndex<=this->rows;rowIndex++)
	{	
		vector<float> columnValues;
		for(unsigned long columnIndex = 0;columnIndex<=this->columns;columnIndex++)
		{
			if(matrixType == "identity")
			{	//Kind of opposite it to parallelization, but useful for prototyping
				if(columnIndex==rowIndex) columnValues.push_back(1.0);//Ensures Creation of Identity Matrix
				else columnValues.push_back(0.0);
			}
			else if((columnIndex <rowIndex) || (matrixType=="zero")) columnValues.push_back(0.0);
			else if(matrixType=="upperTriangular")columnValues.push_back(-(rand()%10)+1);//Ensures all types of matrices are column type
			else {
				cout<<"Incorrect Matrix Type: exiting"<<endl;
				exit(1);
			}
		}
		matrix.push_back(columnValues);
	}

//Matrix Reconditioning For Easier Computation:
	if(matrixType=="upperTriangular")
	{
		//Abs Computation
		vector<vector<float>>::iterator row_iter = this->matrix.begin();
		while(row_iter != this->matrix.end())
		{
			vector<float>::iterator col_iter = row_iter->begin();
			while(col_iter!=row_iter->end())
			{
				*col_iter = abs(static_cast<float>(*col_iter));
				col_iter++;
			}
			row_iter++;
		}
		//cout<<"Finished Doing the absolute calculation"<<endl;
		//Taking the sum across each column
		vector<float> sum_diagonal;
		for(unsigned long _cols = 0;_cols<=this->columns;_cols++)
		{	float temp=0;
			for(unsigned long _rows =0;_rows<=this->rows;_rows++)
			{
				temp+=this->matrix[_rows][_cols];
			}
			sum_diagonal.push_back(temp);
		}
		
		//adding those sums to each of the diagonal matrix.
		for(unsigned long rank=0;rank<=this->rows;rank++)
		{
			this->matrix[rank][rank] = sum_diagonal[rank];
		}
		
	}


}

void MATRIX::PrintMatrix()
{
	assert(this->matrix.size()==((this->rows+1))); 
	assert(this->matrix[0].size()==((this->columns+1)));
	vector<vector<float>>::iterator rowIter = this->matrix.begin();
	
	while(rowIter!= this->matrix.end())
	{
		vector<float>::iterator columnIter = rowIter->begin();
		while(columnIter!=rowIter->end())
		{
			cout<<std::setw(4);
			cout<<(*columnIter)<<" ";
			columnIter++;
		}
		rowIter++;
		cout<<endl;
	}
}

void MATRIX::PushVals(unsigned long rowIndex,unsigned long columnIndex,float value)
{
	assert(rowIndex<=this->rows);
	assert(columnIndex<=this->columns);
	this->matrix[rowIndex][columnIndex]=value;
}

MATRIX Multiply(const MATRIX& mat1,const MATRIX& mat2)
{
	if(mat1.Columns() != mat2.Rows()) throw "Cannot be Multiplied, rows != Columns";
	unsigned long outputRows = mat1.Rows();
	unsigned long outputColumns =  mat2.Columns();
	unsigned long k = mat1.Columns();
	//cout<<"OutputColumns = "<<outputColumns<<endl;
	//cout<<"OutputRows = "<<outputRows<<endl;

	MATRIX outputMatrix("zero",outputRows+1,outputColumns+1);
   for(unsigned long rows =0;rows<=outputRows;rows++)
	{
		for(unsigned long columns = 0;columns<=outputColumns;columns++)
		{	float temp = 0;
			for(unsigned long commonIndex=0;commonIndex<=k;commonIndex++)
			{	
				temp += ((mat1.Element(rows,commonIndex))*(mat2.Element(commonIndex,columns)));
			}
			outputMatrix.PushVals(rows,columns,temp);
		}
	}

return outputMatrix;
}

MATRIX computeInverse(const MATRIX& input)
{	
	assert(input.Rows()==input.Columns());//Check if the input matrix is a square Matrix.
	
	unsigned long n = input.Rows();//Create a Matrix of Size n.
	MATRIX output("zero",n+1,n+1);//Create a Zero Matrix of size of input
	if(n<=512) //If size of Matrix n == 4 it creates an inverse matrix and returns it ,wrt to 0 index.
	{
		output = input.Inverse("parallel");//Finds the inverse of incoming matrix
	}
else
	{	//Eg: 6x6 matrix akak will indexing as 5x5 matrice (Split into 3x3) 
		//Generalize it, since by using unsigned long we can have, 5/2 = 2,go ahead
		#pragma omp parallel 
		{
		#pragma omp single
		{ 
	    unsigned long ni = input.Rows()/2;
		#pragma omp task shared(output)
		{
		output.AssignSubMatrix(0,ni,0,ni,computeInverse(input.SubMatrix(0,ni,0,ni)));//indp task
		}
		#pragma omp task shared(output)
		{
		output.AssignSubMatrix(ni+1,n,ni+1,n,computeInverse(input.SubMatrix(ni+1,n,ni+1,n)));//independent task
	    }
	    #pragma omp taskwait
	    {
		MATRIX R11_inv = -(output.SubMatrix(0,ni,0,ni));
		MATRIX R12 = input.SubMatrix(0,ni,ni+1,n);
		MATRIX R22_inv = output.SubMatrix(ni+1,n,ni+1,n);
		output.AssignSubMatrix(0,ni,ni+1,n,Multiply(Multiply(R11_inv,R12),R22_inv));
	    }
	    }
	}
		}
	return output;
}


MATRIX computeInverse_serial(const MATRIX& input)
{	
	assert(input.Rows()==input.Columns());//Check if the input matrix is a square Matrix.
	
	unsigned long n = input.Rows();//Create a Matrix of Size n.
	MATRIX output("zero",n+1,n+1);//Create a Zero Matrix of size of input
	if(n<=2) //If size of Matrix n == 4 it creates an inverse matrix and returns it ,wrt to 0 index.
	{
		output = input.Inverse("serial");//Finds the inverse of incoming matrix
	}
	else
	{	//Eg: 6x6 matrix akak will indexing as 5x5 matrice (Split into 3x3) 
		//Generalize it, since by using unsigned long we can have, 5/2 = 2,go ahead
		unsigned long ni = input.Rows()/2;
		output.AssignSubMatrix(0,ni,0,ni,computeInverse_serial(input.SubMatrix(0,ni,0,ni)));
		output.AssignSubMatrix(ni+1,n,ni+1,n,computeInverse_serial(input.SubMatrix(ni+1,n,ni+1,n)));
		MATRIX R11_inv = -(output.SubMatrix(0,ni,0,ni));
		MATRIX R12 = input.SubMatrix(0,ni,ni+1,n);
		MATRIX R22_inv = output.SubMatrix(ni+1,n,ni+1,n);
		output.AssignSubMatrix(0,ni,ni+1,n,Multiply(Multiply(R11_inv,R12),R22_inv));

	}
	return output;
}

using namespace::chrono;
int main(int argc ,char **argv)
{	
assert(argc ==3);
unsigned long rank = (unsigned long)atoi(argv[1]);
cout<<"Rank :"<<rank<<endl;
MATRIX A("upperTriangular",rank,rank);


auto start =high_resolution_clock::now();
MATRIX Inv_A =computeInverse(A);
auto stop =high_resolution_clock::now();
auto duration_parallel = duration_cast<milliseconds>(stop-start);
cout<<"Time Taken for Parallel Inverse Computation :"<<duration_parallel.count()<<" ms"<<endl;

auto start_serial =high_resolution_clock::now();
MATRIX Inv_A_serial=computeInverse_serial(A);
auto stop_serial =high_resolution_clock::now();
auto duration_serial = duration_cast<milliseconds>(stop_serial-start_serial);
cout<<"Time Taken for Serial Inverse Computation :"<<duration_serial.count()<<" ms"<<endl;


cout<<endl<<endl;
float speedup = (float)(duration_serial.count())/(duration_parallel.count());
cout<<"Speedup :"<<speedup<<endl;
cout<<"Efficiency :"<<speedup/(atof(argv[2]));
cout<<endl;
cout<<"Do you want to display output? yes =1 no=0"<<endl;
int inp;
cin>>inp;
if(inp = '1') Inv_A_serial.PrintMatrix();
//Multiply(A,Inv_A).PrintMatrix();
return EXIT_SUCCESS;
}
