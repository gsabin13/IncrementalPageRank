#ifndef SPARSE_MATRIX_CSC_H
#define SPARSE_MATRIX_CSC_H

#include "common.h"
#include "sparse_vector.h"
#include "dense_vector.h"

class SparseMatrixCSC{
	vector<int> cols;
	vector<int> rows;
	vector<valtype> matrix_vals;

public:
	SparseMatrixCSC(vector<int> cols1, vector<int> rows1, vector<valtype> matrix_vals1){
		cols = cols1;
		rows = rows1;
		matrix_vals = matrix_vals1;
	}

	size_t getNNZ(){
		return matrix_vals.size();
	}
	
	size_t getNumRows(){
		return cols.size() - 1;
	}

	size_t getNumCols(){
		return this->getNumRows();
	}

	//returns A-B;
	SparseMatrixCSC subtract(SparseMatrixCSC& B);

	//multiplies this matrix with a vector and stores the result in that same vector
	void inPlaceMultiplyWithSparseVector(SparseVector& V);

	DenseVector multiplyWithDenseVector(DenseVector& V);

	void PrintMatrix(){
		
		cout<<"Cols: ";
		for(int i=0;i<cols.size();i++){
			cout<<cols[i]<<" ";
		}
		cout<<endl;

		cout<<"Rows: ";
		for(int i=0;i<rows.size();i++){
			cout<<rows[i]<<" ";
		}
		cout<<endl;
		
		cout<<"vals: ";
		for(int i=0;i<matrix_vals.size();i++){
			cout<<matrix_vals[i]<<" ";
		}
		cout<<endl;		
	}

	


};



#endif