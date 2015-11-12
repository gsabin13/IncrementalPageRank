#include "sparse_matrix_csc.h"

SparseMatrixCSC SparseMatrixCSC::subtract(SparseMatrixCSC& B){
	if(this->getNumRows() != B.getNumRows() || this->getNumCols() != B.getNumCols())
		cout<<"ERROR IN MATRIX SUBTRACTION - rows and cols don't match!!!\n";

	
	vector<int> C_cols;
	vector<int> C_rows;
	vector<valtype> C_matrix_vals;

	C_cols.push_back(0);
	int non_zero_count_C = 0;
	for(int col_num = 0; col_num < this->cols.size()-1;col_num++){
		
		int curr_rowidx_A = this->cols[col_num];
		int curr_rowidx_B = B.cols[col_num];

		while(curr_rowidx_A < this->cols[col_num+1] && curr_rowidx_B < B.cols[col_num+1]){
			int curr_row_A = this->rows[curr_rowidx_A];
			int curr_row_B = B.rows[curr_rowidx_B];

			valtype curr_val_A = this->matrix_vals[curr_rowidx_A];
			valtype curr_val_B = B.matrix_vals[curr_rowidx_B];
			//cout<<"curr_rowidx_A "<< curr_rowidx_A<< " "<<"curr_rowidx_B "<< curr_rowidx_B<<endl;
			//cout<<"curr_row_A "<< curr_row_A<< " "<<"curr_row_B "<< curr_row_B<< " "<<"curr_val_A "<< curr_val_A<< " "<<"curr_val_B "<< curr_val_B<< endl;
			if(curr_row_A == curr_row_B){
				//TODO: don't push back in this case! Pushing back isn't wrong but will not result in an ultra sparse matrix when doing a+delta_a - a. wont make a huge diff performance wise as delta_a is used only once
				//if doing so, remove the below 2 lines and move the non_zero_count_c++ to the else if and else conditions below
				C_rows.push_back(curr_row_A);
				C_matrix_vals.push_back(curr_val_A - curr_val_B);
				curr_rowidx_A++;
				curr_rowidx_B++;
			} else if(curr_row_A < curr_row_B){
				C_rows.push_back(curr_row_A);
				C_matrix_vals.push_back(curr_val_A);
				curr_rowidx_A++;
			} else {
				C_rows.push_back(curr_row_B);
				C_matrix_vals.push_back(-1.0 * curr_val_B);
				curr_rowidx_B++;
			}

			non_zero_count_C++;
		}

		while(curr_rowidx_A < this->cols[col_num+1]){
			int curr_row_A = rows[curr_rowidx_A];
			valtype curr_val_A = matrix_vals[curr_rowidx_A];
			C_rows.push_back(curr_row_A);
			C_matrix_vals.push_back(curr_val_A);
			curr_rowidx_A++;
			non_zero_count_C++;
		}
		while(curr_rowidx_B < B.cols[col_num+1]){
			int curr_row_B = rows[curr_rowidx_B];
			valtype curr_val_B = matrix_vals[curr_rowidx_B];
			C_rows.push_back(curr_row_B);
			C_matrix_vals.push_back(-1.0 * curr_val_B);
			curr_rowidx_B++;
			non_zero_count_C++;
		}

		//push nnz uptill this col
		C_cols.push_back(non_zero_count_C);

	}

	SparseMatrixCSC C(C_cols, C_rows, C_matrix_vals);
	return C;
}

void SparseMatrixCSC::inPlaceMultiplyWithSparseVector(SparseVector& V){
	V.clearDenseVals();
	for(uint i=0; i<V.indexes_sparse.size();i++){
		int v_index = V.indexes_sparse[i];
		valtype v_val = V.vals_sparse[i];

		for(int curr_rowidx_matrix = this->cols[v_index]; curr_rowidx_matrix < this->cols[v_index+1];curr_rowidx_matrix++){
			int curr_row_matrix = this->rows[curr_rowidx_matrix];
			V.vals_dense[curr_row_matrix] += this->matrix_vals[curr_rowidx_matrix]*v_val;
			V.insertIndexCorrespondingToNonZeroVal(curr_row_matrix);
		}
	}
	V.updateValsSparseBasedOnValsDenseAndIndexesSparse();
}

DenseVector SparseMatrixCSC::multiplyWithDenseVector(DenseVector& V){
	vector<valtype> result(V.size(),0);

	for(uint i=0; i<V.size();i++){
		
		valtype v_val = V[i];

		for(int curr_rowidx_matrix = this->cols[i]; curr_rowidx_matrix < this->cols[i+1];curr_rowidx_matrix++){
			int curr_row_matrix = this->rows[curr_rowidx_matrix];
			result[curr_row_matrix] += this->matrix_vals[curr_rowidx_matrix]*v_val;
		}
	}
	DenseVector result_dense(result);
	return result_dense;
}