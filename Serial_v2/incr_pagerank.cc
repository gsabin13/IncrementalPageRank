#include <bits/stdc++.h>

#include "common.h"
#include "sparse_vector.h"
#include "dense_vector.h"
#include "sparse_matrix_csc.h"
#include "incr_pagerank.h"

#define ALPHA 0.15

#define ENABLE_DEBUG 0

using namespace std;


DenseVector PageRank(SparseMatrixCSC& A_mat, DenseVector& x_const, DenseVector& x_0, double alpha, double threshold, int num_iters){

	DenseVector x(x_0);
	//x = x_0;
	DenseVector x_old(x_0.size());

	int i;

	DenseVector alpha_times_xconst(x_const);
	alpha_times_xconst.InPlaceScalarMultiply(alpha);	
	//x.Print();
	for( i=0;i<num_iters;i++ ){
		//cout<<"\niteration "<<i<<" : ";
		//x_old = x
		x_old.VectorCopy(x);
		//cout<<"x: ";
		//x.Print();
		//SpMV op
		DenseVector A_dot_x_times_1_minus_alpha = (A_mat.multiplyWithDenseVector(x));
		//cout<<"A_dot_x: ";
		//A_dot_x_times_1_minus_alpha.Print();
		
		A_dot_x_times_1_minus_alpha.InPlaceScalarMultiply(1-alpha);
		//cout<<"A_dot_x_times_1_minus_alpha: ";	
		//A_dot_x_times_1_minus_alpha.Print();
		//x = alpha*x_const + (1-alpha)*(A*x)
		x.DenseVectorSum(alpha_times_xconst, A_dot_x_times_1_minus_alpha);
		
		
		//cout<<"PR_X_inter after iteration "<<i<<":\n"<<x<<endl;
		if(IsThresholdCrossed(x_old, x, threshold)){
			break;
		}
	}
	cout<<"\ntotal num iters "<<i<<endl; 

	return x;

}


DenseVector IncrementalPageRank(SparseMatrixCSC& A_plus_delta_A, SparseMatrixCSC& delta_A_mat, DenseVector& x_0, double alpha, double threshold, int num_iters){

	SparseVector result(x_0.size());

	//SpMSpV op
	//lambda_ = (1-alpha)*(delta_A_mat*x_0);
	SparseVector lambda_ = x_0;
	delta_A_mat.inPlaceMultiplyWithSparseVector(lambda_);
	lambda_.inPlaceScalarMultiply(1-alpha);
	
	//initial x_sparse contains all 0s
	SparseVector x_sparse(x_0.size());

	SparseVector x_sparse_old(x_0.size());
	int i;

	for(i=0;i<num_iters;i++){
		//x_sparse_old = x_sparse
		x_sparse_old.SparseVectorCopyFrom(x_sparse);
		//SpMSpV op
		//x_sparse = (1-alpha)*(A_plus_delta_A*(lambda_ + x_sparse));
		x_sparse.inPlaceAdd(lambda_);
		A_plus_delta_A.inPlaceMultiplyWithSparseVector(x_sparse);
		x_sparse.inPlaceScalarMultiply(1-alpha);

#if ENABLE_DEBUG		
		cout<<"1-alpha)* A_plus_delta_A*(xsparse + lambda_) after iteration #"<<i<<": ";
		(x_sparse.getDenseVectorForm()).Print();
#endif
		if(IsThresholdCrossed(x_sparse_old, x_sparse, threshold))
			break;

	}
	cout<<"\ntotal num iters "<<i<<endl; 
	
	//result = x_0 + lambda_ + x_sparse;
	result.inPlaceAdd(x_sparse);
	result.inPlaceAdd(lambda_);
	result.inPlaceAdd(x_0);

	return result.getDenseVectorForm();

}

int main(){

	int c1[] = {0, 1,3,6,7,8};
	int r1[] = {1,0, 4,0,1,3,4,2};
	valtype v1[] = {1, 0.5, 0.5, 0.333333333, 0.333333333, 0.333333333, 1,  1};

	vector<int> cols1(c1, c1+6);
	vector<int> rows1(r1, r1+8);
	vector<valtype> vals1(v1, v1+8);

	int c2[] = {0, 1,3,6,8,9};
	int r2[] = {1,0, 4,0,1,3,1,4,2};
	valtype v2[] = {1, 0.5, 0.5, 0.333333333, 0.333333333, 0.333333333, 0.5, 0.5,  1};


	vector<int> cols2 (c2, c2+6);
	vector<int> rows2(r2, r2+9);
	vector<valtype> vals2(v2, v2+9);

	SparseMatrixCSC mat1(cols1, rows1, vals1);
	SparseMatrixCSC mat2(cols2, rows2, vals2);
	
	mat1.PrintMatrix();

	vector<valtype> x_c(5,0.2);
	DenseVector x_const(x_c);
	DenseVector x_0(x_c);
	
	const double threshold = 1e-9;
	const int num_iters = 200;

	DenseVector result = PageRank(mat1, x_const, x_0, ALPHA, threshold, num_iters);

	SparseMatrixCSC deltaA = mat2.subtract(mat1);
	mat2.PrintMatrix();
	deltaA.PrintMatrix();
	result.Print();

	result = IncrementalPageRank(mat2, deltaA, result, ALPHA, threshold, num_iters);

	result.Print();

	DenseVector result2 = PageRank(mat2, x_const, x_0, ALPHA, threshold, num_iters);
	result2.Print();
	cout<<"hello world \n";
}