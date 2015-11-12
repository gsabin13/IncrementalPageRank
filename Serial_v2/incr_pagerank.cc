#include <bits/stdc++.h>

#include "common.h"
#include "sparse_vector.h"
#include "dense_vector.h"
#include "sparse_matrix_csc.h"
#include "incr_pagerank.h"

#include "COO.h"
#include "tools/ntimer.h"

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
		//x_old = x
		x_old.VectorCopy(x);
		
		//x = alpha*x_const + (1-alpha)*(A*x)
		DenseVector A_dot_x_times_1_minus_alpha = (A_mat.multiplyWithDenseVector(x));
		A_dot_x_times_1_minus_alpha.InPlaceScalarMultiply(1-alpha);
		x.DenseVectorSum(alpha_times_xconst, A_dot_x_times_1_minus_alpha);
		
		
		//cout<<"PR_X_inter after iteration "<<i<<":\n"<<x<<endl;
		if(IsThresholdCrossed(x_old, x, threshold)){
			break;
		}
	}
	cout<<"\ntotal num iters "<<i<<endl; 

	return x;

}


DenseVector IncrementalPageRankSpMSpV(SparseMatrixCSC& A_plus_delta_A, SparseMatrixCSC& delta_A_mat, DenseVector& x_0, double alpha, double threshold, int num_iters){

	double copy_time =0.0, inplace_add_time=0.0, mul_time=0.0, inplace_scalar_mul_time=0.0, thresholdcheck_time=0.0;

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
			
        double pre = time_in_mill_now();
		x_sparse_old.SparseVectorCopyFrom(x_sparse);
		double post = time_in_mill_now();
		copy_time += (pre-post);
		//SpMSpV op
		//x_sparse = (1-alpha)*(A_plus_delta_A*(lambda_ + x_sparse));
		pre = time_in_mill_now();
		x_sparse.inPlaceAdd(lambda_);
		post = time_in_mill_now();
		inplace_add_time += (pre-post);
		
		pre = time_in_mill_now();
		A_plus_delta_A.inPlaceMultiplyWithSparseVector(x_sparse);
		post = time_in_mill_now();
		mul_time += (pre-post);

		pre = time_in_mill_now();
		x_sparse.inPlaceScalarMultiply(1-alpha);
		post = time_in_mill_now();
		inplace_scalar_mul_time += (pre-post);

#if ENABLE_DEBUG		
		cout<<"1-alpha)* A_plus_delta_A*(xsparse + lambda_) after iteration #"<<i<<": ";
		(x_sparse.getDenseVectorForm()).Print();
#endif
		pre = time_in_mill_now();
		if(IsThresholdCrossed(x_sparse_old, x_sparse, threshold))
			break;
		post = time_in_mill_now();
		thresholdcheck_time += (pre-post);
	}
	cout<<"\ntotal num iters "<<i<<endl; 
	
	cout<<" copy_time" << copy_time<<endl
		<<" inplace_add_time" << inplace_add_time<<endl
		<<" mul_time" << mul_time<<endl
		<<" inplace_scalar_mul_time" << inplace_scalar_mul_time<<endl
		<<" thresholdcheck_time" << thresholdcheck_time<<endl
		;

	//result = x_0 + lambda_ + x_sparse;
	result.inPlaceAdd(x_sparse);
	result.inPlaceAdd(lambda_);
	result.inPlaceAdd(x_0);

	return result.getDenseVectorForm();

}



int main(int argc, char *argv[])
{
	if(argc != 3){
		cout<<"Incorrect format. Format should be: ./a.out <filename> <output_directory>"<<endl;
		return 0;
	}
	
/********************************************************
Creation of files for multiple outputs that are printed during the program 
********************************************************/
	char* filepath = argv[1];
	char * output_directory_path = argv[2];
	string filepath_st(filepath);
	string file_name_st = string(output_directory_path) + "/" + filepath_st.substr(filepath_st.find_last_of("/")+1, filepath_st.length());
	
	char timing_file[100];
	strcpy(timing_file, file_name_st.c_str());
	strcat(timing_file, ".timing");


	string orig_file_name_st = file_name_st +".orig";
	char orig_intermediate_edgelist_file[100];
	strcpy(orig_intermediate_edgelist_file, orig_file_name_st.c_str());
	strcat(orig_intermediate_edgelist_file,".intermediate.edgelist");
	char orig_output_file[100];
	strcpy(orig_output_file, orig_file_name_st.c_str());
	strcat(orig_output_file,".out");

	string incr_file_name_st = file_name_st + ".incr";
	char incr_intermediate_edgelist_file[100];
	strcpy(incr_intermediate_edgelist_file, incr_file_name_st.c_str());
	strcat(incr_intermediate_edgelist_file,".intermediate.edgelist");
	char incr_output_file[100];
	strcpy(incr_output_file, incr_file_name_st.c_str());
	strcat(incr_output_file,".out");

	ofstream fout(timing_file);
	const int NUM_ITERS = 300;
	int NUMROWS, NUMCOLS;
	const double THRESHOLD = 0.0000000001;

/********************************************************
Reading from input matrix and creating an adjacency list 
********************************************************/
	COO cooAt;
	
  	cooAt.readSNAPFile(filepath, false);
	//printf("cooAt before removing Duplicates rows=%d cols=%d nnz=%d\n", cooAt.rows, cooAt.cols, cooAt.nnz);
	cooAt.orderedAndDuplicatesRemoving();
	//printf("cooAt after removing Duplicates rows=%d cols=%d nnz=%d\n", cooAt.rows, cooAt.cols, cooAt.nnz);

	vector<Node> adj_list = CreateAdjacencyListFromEdgeList(cooAt);

	PrintEdgeListFromAdjListToFile(adj_list,adj_list.size(), cooAt.nnz, orig_intermediate_edgelist_file);

/********************************************************
Creation of SparseMatrix and initial vector objects 
********************************************************/
	SparseMatrixCSC A_mat = CreateMatrixFromAdjList(NUMROWS, NUMCOLS, adj_list);
	vector<valtype> x_init_vec(NUMROWS,1.0/double(NUMROWS));
	DenseVector x_init(x_init_vec);
	DenseVector x_0(x_init_vec);

#if ENABLE_PRINT
	cout<<"Original matrix transpose: \n";
	A_mat.PrintMatrix();

	std::cout << "Initial PageRank vector: ";
	x_init.Print();
#endif
	

/********************************************************
PageRank() call - does PageRank on the original graph 
********************************************************/

    double orig_pre = time_in_mill_now();
	DenseVector result = PageRank(A_mat, x_init, x_0, ALPHA, THRESHOLD, NUM_ITERS);
	double orig_post = time_in_mill_now();
    fout << file_name_st<<", "
              << orig_post-orig_pre
              << "msec, ";

	PrintOutputToFile(orig_output_file, result);

/********************************************************
Modification of graph and creation of new SparseMatrix objects
********************************************************/

	//Function that modifies graph so that we can try out incremental PR on it
	int num_edges_added = ModifyGraph(adj_list);
	PrintEdgeListFromAdjListToFile(adj_list, adj_list.size(), cooAt.nnz + num_edges_added, incr_intermediate_edgelist_file);

	SparseMatrixCSC A_plus_delta_A = CreateMatrixFromAdjList(NUMROWS, NUMCOLS, adj_list);
	SparseMatrixCSC deltaA = A_plus_delta_A.subtract(A_mat);


#if ENABLE_PRINT
	cout<<"New matrix transpose: \n";
	A_plus_delta_A.PrintMatrix();
	
	cout<<"Delta matrix transpose: \n";
	deltaA.PrintMatrix();
#endif
	
/********************************************************
IncrementalPageRankSpMSpV() call - does SpmSpV based PageRank on the modified graph
********************************************************/
	double incr_pre = time_in_mill_now();
	DenseVector result_incremental = IncrementalPageRankSpMSpV(A_plus_delta_A, deltaA, result, ALPHA, THRESHOLD, NUM_ITERS);
	double incr_post = time_in_mill_now();
	fout << incr_post-incr_pre
              << "msec, ";

	PrintOutputToFile(incr_output_file, result_incremental);

/********************************************************
PageRank() call - does PageRank on the modified graph from scratch
********************************************************/
	double incr_scratch_pre = time_in_mill_now();
	DenseVector result_scratch = PageRank(A_plus_delta_A, x_init, result, ALPHA, THRESHOLD, NUM_ITERS);
	double incr_scratch_post = time_in_mill_now();
	fout << incr_scratch_post-incr_scratch_pre
              << "msec\n";

	strcat(incr_output_file,".scratch");
	PrintOutputToFile(incr_output_file, result_scratch);	

}





/*
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

	DenseVector result2 = IncrementalPageRank(mat2, deltaA, result, ALPHA, threshold, num_iters);

	result2.Print();

	DenseVector result3 = PageRank(mat2, x_const, result, ALPHA, threshold, num_iters);
	result3.Print();
	cout<<"hello world \n";
}
*/