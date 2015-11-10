#include <iostream>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <vector>
#include <bits/stdc++.h>

#include "spmv.h"
#include "COO.h"
#include "tools/ntimer.h"

#define ENABLE_PRINT 0

#define ALPHA 0.15

using namespace Eigen;
using namespace std;


VectorXd PageRank(MyMatrix& A_mat, VectorXd& x_const, VectorXd& x_0, double alpha, double threshold, int num_iters){

	VectorXd x(x_0.rows());
	x = x_0;
	VectorXd x_old(x_0.rows());

	int i;
	for( i=0;i<num_iters && !IsThresholdCrossed(x_old,x, threshold);i++ ){
		x_old = x;
		//SpMV op
		x = (alpha*x_const) + (1-alpha)*(A_mat*x);
		//cout<<"PR_X_inter after iteration "<<i<<":\n"<<x<<endl;
	}
	cout<<"\ntotal num iters "<<i<<endl; 

	return x;

}

VectorXd IncrementalPageRank(MyMatrix& A_plus_delta_A, MyMatrix& delta_A_mat, VectorXd& x_0, double alpha, double threshold, int num_iters){

	VectorXd result(x_0.rows());
	VectorXd lambda_(x_0.rows());

	//SpMSpV op
	lambda_ = (1-alpha)*(delta_A_mat*x_0);

	VectorXd x_sparse(x_0.rows());
	for(int i=0;i<x_sparse.rows();i++){
		x_sparse(i) = 0;
	}

	VectorXd x_sparse_old(x_0.rows());
	int i;

	for(i=0;i<num_iters;i++){
		//SpMSpV op
		x_sparse = (1-alpha)*(A_plus_delta_A*(lambda_ + x_sparse));		
		if(IsThresholdCrossed(x_sparse_old, x_sparse, threshold))
			break;
		x_sparse_old = x_sparse;

	}
	cout<<"\ntotal num iters "<<i<<endl; 
	
	result = x_0 + lambda_ + x_sparse;
	return result;

}
	
void PrintEdgeListFromAdjList(vector<Node>& adj_list, int num_nodes, int num_edges, char filename[]){

	ofstream fout(filename);

	fout<<num_nodes<<" "<<num_edges<<endl;
	for(int i=0;i<num_nodes;i++){
		for(int j=0;j<adj_list[i].neighbours.size();j++){
			fout<<i<<" "<<adj_list[i].neighbours[j]<<endl;
		}
	}
}

void PrintOutput(char filename[], VectorXd& result){
	ofstream fout(filename);
	fout << result << endl;
}

int main(int argc, char *argv[])
{
	if(argc != 3){
		cout<<"Incorrect format. Format should be: ./a.out <filename> <output_directory>"<<endl;
		return 0;
	}
	
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

	//cout<<orig_intermediate_edgelist_file<<endl;

	COO cooAt;
	
  	cooAt.readSNAPFile(filepath, false);
	//printf("cooAt before removing Duplicates rows=%d cols=%d nnz=%d\n", cooAt.rows, cooAt.cols, cooAt.nnz);
	cooAt.orderedAndDuplicatesRemoving();
	//printf("cooAt after removing Duplicates rows=%d cols=%d nnz=%d\n", cooAt.rows, cooAt.cols, cooAt.nnz);

	vector<Node> adj_list = CreateAdjacencyListFromEdgeList(cooAt);


	PrintEdgeListFromAdjList(adj_list,adj_list.size(), cooAt.nnz, orig_intermediate_edgelist_file);

	int NUM_ITERS = 300;
	int NUMROWS, NUMCOLS;
	const double THRESHOLD =0.0000000001;
	
	MyMatrix mat = CreateMatrixFromAdjList(NUMROWS, NUMCOLS, adj_list);
	
	MyMatrix mat_transpose = MyMatrix(mat.transpose());

#if ENABLE_PRINT
	cout<<"Original matrix transpose: \n";
	PrintMatrix(mat_transpose);
#endif

	//d because double - need to change this if changing type
	VectorXd x(NUMROWS), x_init(NUMROWS);

	//initializing pagerank vector
	for(int i=0;i<NUMROWS;i++){
		x_init(i) = 1.0/double(NUMROWS);
	}

#if ENABLE_PRINT
	std::cout << "Initial PageRank vector:\n" << x_init << std::endl;
#endif
	x = x_init;	

	ofstream fout(timing_file);

    double orig_pre = time_in_mill_now();
	VectorXd result = PageRank(mat_transpose, x_init,x, ALPHA, THRESHOLD, NUM_ITERS);
	double orig_post = time_in_mill_now();
    fout << file_name_st<<", "
              << orig_post-orig_pre
              << "msec, ";

	PrintOutput(orig_output_file, result);


	//Function that modifies graph so that we can try out incremental PR on it
	int num_edges_added = ModifyGraph(adj_list);
	PrintEdgeListFromAdjList(adj_list, adj_list.size(), cooAt.nnz + num_edges_added, incr_intermediate_edgelist_file);

//	char* filename_delta = "delta_edgelist.in";
	MyMatrix mat_new = CreateMatrixFromAdjList(NUMROWS, NUMCOLS, adj_list);
	MyMatrix mat_new_transpose = MyMatrix(mat_new.transpose());

#if ENABLE_PRINT
	cout<<"New matrix transpose: \n";
	PrintMatrix(mat_new_transpose);
#endif

	MyMatrix delta_mat = mat_new_transpose - mat_transpose;

#if ENABLE_PRINT	
	cout<<"Delta matrix transpose: \n";
	PrintMatrix(delta_mat);
#endif
	
	double incr_pre = time_in_mill_now();
	VectorXd result_incremental = IncrementalPageRank(mat_new_transpose, delta_mat, result, ALPHA, THRESHOLD, NUM_ITERS);
	double incr_post = time_in_mill_now();
	fout << incr_post-incr_pre
              << "msec, ";

	PrintOutput(incr_output_file, result_incremental);

	double incr_scratch_pre = time_in_mill_now();
	VectorXd result_scratch = PageRank(mat_new_transpose, x_init, result, ALPHA, THRESHOLD, NUM_ITERS);
	double incr_scratch_post = time_in_mill_now();
	fout << incr_scratch_post-incr_scratch_pre
              << "msec\n";

	strcat(incr_output_file,".scratch");
	PrintOutput(incr_output_file, result_scratch);	


	//TODO: write an output validator between scratch and incremental


/*
	//In case the delta is provided as just modified edges and not the entire new graph

	vector<Node> delta_adj_list = CreateAdjacencyListFromEdgeList(filename_delta);

	std::vector<T> tripletList_delta_A;
	//tripletList.reserve();
	for(int i=0;i<nnz;i++)
	{
		tripletList_delta_A.push_back(T(rows[i],cols[i],1.0/double(adj_list[rows[i]].neighbours.size())));
	}

	for(int i=0;i<NUMROWS;i++){
		cout<<i<<": ";
		//only make changes to the matrix if there is an outer link in delta
		if(delta_adj_list[i].neighbours.size() != 0){
			//if dangling (a node that has no outlinks), then the adjacency list needs to be modified. The adjacency list for a dangling node is the list of all nodes. So we need to erase that and replace it with just the new nodes in the delta_adj_list
			if(adj_list[i].isDangling){
				adj_list[i].isDangling = false;
				adj_list[i].neighbours.erase(adj_list[i].neighbours.begin(), adj_list[i].neighbours.end());
				for(int j=0;j<delta_adj_list[i].neighbours.size();j++){
					adj_list[i].neighbours.push_back(j);
					//rows.push_back(i);cols.push_back(adj_list[i].neighbours[j]);
					cout<<adj_list[i].neighbours[j]<<" ";
				}

			}
			for(int j=0;j<adj_list[i].neighbours.size();j++){
				rows.push_back(i);cols.push_back(adj_list[i].neighbours[j]);
				cout<<adj_list[i].neighbours[j]<<" ";
			}
		} 

		cout<<endl;
	}
*/
	
}