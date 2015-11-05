#include <iostream>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <vector>
#include <bits/stdc++.h>

#include <fstream>

#define ENABLE_PRINT 0

#define ALPHA 0.15

using namespace Eigen;
using namespace std;

typedef double matrixValType;

typedef SparseMatrix<matrixValType> MyMatrix;

void PrintMatrix(MyMatrix& mat){
	std::cout<<"\n Number of Rows = "<<mat.outerSize()<<"\n";
	for (int k=0; k<mat.outerSize(); ++k){
		std::cout<<"Column "<<k<<" : ";
		for (MyMatrix::InnerIterator it(mat,k); it; ++it)
		{
			std::cout<<it.value()<<" row#: "<<it.row()<< " col#: "<<it.col()<<" index#: "<<it.index()<<"\n";
		}
		std::cout<<"\n";
	}
}

matrixValType getSumOfVectorElems(VectorXd& x){
	matrixValType sum = 0.0;
	for(int i=0;i<x.rows();i++){
		sum+=x(i);
	}
	return sum;
}

bool IsThresholdCrossed(VectorXd& x_old, VectorXd& x_new, double threshold){
	for(int i=0;i<x_old.rows();i++){
		//cout<<"abs " << abs(x_old(i)-x_new(i))<<" ";
		if(abs(x_old(i)-x_new(i)) > threshold){

			return false;
		}
	}
	//cout<<"\nTHRESHOLD CROSSED! "<<threshold<<"\n";
	return true;

}


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

struct Node{
	int id;
	//boolean to indicate whether it is a dangling node
	bool isDangling;
	vector<int> neighbours;

	Node(): isDangling(false){}
};

vector<Node> CreateAdjacencyListFromEdgeList(string filename){
	ifstream fin(filename);

	int num_nodes, num_edges;
	fin>>num_nodes>>num_edges;

	vector<Node> adj_list(num_nodes);

	for(int i =0;i<num_edges;i++){
		int a,b;
		fin>>a>>b;
		adj_list[a].neighbours.push_back(b);
	}

	return adj_list;
}


MyMatrix CreateMatrixFromFile(int& NUMROWS, int& NUMCOLS, string filename){

	vector<Node> adj_list = CreateAdjacencyListFromEdgeList(filename);
	
	NUMROWS = adj_list.size();
	NUMCOLS = NUMROWS;
	
	std::vector<int> rows;
	std::vector<int> cols;

	for(int i=0;i<NUMROWS;i++){
		//cout<<i<<": ";
		//if no outer link, assume it links to all - this is done in original PR as well!
		if(adj_list[i].neighbours.size() == 0){
			adj_list[i].isDangling = true;
			for(int j=0;j<NUMROWS;j++){
				adj_list[i].neighbours.push_back(j);
			}
		}
		for(int j=0;j<adj_list[i].neighbours.size();j++){
			rows.push_back(i);cols.push_back(adj_list[i].neighbours[j]);
			//cout<<adj_list[i].neighbours[j]<<" ";
		}

		//cout<<endl;
	}

	int nnz = rows.size();

	int estimation_of_entries = nnz;
	typedef Eigen::Triplet<matrixValType> T;
	std::vector<T> tripletList;
	tripletList.reserve(estimation_of_entries);
	for(int i=0;i<nnz;i++)
	{
		tripletList.push_back(T(rows[i],cols[i],1.0/double(adj_list[rows[i]].neighbours.size())));
	}

	MyMatrix mat(NUMROWS, NUMCOLS);
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
	
	return mat;
}



int main()
{

	int NUM_ITERS = 300;
	int NUMROWS, NUMCOLS;
	const double THRESHOLD =0.0000000001;

	string filename = "edgelist.in";
	MyMatrix mat = CreateMatrixFromFile(NUMROWS, NUMCOLS, filename);
	
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

	x = x_init;

	std::cout << "Initial PageRank vector:\n" << x_init << std::endl;
	
	VectorXd result = PageRank(mat_transpose, x_init,x, ALPHA, THRESHOLD, NUM_ITERS);

	std::cout << "Here is the PageRank vector after initial PageRank:\n" << result << std::endl;

	string filename_delta = "delta_edgelist.in";
	MyMatrix mat_new = CreateMatrixFromFile(NUMROWS, NUMCOLS, filename_delta);
	
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

	VectorXd result_incremental = IncrementalPageRank(mat_new_transpose, delta_mat, result, ALPHA, THRESHOLD, NUM_ITERS);
	std::cout << "PageRank vector after incremental PR:\n" << result_incremental << std::endl;

	VectorXd result_scratch = PageRank(mat_new_transpose, x_init, result, ALPHA, THRESHOLD, NUM_ITERS);
	std::cout << "PageRank vector after PR from scratch:\n" << result_scratch << std::endl;

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