#ifndef INCR_PAGERANK_H
#define INCR_PAGERANK_H

#include "common.h"
#include "sparse_vector.h"
#include "dense_vector.h"
#include "sparse_matrix_csc.h"
#include "incr_pagerank.h"

#include "COO.h"
#include "tools/ntimer.h"

struct Node{
	int id;
	//boolean to indicate whether it is a dangling node
	bool isDangling;
	vector<int> neighbours;

	Node(): isDangling(false){}
};


bool IsThresholdCrossed(DenseVector& x_old, DenseVector& x_new, double threshold){
	for(int i=0;i<x_old.size();i++){
		//cout<<"abs " << abs(x_old(i)-x_new(i))<<" ";
		if(abs(x_old[i]-x_new[i]) > threshold){

			return false;
		}
	}
	cout<<"\nTHRESHOLD CROSSED! "<<threshold<<"\n";
	return true;

}

bool IsThresholdCrossed(SparseVector& x_sparse_old, SparseVector& x_sparse_new, double threshold){
	//case when x_sparse starts off as 0 vector
	if(x_sparse_old.indexes_sparse.size() == 0){
		return false;
	}
	for(int i=0;i<x_sparse_old.indexes_sparse.size(); i++){
		int index_being_compared = x_sparse_old.indexes_sparse[i];
		if(abs(x_sparse_old.vals_dense[index_being_compared]-x_sparse_new.vals_dense[index_being_compared]) > threshold){

			return false;
		}
	}
	cout<<"\nTHRESHOLD CROSSED! "<<threshold<<"\n";
	return true;
}


vector<Node> CreateAdjacencyListFromEdgeList(COO& cooAt){
	
	//ifstream fin(filename);

	int num_nodes, num_edges;
	//fin>>num_nodes>>num_edges;
	num_nodes = cooAt.rows;
	num_edges = cooAt.nnz;

	vector<Node> adj_list(num_nodes);


	for(int i =0;i<num_edges;i++){
		int a,b;
		a = cooAt.cooRowIndex[i]; b = cooAt.cooColIndex[i];
		adj_list[a].neighbours.push_back(b);
	}

	return adj_list;
}

//Note that the matrix that is being returned is the transpose of the (PageRank) adjacency matrix of the graph. It is in CSC format
SparseMatrixCSC CreateMatrixFromAdjList(int& NUMROWS, int& NUMCOLS, vector<Node>& adj_list){
	
	NUMROWS = adj_list.size();
	NUMCOLS = NUMROWS;
	
	vector<int> cols(NUMROWS+1);
	vector<int> rows;
	vector<valtype> matrix_vals;

	cols[0] = 0;
	int current_cols_idx = 1;

	for(int i=0;i<NUMROWS;i++){
		//if no outer link, assume it links to all - this is done in original PR as well!
		if(adj_list[i].neighbours.size() == 0){
			adj_list[i].isDangling = true;
			for(int j=0;j<NUMROWS;j++){
				adj_list[i].neighbours.push_back(j);
			}
		}
		//sorting so that insertion into cols, rows, vals vectors of SparseMatrixCSC is more convenient
		sort(adj_list[i].neighbours.begin(), adj_list[i].neighbours.end());
		
		int num_neighbours = adj_list[i].neighbours.size();
		for(int j=0;j<num_neighbours;j++){
			rows.push_back(adj_list[i].neighbours[j]);
			matrix_vals.push_back(1.0/valtype(num_neighbours));

		}
		cols[current_cols_idx] = cols[current_cols_idx-1] + num_neighbours;
		current_cols_idx++;

	}

	SparseMatrixCSC mat(cols, rows, matrix_vals);
	
	return mat;
}

bool RandRowColExistsInGraph(int rand_row, int rand_col, vector<Node>& adj_list){

	return adj_list[rand_row].neighbours.size() != 0 && 
			(std::find(adj_list[rand_row].neighbours.begin(), adj_list[rand_row].neighbours.end(), rand_col) 
				!= adj_list[rand_row].neighbours.end());

}

//TODO: CHANGE THIS TO ACCEPT MULTIPLE TYPES OF MODIFICATION (like the percentage of rows that are being changed)
int ModifyGraph(vector<Node>& adj_list){
	int num_edges_added = 0;
	//right now, adding an edge to a random 10% of nodes
	int num_rows = adj_list.size();
	for(int i=0;i<num_rows/10;i++){
		int rand_row = rand()%num_rows;
		int rand_col = rand()%num_rows;
		while(RandRowColExistsInGraph(rand_row, rand_col, adj_list)){
			rand_row = rand()%num_rows;
			rand_col = rand()%num_rows;
		}
		adj_list[rand_row].neighbours.push_back(rand_col);
		num_edges_added++;
	}
	return num_edges_added;
}


void PrintEdgeListFromAdjListToFile(vector<Node>& adj_list, int num_nodes, int num_edges, char filename[]){

	ofstream fout(filename);

	fout<<num_nodes<<" "<<num_edges<<endl;
	for(int i=0;i<num_nodes;i++){
		for(int j=0;j<adj_list[i].neighbours.size();j++){
			fout<<i<<" "<<adj_list[i].neighbours[j]<<endl;
		}
	}
}

void PrintOutputToFile(char filename[], DenseVector& result){
	ofstream fout(filename);
	for(int i=0;i<result.size();i++){
			fout<<result[i]<<"\n";
	}
	fout<<endl;
}



#endif //INCR_PAGERANK_H