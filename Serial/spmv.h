#ifndef SPMV_H
#define SPMV_H

#include <iostream>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <vector>
#include <bits/stdc++.h>

#include <fstream>

#define ENABLE_PRINT 0

using namespace Eigen;
using namespace std;


typedef double matrixValType;

typedef SparseMatrix<matrixValType> MyMatrix;

struct Node{
	int id;
	//boolean to indicate whether it is a dangling node
	bool isDangling;
	vector<int> neighbours;

	Node(): isDangling(false){}
};

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
	
	typedef Eigen::Triplet<matrixValType> T;
	std::vector<T> tripletList;

	for(int i=0;i<NUMROWS;i++){
		//if no outer link, assume it links to all - this is done in original PR as well!
		if(adj_list[i].neighbours.size() == 0){
			adj_list[i].isDangling = true;
			for(int j=0;j<NUMROWS;j++){
				adj_list[i].neighbours.push_back(j);
			}
		}
		for(int j=0;j<adj_list[i].neighbours.size();j++){
			tripletList.push_back(T(i,adj_list[i].neighbours[j],1.0/double(adj_list[i].neighbours.size())));
		}

	}

	MyMatrix mat(NUMROWS, NUMCOLS);
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
	
	return mat;
}

#endif