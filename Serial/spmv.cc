#include <iostream>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <vector>
#include <bits/stdc++.h>

#include <fstream>

#define ALPHA 0.15

using namespace Eigen;
using namespace std;

typedef double matrixValType;

typedef SparseMatrix<matrixValType> MyMatrix;

void PrintMatrix(MyMatrix& mat){
	//SparseMatrix<double> mat(rows,cols);
	std::cout<<"\nmat.outersize() = "<<mat.outerSize()<<"\n";
	for (int k=0; k<mat.outerSize(); ++k){
		std::cout<<"row "<<k<<" : ";
		for (MyMatrix::InnerIterator it(mat,k); it; ++it)
		{
			std::cout<<it.value()<<" "<<it.row()<< " "<<it.col()<<" "<<it.index()<<"\n";
		    //it.row();   // row index
		   	//it.col();   // col index (here it is equal to k)
		    //it.index(); // inner index, here it is equal to it.row()
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
	//cout<<"x_old.rows()"<<x_old.rows()<<endl;
	for(int i=0;i<x_old.rows();i++){
		//cout<<"abs " << abs(x_old(i)-x_new(i))<<" ";
		if(abs(x_old(i)-x_new(i)) > threshold){

			return false;
		}
	}
	cout<<"\nTHRESHOLD CROSSED! "<<threshold<<"\n";
	return true;

}


VectorXd PageRank(MyMatrix& A_mat, VectorXd& x_const, VectorXd& x_0, double alpha, double threshold, int num_iters){

	VectorXd x(x_0.rows());
	x = x_0;
	VectorXd x_old(x_0.rows());
	//x_old = x_0;
	//double diff = threshold + 1;

	//matrixValType sum_old = getSumOfVectorElems(x_0);
	int i;
	for( i=0;i<num_iters && !IsThresholdCrossed(x_old,x, threshold);i++ ){
		//do mat-vec
		x_old = x;
		//x = (alpha*x_const) + (1-alpha)*(A_mat*x);
		x = (alpha*x_const) + (1-alpha)*(A_mat*x);
//		cout<<"intermed x = "<<x<<endl;
//		cout<<"sum = "<<getSumOfVectorElems(x)<<endl;
		//matrixValType sum_new = getSumOfVectorElems(x);
		
		//diff = abs(sum_new-sum_old);
		//sum_old = sum_new;

	}
	cout<<"\ntotal num iters "<<i<<endl; 

	return x;

}

void IncrementalPageRank(MyMatrix& A_new_mat, MyMatrix& delta_A_mat, VectorXd& x_const, VectorXd& x_0, double alpha, double threshold, int num_iters){

	//VectorXd 
}

struct Node{
	int id;
	vector<int> neighbours;
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

int main()
{
	MatrixXd m(2,2);
	m(0,0) = 3;
	m(1,0) = 2.5;
	m(0,1) = -1;
	m(1,1) = m(1,0) + m(0,1);
	std::cout << m << std::endl;

	//std::vector<double> vals = {22, 7, 3, 5, 14, 1, 17, 8};

	string filename = "edgelist.in";
	vector<Node> adj_list = CreateAdjacencyListFromEdgeList(filename);
	const int NUMROWS = adj_list.size();
	const int NUMCOLS = NUMROWS;
	const int NUM_ITERS = 300;

	std::vector<int> rows;// = {1, 2, 0, 2, 4, 2, 1, 4};
	std::vector<int> cols;// = {0, 0, 1, 1, 2, 3, 4, 4};

	for(int i=0;i<NUMROWS;i++){
		cout<<i<<": ";
		//if no outer link, assume it links to all - this is done in original PR as well!
		if(adj_list[i].neighbours.size() == 0){
			for(int j=0;j<NUMROWS;j++){
				adj_list[i].neighbours.push_back(j);
			}
		}
		for(int j=0;j<adj_list[i].neighbours.size();j++){
			rows.push_back(i);cols.push_back(adj_list[i].neighbours[j]);
			cout<<adj_list[i].neighbours[j]<<" ";
		}


		cout<<endl;
	}

	




	int nnz = rows.size();
/*	
	SparseMatrix<double> mat(NUMROWS, NUMCOLS);
	mat.reserve(VectorXi::Constant(NUMCOLS,6));

	for(int i=0;i<vals.size();i++){
		mat.insert(rows[i],cols[i]) = vals[i];
	}
*/
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

	PrintMatrix(mat);

	//d because double - need to change this if changing type
	VectorXd x(NUMROWS), x_init(NUMROWS);

	//initializing pagerank vector
	for(int i=0;i<NUMROWS;i++){
		x_init(i) = 1.0/double(NUMROWS);
	}

	x = x_init;
	//std::cout << "Here is the vector x:\n" << x << std::endl;
	std::cout << "num rows and cols\n" << x.rows() << " "<<x.cols()<<std::endl;
	std::cout << "Here is the vector x_init:\n" << x_init << std::endl;


	//SparseMatrix<matrixValType>& A_mat, VectorXd& x_const, VectorXd& x_0, double alpha, double threshold, int num_iters
	double threshold =0.0000000001;
	cout<<"thresh:" <<threshold<<endl;
	MyMatrix mat_transpose = MyMatrix(mat.transpose());
	VectorXd result = PageRank(mat_transpose, x_init,x, ALPHA, threshold, NUM_ITERS);

	std::cout << "Here is the vector result:\n" << result << std::endl;
	//x = mat*x;
	
	std::cout << "Here is the vector x:\n" << x << std::endl;
}