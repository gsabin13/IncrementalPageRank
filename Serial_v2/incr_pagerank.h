#ifndef INCR_PAGERANK_H
#define INCR_PAGERANK_H

#include "common.h"
#include "sparse_vector.h"
#include "dense_vector.h"
#include "sparse_matrix_csc.h"
#include "incr_pagerank.h"


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
#endif //INCR_PAGERANK_H