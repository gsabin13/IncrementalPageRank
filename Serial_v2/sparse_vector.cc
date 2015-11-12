#include "sparse_vector.h"

int SparseVector::getNNZ(){
	return indexes_sparse.size();
}

bool SparseVector::isIndexNonZero(int i){
	return indexes_dense[i];
}

DenseVector SparseVector::getDenseVectorForm(){
	DenseVector dense_vect(vals_dense);
	return dense_vect;
}

void SparseVector::SparseVectorCopyFrom(SparseVector src){
	int i;
	for(i=0;i<src.indexes_sparse.size() && i < (this->indexes_sparse.size());i++){
		this->indexes_sparse[i]=src.indexes_sparse[i];
		this->vals_sparse[i] = src.vals_sparse[i];
		this->indexes_dense[src.indexes_sparse[i]] = true;
		this->vals_dense[src.indexes_sparse[i]] = src.vals_sparse[i];
	}
	for(;i<src.indexes_sparse.size();i++){
		this->indexes_sparse.push_back(src.indexes_sparse[i]);
		this->vals_sparse.push_back(src.vals_sparse[i]);
		this->indexes_dense[src.indexes_sparse[i]] = true;
		this->vals_dense[src.indexes_sparse[i]] = src.vals_sparse[i];
	}
}

//inserts index into indexes_dense and indexes_sparse if not already non-zero. Does *not* insert value corresponding to index
void SparseVector::insertIndexCorrespondingToNonZeroVal(int i){
	
	//cout<<"inserting "<<i<<endl;
	if(!(this->isIndexNonZero(i))){
		//cout<<"actually inserting "<<i<<endl;
		indexes_dense[i] = true;
		indexes_sparse.push_back(i);
	}
}

//assumes vals_dense, index_dense, index_sparse have been updated, and now updates vals_sparse based on those 3 arrays
void SparseVector::updateValsSparseBasedOnValsDenseAndIndexesSparse(){
	uint i;
	//replace existing vals_sparse elements with new vals
	for(i=0;i<vals_sparse.size();i++){
		vals_sparse[i] = vals_dense[indexes_sparse[i]];
	}
	//add newly added vals_sparse elements
	for(;i<indexes_sparse.size();i++){
		vals_sparse.push_back(vals_dense[indexes_sparse[i]]);
	}
}

//method that adds rhs_vect to the existing object
void SparseVector::inPlaceAdd(SparseVector rhs_vect){
	//SparseVector result(rhs_vect);
	for(uint i=0;i<rhs_vect.indexes_sparse.size();i++){
		//update indexes_sparse and dense arrays
		this->insertIndexCorrespondingToNonZeroVal(rhs_vect.indexes_sparse[i]);
		vals_dense[rhs_vect.indexes_sparse[i]] += rhs_vect.vals_sparse[i];
	}
	
	updateValsSparseBasedOnValsDenseAndIndexesSparse();
	
}

void SparseVector::inPlaceScalarMultiply(valtype scalar){
	for(uint i=0;i<indexes_sparse.size();i++){
		vals_sparse[i]*=scalar;
		vals_dense[indexes_sparse[i]]=vals_sparse[i];
	}
}

//sets all non-zero values to 0.
//TODO: check perf of erase and resize vs fill all vs filling out each non-zero individually?
void SparseVector::clearDenseVals(){
	for(uint i =0; i<indexes_sparse.size();i++){
		vals_dense[indexes_sparse[i]] = 0;
	}
}
