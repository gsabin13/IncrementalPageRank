#ifndef SPARSE_VECTOR_H
#define SPARSE_VECTOR_H

#include "common.h"
#include "dense_vector.h"

class SparseVector{

public:
	//vector containing all indices with non-zero vals (this vector is unordered, i.e. indices are not stored in increasing order)
	vector<int> indexes_sparse;

	//vector containing all corresponding values.
	vector<valtype> vals_sparse;


	//(ordered) dense representation of indices. if indexes_dense[i] == true, then i is contained in indexes_sparse
	vector<bool> indexes_dense;

	//Maybe this should be a DenseVector?
	//(ordered) dense representation of vals. if indexes_dense[i] == true, then vals_dense[i] is contained in vals_sparse. else vals_dense[i] = 0
	vector<valtype> vals_dense;

	SparseVector(int num_elements){
		indexes_dense.resize(num_elements,false);
		vals_dense.resize(num_elements,0);
	}

	SparseVector(DenseVector& dense_vect){
		indexes_dense.resize(dense_vect.size(),false);
		vals_dense.resize(dense_vect.size(),0);

		for(uint i=0;i<dense_vect.size();i++){	
			if(dense_vect[i] != 0){
				indexes_sparse.push_back(i);
				vals_sparse.push_back(dense_vect[i]);

				indexes_dense[i] = true;
				vals_dense[i] = dense_vect[i];
			}
		}
	}

	void printDenseRepresentation(){
		for(uint i=0;i<vals_dense.size();i++){
			cout<<vals_dense[i]<<" ";
		}
		cout<<endl;
	}

	void printSparseRepresentation(){
		for(uint i=0;i<vals_sparse.size();i++){
			cout<<indexes_sparse[i]<<" "<<vals_sparse[i]<<endl;
		}
		cout<<endl;
	}

	int getNNZ();

	bool isIndexNonZero(int i);

	//inserts index if not already non-zero. Does *not* insert value corresponding to index
	void insertIndexCorrespondingToNonZeroVal(int i);

	//assumes vals_dense, index_dense, index_sparse have been updated, and now updates vals_sparse based on those 3 arrays
	void updateValsSparseBasedOnValsDenseAndIndexesSparse();

	//method that adds rhs_vect to the existing object
	void incrementalAdd(SparseVector rhs_vect);

	//sets all non-zero values to 0 in the vals_dense array.
	//TODO: check perf of erase and resize vs fill all vs filling out each non-zero individually? currently doing the latter
	void clearDenseVals();
};

#endif //SPARSE_VECTOR_H