#ifndef DENSE_VECTOR_H
#define DENSE_VECTOR_H

#include "common.h"

class DenseVector{
	//vector containing all indices with non-zero vals (this vector is unordered, i.e. indices are not stored in increasing order)
	vector<valtype> vals;

public:
	DenseVector(int num_elements){
		vals.resize(num_elements,0);
	}

	DenseVector(vector<valtype> dense_vect){
		vals = dense_vect;
	}

	size_t size(){
		return vals.size();
	}

	//TODO: double check overloading format - const and all?
	valtype operator [](uint i) const    {
		if (i<0 || i>= vals.size()){
			cout<<"ERROR in ACCESSING DenseVector element!\n";
		}
		return vals[i];
	}
    valtype& operator [](uint i) {
		if (i<0 || i>= vals.size()){
			cout<<"ERROR in ACCESSING DenseVector element!\n";
		}return vals[i];}


};

#endif //DENSE_VECTOR_H