#ifndef DENSE_VECTOR_H
#define DENSE_VECTOR_H

#include "common.h"

class DenseVector{
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

	void InPlaceScalarMultiply(valtype scalar){
		for(uint i=0;i<vals.size();i++){
			vals[i]*=scalar;
		}

	}

	void VectorCopy(DenseVector& src){
		if(src.size() != this->size()){
			cout<<"ERROR in VectorCopy!\n";
		}
		for(int i=0;i<src.size();i++){
			vals[i] = src[i];
		}

	}

	void DenseVectorSum(DenseVector left, DenseVector right){
		if(left.size() != this->size() || right.size() != this->size()){
			cout<<"ERROR in VectorCopy!\n";
		}
		for(int i=0;i<this->size();i++){
			vals[i] = left[i] + right[i];
		}
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
		}
		return vals[i];
	}

	void Print(){
		cout<<"Dense vector: ";
		for(int i=0;i<vals.size();i++){
			cout<<vals[i]<<" ";
		}
		cout<<endl;
	}



};

#endif //DENSE_VECTOR_H