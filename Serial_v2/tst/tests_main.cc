#include <bits/stdc++.h>
#include "../dense_vector.h"
#include "../sparse_vector.h"

void testDenseVector(){
	vector<valtype> x(10,1);
	DenseVector d(x);
	if(d.size()!=10){
		cout<<"testDenseVector TEST FAILED!"<<endl;	
	}
	
	for(int i=0;i<d.size();i++){
		if(d[i] != x[i]){
			cout<<"testDenseVector TEST FAILED!"<<endl;	
		}
	}
}

void testSparseVector(){
	vector<valtype> x(10,0);
	
	x[0] = 1.2; x[5] = 9.2;
	DenseVector d1(x);
	SparseVector s1(d1);

	if(s1.isIndexNonZero(3) || !s1.isIndexNonZero(5)){
		cout<<"testSparseVector TEST FAILED!"<<endl;	
	}

	x[0] = 0; x[2]=1.3; x[5]=2.0;
	DenseVector d2(x);
	SparseVector s2(d2);

	s1.incrementalAdd(s2);
	s1.printDenseRepresentation();
	s1.printSparseRepresentation();

}

int main(){
	
	testDenseVector();

	testSparseVector();



	cout<<"tests ended\n";
	
	return 1;
}