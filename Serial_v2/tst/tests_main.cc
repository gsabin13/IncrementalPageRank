#include <bits/stdc++.h>
#include "../dense_vector.h"
#include "../sparse_vector.h"
#include "../sparse_matrix_csc.h"

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

	s1.inPlaceAdd(s2);
	s1.printDenseRepresentation();
	s1.printSparseRepresentation();

}

void testSpMMSubtract(){

	int c1[] = {0, 2, 3, 6};
	int r1[] = {0, 2, 2, 0, 1, 2};
	int v1[] = {1, 4, 5, 2, 3, 6};

	int c2[] = {0, 1, 3, 4};
	int r2[] = {1, 1, 2, 0};
	int v2[] = {3, 1, 6, 1};
	vector<int> cols1(c1, c1+4);
	vector<int> cols2 (c2, c2+4);

	vector<int> rows1(r1, r1+6);
	vector<int> rows2(r2, r2+4);

	vector<valtype> vals1(v1, v1+6);
	vector<valtype> vals2(v2, v2+4);

	SparseMatrixCSC mat1(cols1, rows1, vals1);
	SparseMatrixCSC mat2(cols2, rows2, vals2);
	mat1.PrintMatrix();
	mat2.PrintMatrix();

	SparseMatrixCSC mat3 = mat1.subtract(mat2);
	mat3.PrintMatrix();

	int c4[] = {0, 1,1,1};
	int r4[] = {1};
	int v4[] = {3};
	vector<int> cols4 (c4, c4+4);
	vector<int> rows4(r4, r4+1);
	vector<valtype> vals4(v4, v4+1);
	SparseMatrixCSC mat4(cols4, rows4, vals4);
	mat3 = mat1.subtract(mat4);
	mat3.PrintMatrix();



}

void testSpMSpV(){
	int c1[] = {0, 2, 3, 6};
	int r1[] = {0, 2, 2, 0, 1, 2};
	int v1[] = {1, 4, 5, 2, 3, 6};

	int c2[] = {0, 1, 2, 4};
	int r2[] = {0, 2, 0, 2};
	int v2[] = {3, 1, 1, 6};
	vector<int> cols1(c1, c1+4);
	vector<int> cols2 (c2, c2+4);

	vector<int> rows1(r1, r1+6);
	vector<int> rows2(r2, r2+4);

	vector<valtype> vals1(v1, v1+6);
	vector<valtype> vals2(v2, v2+4);

	SparseMatrixCSC mat1(cols1, rows1, vals1);
	SparseMatrixCSC mat2(cols2, rows2, vals2);
	mat1.PrintMatrix();
	
	vector<valtype> x(3,0);
	
	x[0] = 1; x[2] = 9;
	DenseVector d1(x);
	SparseVector s1(d1);
	SparseVector s2(s1);

	mat1.inPlaceMultiplyWithSparseVector(s1);

	s1.printDenseRepresentation();
	s1.printSparseRepresentation();
	mat2.PrintMatrix();
	mat2.inPlaceMultiplyWithSparseVector(s2);
	s2.printDenseRepresentation();
	s2.printSparseRepresentation();
	mat2.inPlaceMultiplyWithSparseVector(s2);
	s2.printDenseRepresentation();
	s2.printSparseRepresentation();

	s2.inPlaceAdd(s1);																																				
	s2.printDenseRepresentation();
	s2.printSparseRepresentation();

	mat2.inPlaceMultiplyWithSparseVector(s2);
	s2.printDenseRepresentation();
	s2.printSparseRepresentation();
}

void testSpMV(){

	int c1[] = {0, 2, 3, 6};
	int r1[] = {0, 2, 2, 0, 1, 2};
	int v1[] = {1, 4, 5, 2, 3, 6};

	int c2[] = {0, 1, 2, 4};
	int r2[] = {0, 2, 0, 2};
	int v2[] = {3, 1, 1, 6};
	vector<int> cols1(c1, c1+4);
	vector<int> cols2 (c2, c2+4);

	vector<int> rows1(r1, r1+6);
	vector<int> rows2(r2, r2+4);

	vector<valtype> vals1(v1, v1+6);
	vector<valtype> vals2(v2, v2+4);

	SparseMatrixCSC mat1(cols1, rows1, vals1);
	SparseMatrixCSC mat2(cols2, rows2, vals2);
	mat1.PrintMatrix();

	vector<valtype> x(3,2);
	DenseVector d1(x);
	d1.Print();

	DenseVector d2 = mat1.multiplyWithDenseVector(d1);
	d2.Print();

}

int main(){
	
	testDenseVector();

	testSparseVector();

	testSpMMSubtract();

	testSpMSpV();

	testSpMV();

	cout<<"tests ended\n";
	
	return 1;
}