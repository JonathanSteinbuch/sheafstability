//============================================================================
// Name        : SparsePolyMatrix.hpp
// Author      : Jonathan Steinbuch
//============================================================================


#ifndef SRC_MATRICES_SPARSEPOLYMATRIX_HPP_
#define SRC_MATRICES_SPARSEPOLYMATRIX_HPP_

#include <iomanip>
#include <vector>

#include "../ring/Poly.hpp"
#include "DensePolyMatrix.hpp"
#include "SparseScalarMatrix.hpp"


using namespace std;

template <typename _Scalar>
class Poly;

template <typename _Scalar>
class DensePolyMatrix;


template<typename _Scalar>
using spPRow = boost::container::flat_map<unsigned int, Poly<_Scalar> >;

using spPCol = boost::container::flat_set<unsigned int>;

template <typename _Scalar>
class SparsePolyMatrix{
	PolyRing<_Scalar>* Ring;

	vector <int> indegrees;
	vector <int> outdegrees;

	vector <string> inlabels;
	vector <string> outlabels;

	vector<spPRow<_Scalar> > rows;
	vector<spPCol > cols;

	bool computeDegrees(const vector <int> &outdegs);

	void computeCols(unsigned int numcols);
public:
//Constructors:
	SparsePolyMatrix(): Ring(NULL) {}

	SparsePolyMatrix(PolyRing<_Scalar>* baseRing,const vector<spPRow<_Scalar> > &rows,const vector<spPCol > &cols,const vector <int> &outdegs,const vector <int> &indegs);

	SparsePolyMatrix(PolyRing<_Scalar>* baseRing,const vector<spPRow<_Scalar> > &rows,const vector<spPCol > &cols,const vector <int> &outdegs);

	SparsePolyMatrix(PolyRing<_Scalar>* baseRing,const vector<spPRow<_Scalar> > &rows,unsigned int numcols,const vector <int> &outdegs);

	SparsePolyMatrix(PolyRing<_Scalar>* baseRing,const vector<spPRow<_Scalar> > &rows,unsigned int numcols,const vector <int> &outdegs,const vector <int> &indegs);


	SparsePolyMatrix(PolyRing<_Scalar>* baseRing,const SparseScalarMatrix<_Scalar>& inmatrix,const unsigned int sx,const unsigned int sy, const vector <int> &outdegs, const int degree, bool homogeneous = true);

	SparsePolyMatrix(const DensePolyMatrix<_Scalar> &old);


//Matrix operations:
	SparseScalarMatrix<_Scalar> degreeMatrix(const int degree,bool homogeneous = true);

	SparsePolyMatrix<_Scalar> symMatrix(const unsigned int power) const;

	SparsePolyMatrix<_Scalar> extMatrix(const unsigned int power) const;

	SparsePolyMatrix<_Scalar> powerMatrix(const unsigned int power) const;
//Computations:
	SparsePolyMatrix<_Scalar> operator*(const SparsePolyMatrix<_Scalar>& rhs) const;

	SparsePolyMatrix<_Scalar> kernel(unsigned int degreeOfInterest, bool homogeneous = true, const program_options & opt = noOptions);
//Getters and Setters:
	void setLabels(const vector<string> &inlabels,const vector<string> &outlabels) {
		this->inlabels = inlabels;
		this->outlabels = outlabels;
	}

	int indegree();

	int outdegree();

	const unsigned int getCols() const {
		return cols.size();
	}

	const unsigned int getRows() const {
		return rows.size();
	}

	const vector<int>& getIndegrees() const {
		return indegrees;
	}

	const vector<int>& getOutdegrees() const {
		return outdegrees;
	}

	PolyRing<_Scalar>* getRing() const
	{
		return Ring;
	}

	const Poly<_Scalar> randomreadaccess(unsigned int r, unsigned int c) const;

//Output:
	template <typename tripScalar>
	void degreeTriplets(const int degree, vector<rcvTriplet<tripScalar> > & tripletList, unsigned int& trows, unsigned int& tcols, bool homogeneous = true) const;

	template <typename Scalar>
	friend std::ostream &operator<<(std::ostream &os, SparsePolyMatrix<Scalar> const &m);

	void printForMacaulay2(std::ostream &os);
	void printForLatex(std::ostream &os);



};



#endif /* SRC_MATRICES_SPARSEPOLYMATRIX_HPP_ */
