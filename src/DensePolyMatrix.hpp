//============================================================================
// Name        : DensePolyMatrix.hpp
// Author      : Jonathan Steinbuch
//============================================================================


#ifndef SRC_DENSEPOLYMATRIX_HPP_
#define SRC_DENSEPOLYMATRIX_HPP_

#include <iomanip>
#include <vector>

#include "Poly.hpp"
#include "SparseScalarMatrix.hpp"
#include "SparsePolyMatrix.hpp"

const unsigned int ERROR = UINT_MAX;

template <typename _Scalar>
class Poly;

template <typename _Scalar>
class SparsePolyMatrix;

using namespace std;

//Implements a densely represented matrix with entries in the given Ring.
template <typename _Scalar>
class DensePolyMatrix{
	PolyRing<_Scalar>* Ring; //base ring

	unsigned int rows, cols; //dimensions of the matrix

	vector <int> indegrees; //Twists of the input
	vector <int> outdegrees; //Twists of the output

	vector <string> inlabels; //Labels for the input
	vector <string> outlabels; //Labels for the output

	vector<vector<Poly<_Scalar> > > entries; //The main matrix data

	bool computeDegrees(const vector <int> &outdegs);

	template <typename tripScalar>
	vector<rcvTriplet<tripScalar> > degreeTriplets(const int degree, unsigned int& trows, unsigned int& tcols, bool homogeneous = true);

public:
//Constructors:
	DensePolyMatrix() {}

	DensePolyMatrix(PolyRing<_Scalar>* baseRing,const vector< vector <Poly<_Scalar> > > &entries);

	DensePolyMatrix(PolyRing<_Scalar>* baseRing,const vector< vector <string > > &entries);

	DensePolyMatrix(PolyRing<_Scalar>* baseRing,const unsigned int sx, const unsigned int sy,const vector< vector <Poly<_Scalar> > > &entries);

	DensePolyMatrix(PolyRing<_Scalar>* baseRing,const unsigned int sx, const unsigned int sy,const vector< vector <Poly<_Scalar> > > &entries, const vector <int> &outdegs);

	DensePolyMatrix(PolyRing<_Scalar>* baseRing,const unsigned int sx, const unsigned int sy,const vector< vector <Poly<_Scalar> > > &entries, const vector <int> &outdegs, const vector <int> &indegs);

	template <unsigned int sx, unsigned int sy>
	DensePolyMatrix(PolyRing<_Scalar>* baseRing,const Poly<_Scalar> (&entries)[sx][sy]):
			Ring(baseRing), rows(sx), cols(sy)
	{
		for (unsigned int i = 0; i < sx; i++)
		{
			vector<Poly<_Scalar> > row;
			row.clear();
			for (unsigned int j = 0; j < sy; j++)
			{
				row.push_back(entries[i][j]);
			}
			this->entries.push_back(row);
		}
		computeDegrees(vector<int>());
	}

	DensePolyMatrix(PolyRing<_Scalar>* baseRing,const SparseScalarMatrix<_Scalar>& inmatrix,const unsigned int sx,const unsigned int sy, const vector <int> &outdegs, const int degree, bool homogeneous = true);

	DensePolyMatrix(const SparsePolyMatrix<_Scalar> &old);

//Matrix operations:
	SparseScalarMatrix<_Scalar> degreeMatrix(const int degree,bool homogeneous = true);

	template <typename Scalar>
	friend DensePolyMatrix<Scalar> directsum(const vector<DensePolyMatrix<Scalar> > &summands);
	template <typename Scalar>
	friend DensePolyMatrix<Scalar> directsum(const vector<DensePolyMatrix<Scalar>* > &summands);

	template <typename Scalar>
	friend DensePolyMatrix<Scalar> sum(const vector<const DensePolyMatrix<Scalar>* > &summands);

	DensePolyMatrix<_Scalar> ntimes(const unsigned int power) const;

	static DensePolyMatrix<_Scalar> irrelevantPower(PolyRing<_Scalar> &mB, int degree);

	DensePolyMatrix<_Scalar> symMatrix(const unsigned int power) const;

	DensePolyMatrix<_Scalar> extMatrix(const unsigned int power) const;

	DensePolyMatrix<_Scalar> powerMatrix(const unsigned int power) const;

//Computations:
	DensePolyMatrix<_Scalar> operator*(const DensePolyMatrix<_Scalar>& rhs) const;

	DensePolyMatrix<_Scalar> operator%(const DensePolyMatrix<_Scalar>& rhs) const;

	DensePolyMatrix<_Scalar> kernel(unsigned int degreeOfInterest, bool homogeneous = true, const program_options & opt = noOptions);

	DensePolyMatrix<_Scalar> jacobiMatrix();

	static DensePolyMatrix<_Scalar> jacobiTaylor(PolyRing<_Scalar>* mB, const Poly<_Scalar>& poly, unsigned int n);

//Getters and Setters:
	void setLabels(const vector<string> &inlabels,const vector<string> &outlabels) {
		this->inlabels = inlabels;
		this->outlabels = outlabels;
	}

	int indegree();

	int outdegree();

	const unsigned int getCols() const {
		return cols;
	}

	const unsigned int getRows() const {
		return rows;
	}

	const vector<int>& getIndegrees() const {
		return indegrees;
	}

	const vector<int>& getOutdegrees() const {
		return outdegrees;
	}

	const Poly<_Scalar> getEntry(unsigned int row, unsigned int col) const {
		return entries[row][col];
	}

	const vector<Poly<_Scalar> > getRowEntries(unsigned int row) const {
		return entries[row];
	}

	const vector<Poly<_Scalar> > computeMinors(unsigned int n) const;

	PolyRing<_Scalar>* getRing() const
	{
		return Ring;
	}

//Output:
	template <typename Scalar>
	friend std::ostream &operator<<(std::ostream &os, DensePolyMatrix<Scalar> const &m);

	void printAsOperator(std::ostream &os);

	void printForMacaulay2(std::ostream &os);


	friend class SparsePolyMatrix<_Scalar>;
};



#endif /* SRC_DENSEPOLYMATRIX_HPP_ */
