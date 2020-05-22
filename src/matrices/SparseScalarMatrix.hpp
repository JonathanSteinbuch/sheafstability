//============================================================================
// Name        : SparseScalarMatrix.hpp
// Author      : Jonathan Steinbuch
//============================================================================

#ifndef SRC_MATRICES_SPARSESCALARMATRIX_HPP_
#define SRC_MATRICES_SPARSESCALARMATRIX_HPP_
#include <cassert>
#include <vector>
#include <numeric>
#include <map>
#include <set>
#include <iostream>

#include <gmpxx.h>

#include <boost/container/flat_set.hpp>
#include <boost/container/flat_map.hpp>

#include "../helpers/Helpers.hpp"

using namespace std;

template<typename _Scalar>
class SparseScalarMatrix;

template<typename Scalar>
using spRow = boost::container::flat_map<unsigned int, Scalar>;

using spCol = boost::container::flat_set<unsigned int>;

/*template<typename Scalar>
using spRow = map<unsigned int, Scalar>;

using spCol = set<unsigned int>;*/

class Permutation
{

public:

	vector<unsigned int> _image;
	vector<unsigned int> _preimage;

	Permutation()
	{
	}

	Permutation(unsigned int n)
	{
		for (unsigned int i = 0; i < n; i++)
		{
			_image.push_back(i);
			_preimage.push_back(i);
		}
	}

	Permutation(unsigned int n,const vector<unsigned int>& order, bool isImage = true)
	{
		if (isImage)
		{
			_preimage = vector<unsigned int>(n, UINT_MAX);
			for (unsigned int i = 0; i < n; i++)
			{
				unsigned int im_i = i < order.size() ? order[i] : i;
				_image.push_back(im_i);
				_preimage[im_i] = i;
			}
		}
		else
		{
			_image = vector<unsigned int>(order.size(), UINT_MAX);

			for (unsigned int i = 0; i < n; i++)
			{
				unsigned int preim_i = i < order.size() ? order[i] : i;
				_preimage.push_back(preim_i);
				_image[preim_i] = i;
			}
		}
	}

	Permutation(vector<unsigned int>& order, bool isImage = true) :
			Permutation(order.size(), order, isImage)
	{

	}

	const unsigned int size() const
	{
		return _image.size();
	}

	Permutation& operator*=(const Permutation& rhs)
	{
		*this = (*this) * rhs;
		return *this;
	}

	bool checkintegrity()
	{
		for (unsigned int i = 0; i < _image.size(); i++)
		{
			if (_preimage[_image[i]] != i)
			{
				return false;
			}
		}

		for (unsigned int i = 0; i < _preimage.size(); i++)
		{
			if (_image[_preimage[i]] != i)
			{
				return false;
			}
		}
		return true;
	}

	friend Permutation operator*(const Permutation& lhs, const Permutation& rhs)
	{
		vector<unsigned int> order;
		assert(lhs.size() == rhs.size());
		for (unsigned int i = 0; i < lhs.size(); i++)
		{
			unsigned int x = lhs.image(rhs.image(i));
			order.push_back(x);
		}
		Permutation nPerm(order);
		return nPerm;
	}

	void swap(unsigned int a, unsigned int b)
	{
		unsigned int pa = _preimage[a];
		unsigned int pb = _preimage[b];
		_image[pa] = b;
		_image[pb] = a;
		_preimage[a] = pb;
		_preimage[b] = pa;
	}

	void swap(vector<unsigned int> sequence)
	{
		for (unsigned int i = 0; i < sequence.size() - 1; i++)
		{
			_image[_preimage[sequence[i]]] = sequence[i + 1];
		}
		_image[_preimage[sequence[sequence.size() - 1]]] = sequence[0];

		unsigned int ptemp = _preimage[sequence[sequence.size() - 1]];
		for (unsigned int i = sequence.size() - 1; i >= 1; i--)
		{
			_preimage[sequence[i]] = _preimage[sequence[i - 1]];
		}
		_preimage[sequence[0]] = ptemp;
	}

	const unsigned int image(unsigned int i) const
	{
		return _image[i];
	}

	const unsigned int preimage(unsigned int i) const
	{
		return _preimage[i];
	}

	template<typename Scalar>
	SparseScalarMatrix<Scalar> matrix(bool invert = false)
	{
		vector<spRow<Scalar> > rows;
		vector<spCol > cols;

		assert(_preimage.size() == _image.size());
		unsigned int s = _image.size();

		for (unsigned int i = 0; i < s; i++)
		{
			spRow<Scalar> row;
			unsigned int pos = invert ? _image[i] : _preimage[i];
			row.insert(pair<unsigned int, Scalar>(pos, 1));
			rows.push_back(row);
		}

		for (unsigned int i = 0; i < s; i++)
		{
			spCol col;
			unsigned int pos = invert ? _preimage[i] : _image[i];
			col.insert(pos);
			cols.push_back(col);
		}
		return SparseScalarMatrix<Scalar>(rows, cols);
	}
};

/*template<typename Scalar>
 Scalar gcd(Scalar u, Scalar v)
 {
 while (v != 0)
 {
 Scalar r = u % v;
 u = v;
 v = r;
 }
 return u;
 }*/
const mpq_class gcd(const mpq_class &u, const mpq_class &v);
const mpq_class lcm(const mpq_class &u, const mpq_class &v);

template<typename Scalar>
struct rcvTriplet{
	unsigned int row;
	unsigned int col;
	Scalar value;

	rcvTriplet(unsigned int row, unsigned int col, Scalar value) : row(row), col(col), value(value)
	{

	}
};

template<typename Scalar>
class SparsePolyMatrix;

template<typename Scalar>
class SparseScalarMatrix
{
	Permutation rowPerm;
	Permutation colPerm;

	vector<spRow<Scalar> > rows;
	vector<spCol > cols; // we only store positions of nonzeros in the columns

	spCol outrows;

	Scalar maxabs = 0;

//Basic Linear Algebra:
	Scalar rowproduct(const spRow<Scalar>& v1, const spRow<Scalar>& v2);

	void add(unsigned int r, unsigned int c,const Scalar &value);

	// moves the rows specified in sequence backwards
	void swapRows(const vector<unsigned int>& sequence);

	void swapRows(unsigned int a, unsigned int b);

	void swapCols(unsigned int a, unsigned int b);

	void symmSwaps(vector<unsigned int>& newpos, bool invert = false);

	void multiply(spRow<Scalar>& v, Scalar factor);

	void addrows(const unsigned int targetrow,const unsigned int otherrow,const Scalar &factortarget,const Scalar &factorother);
	void addrowsslow(const unsigned int targetrow,const unsigned int otherrow,const Scalar &factortarget,const Scalar &factorother);

	void addcols(unsigned int targetcol, unsigned int othercol, Scalar factortarget, Scalar factorother);

	void simplifyrow(unsigned int row);

	unsigned int rowcount(unsigned int row) const;

	unsigned int rowcount(unsigned int row,const pair<unsigned int, unsigned int>& range) const;

	SparseScalarMatrix<Scalar> colSelect(vector<unsigned int> selection);

	SparseScalarMatrix<Scalar> colUnSelect(unsigned int unselectionCol, unsigned int unselectionRow = INT_MAX);

	void fillcolumns();

	void reduce(unsigned int &defect, const program_options &opt = noOptions);

	void tarjan(unsigned int maxr, vector<pair<unsigned int, unsigned int> >& blocks);

	bool checkintegrity();


public:
//Constructors:
	SparseScalarMatrix(const vector<rcvTriplet<Scalar> > &triplets, unsigned int nrows, unsigned int ncols);

	SparseScalarMatrix(const SparsePolyMatrix<Scalar> &inputMatrix, const int degree, bool homogeneous = true);

	SparseScalarMatrix(const vector<spRow<Scalar> > &rows,const vector<spCol > &cols);

	SparseScalarMatrix()
	{
		outrows.clear();
	}

	SparseScalarMatrix(const vector<spRow<Scalar>>& vecs, unsigned int ncols, bool rowMajor = true);

//Arithmetic:
	SparseScalarMatrix<Scalar>& operator*=(const SparseScalarMatrix<Scalar>& rhs);

	bool operator==(const SparseScalarMatrix<Scalar>& rhs) const;

	SparseScalarMatrix<Scalar> operator*(const SparseScalarMatrix<Scalar>& rhs) const;

//Access:
	const Scalar randomreadaccess(unsigned int r, unsigned int c) const;

	size_t getRows() const
	{
		return rows.size();
	}

	size_t getCols() const
	{
		return cols.size();
	}

//Output:
	void printblock(unsigned int sx, unsigned int sy, unsigned int ex, unsigned int ey, std::ostream &os = cout);

	template<typename _Scalar>
	friend std::ostream &operator<<(std::ostream &os, SparseScalarMatrix<_Scalar> const &m);

//Linear Algebra:
	unsigned int analyze( const program_options & opt = noOptions);

	SparseScalarMatrix<Scalar> modulo(unsigned int lhscols);

	SparseScalarMatrix<Scalar> kernel( const program_options &opt = noOptions);
};

#endif /* SRC_MATRICES_SPARSESCALARMATRIX_HPP_ */
