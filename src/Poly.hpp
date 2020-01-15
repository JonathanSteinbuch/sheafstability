//============================================================================
// Name        : Poly.hpp
// Author      : Jonathan Steinbuch
//============================================================================


#ifndef SRC_POLY_HPP_
#define SRC_POLY_HPP_

#include <gmpxx.h>

template<typename _Scalar>
class PolyRing;

#include "SparseScalarMatrix.hpp"
#include "PolyRing.hpp"
#include "Helpers.hpp"

template<typename Scalar>
using polyContainer = map<mBmonomial, Scalar>;

//This class implements a polynomial
template<typename _Scalar>
class Poly
{
private:
	PolyRing<_Scalar>* Ring; //The ring of which the polynomial is an element

	polyContainer<_Scalar> poly; //The actual polynomial data

	Poly(PolyRing<_Scalar>* baseRing, string::const_iterator& it, const string& input, bool bracket);

public:
//Constructors:
	Poly(const Poly<_Scalar>& old) :
			Ring(old.Ring), poly(old.poly)
	{
	}

	Poly(PolyRing<_Scalar>* baseRing, const mBmonomial index, const _Scalar coefficient) :
			Ring(baseRing)
	{
		poly.clear();
		poly.insert(pair<mBmonomial, _Scalar>(index, coefficient));
	}

	Poly(PolyRing<_Scalar>* baseRing, const _Scalar coefficient) :
			Ring(baseRing)
	{
		poly.clear();
		poly.insert(pair<mBmonomial, _Scalar>(0, coefficient));
	}

	Poly(PolyRing<_Scalar>* baseRing) :
			Ring(baseRing)
	{
	}

	Poly(PolyRing<_Scalar>* baseRing, const string& input);

	static Poly<_Scalar> variable(PolyRing<_Scalar>* baseRing, const unsigned int i);

	static Poly<_Scalar> generic(PolyRing<_Scalar>* baseRing, const unsigned int degree);

	static Poly<_Scalar> zero(PolyRing<_Scalar>* baseRing)
	{
		return Poly<_Scalar>(baseRing);
	}

//Order and Monomials:
	const unsigned int supportSize() const{
		return poly.size();
	}

	const unsigned int degree() const;

	const pair<mBmonomial, _Scalar> leadingTerm() const;

	const pair<mBmonomial, _Scalar> trailingTerm() const;

	const mBmonomial leadingMonomial() const;

	const mBmonomial trailingMonomial() const;

	const _Scalar leadingCoefficient() const;

	const _Scalar trailingCoefficient() const;

	bool popLeadingTerm();

	bool isHomogeneous() const;

	const bool isNull() const;

	void homogenize(unsigned int variable);

	Poly<_Scalar>& trim();

	void simplify(); //divide out common integer factors
	bool reduce(Poly<_Scalar>& reducer);

//Arithmetic:
	Poly<_Scalar>& monMultEq(const mBmonomial monomial);

	Poly<_Scalar> monMult(const mBmonomial monomial) const;

	Poly<_Scalar> monAdd(const mBmonomial rhs) const;

	Poly<_Scalar>& operator*=(const pair<mBmonomial, _Scalar>& rhs);
	Poly<_Scalar> operator*(const pair<mBmonomial, _Scalar>& rhs) const;

	Poly<_Scalar>& operator*=(const _Scalar &rhs);
	Poly<_Scalar>& operator/=(const _Scalar &rhs);

	Poly<_Scalar> operator*(const _Scalar &rhs) const;

	template<typename Scalar>
	friend Poly<Scalar> operator*(const Scalar &lhs,Poly<Scalar> rhs);

	Poly<_Scalar>& operator+=(const pair<mBmonomial, _Scalar>& rhs);

	Poly<_Scalar>& operator+=(const Poly<_Scalar>& rhs);
	Poly<_Scalar> operator+(const Poly<_Scalar>& rhs) const;

	Poly<_Scalar>& operator-=(const Poly<_Scalar>& rhs);

	bool operator==(const Poly<_Scalar>& rhs);

	Poly<_Scalar> operator-(const Poly<_Scalar>& rhs) const;

	Poly<_Scalar>& operator*=(Poly<_Scalar> rhs);

	Poly<_Scalar> operator*(const Poly<_Scalar>& rhs) const;

	Poly<_Scalar> taylorDerive(mBmonomial monomial) const;

	Poly<_Scalar>& operator^=(const unsigned int exponent);

	Poly<_Scalar> operator^(const unsigned int exponent) const;

//Output and Getters:
	template<typename tripScalar>
	vector<rcvTriplet<tripScalar> >& appendMatrix(vector<rcvTriplet<tripScalar> >& ret, const int outdegree, const int indegree, unsigned int rows, unsigned int cols, const unsigned int rowoffset = 0, const unsigned int coloffset = 0, bool homogeneous = true) const;

	vector<rcvTriplet<_Scalar> > getMatrix(const int outdegree, const int indegree, unsigned int rows, unsigned int cols, const unsigned int rowoffset = 0, const unsigned int coloffset = 0);

	template<typename Scalar>
	friend std::ostream &operator<<(std::ostream &os, const Poly<Scalar> &m);

	void print(std::ostream &os, bool longform=false) const;

	PolyRing<_Scalar>* getRing() const
	{
		return Ring;
	}

	//Ring change. Only possible if base rings coincide.
	bool changeRing(PolyRing<_Scalar>* newRing)
	{
		if(newRing->getBaseRing() == Ring->getBaseRing())
		{
			Ring = newRing;
			return true;
		} else {
			return false;
		}
	}
};

#endif /* SRC_POLY_HPP_ */
