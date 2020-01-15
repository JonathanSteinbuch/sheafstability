//============================================================================
// Name        : PolyRing.hpp
// Author      : Jonathan Steinbuch
//============================================================================


#ifndef SRC_POLYRING_HPP_
#define SRC_POLYRING_HPP_

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <array>

#include <gmpxx.h>

#include "DensePolyMatrix.hpp"
#include "Helpers.hpp"
#include "Poly.hpp"


using namespace std;


template <typename _Scalar>
class Poly;

//This is a struct to store all data relating monomial IDs to the actual monomials and their degrees
struct monlookups{

	lutable lu;
	unsigned int maxdegree; //maximal degree to which we computed. Will change dynamically
	unsigned int vars; //number of variables
	vector<string> printvar; //print strings for the variables
	const unsigned int order; //deglex = 0 or degrevlex = 1;
	vector<vector<unsigned int> > list; //exponent-list for each monomial
	unsigned int realvars; //in the case of helper variables to compute homogeneously for non-homogeneous polynomials
	vector<int> weights; //weights of the variables; standard is 1
	vector<vector <unsigned int> > monomOffset; //internal offsets to compute the monomial ID from its exponent-list
	vector<unsigned int> monN;//monN[N] contains the number of monomials up to degree N-1

	monlookups(const unsigned int vars, const vector<string> (&printvar), const unsigned int order, int realvars = -1, const vector<int> &weights = {}, bool invert=false);
};

//This class implements (a quotient of) a polynomial ring.
template <typename _Scalar>
class PolyRing {

	PolyRing* BaseRing; //Base Ring. Make sure that the actual base ring has itself as base ring.

	monlookups* mons; //The monomial ID structure.

	vector<Poly<_Scalar> > gB; //The Gröbner Basis of the defining ideal
	vector<Poly<_Scalar> > lookup; //Lookup table for each monomial after applying reduction by the Gröbner Basis.
	vector<mBmonomial> basis; //Monomial Basis with respect to the Gröbner Basis
	vector<unsigned int> basisdegstart; //Where in the Basisvector each degree starts.

	void reCalcRelations(const unsigned int minDegree);

	void reSizeBasis();
	void reSize(const unsigned int neededDegree);
	void reSizeToIndex(const unsigned int neededIndex);


	unsigned int getDegree(const unsigned int index) const;

	unsigned int getIndex(const vector<unsigned int> &x);

	unsigned int getIndex(const unsigned int degree, const vector<unsigned int> &x);

	unsigned int getBasisIndexfromDegBasisIndex(const unsigned int degree, const unsigned int degBasisIndex) const;

	const Poly<_Scalar> &getLookup(const mBmonomial monomial);

	void monpoly(const int deg,const unsigned int vars, Poly<mpq_class> &retVal,const Poly<mpq_class> &n); //Helper for Hilbert Polynomial

	Poly<_Scalar> groebnerS(const Poly<_Scalar> &p, const Poly<_Scalar> &q);

	void setRelations(const vector<Poly<_Scalar> > &groebnerBasis);
public:
	static const unsigned int mB_ERROR = UINT_MAX;
	const program_options* opt;

//Constructor:
	PolyRing(const unsigned int vars, const vector<string> (&printvar), const unsigned int order, int realvars = -1, const vector<int> &weights = {}, bool invert=false);

	PolyRing(PolyRing* BaseRing, const vector<Poly<_Scalar> > &relations);

	PolyRing(PolyRing* BaseRing, const vector<string> &relations);

	~PolyRing();

	void setOptions(const program_options& opt)
	{
		this->opt = &opt;
		if(BaseRing != this)
		{
			BaseRing->setOptions(opt);
		}
	}

//Getters:
	vector<unsigned int> getHilbert();

	unsigned int getDegree(const vector<unsigned int> &x) const;

	unsigned int getDegree(const mBmonomial monomial) const;

	mBmonomial getMonomial(const vector<unsigned int> &x);

	unsigned int getVar(string var, unsigned int pos = 0);

	mBmonomial getMonomial(const unsigned int degree, const vector<unsigned int> &x);

	mBmonomial getMonomialfromBasis(const unsigned int degree, const unsigned int degBasisIndex) const;

	unsigned int getDegBasisIndex(const mBmonomial monomial) const;

	const vector<unsigned int> &getExponents(const mBmonomial monomial);

	unsigned int getDegBasisIndex(const unsigned int degree, const mBmonomial monomial) const;

	const Poly<_Scalar> getLookup(const mBmonomial monomial,const _Scalar coefficient);

	unsigned int getBasisSize(const int degree);

	const unsigned int getMonN(unsigned int degreePlus1)
	{
		reSize(degreePlus1);
		return mons->monN[degreePlus1];
	}

	const mBmonomial& getBasisElem(unsigned int index)
	{
		reSizeToIndex(index);
		return basis[index];
	}

	const unsigned int getIndexOfBasisDegStart(unsigned int degree)
	{
		reSize(degree);
		return basisdegstart[degree];
	}

	const vector<string>& getPrintVar() const
	{
		return mons->printvar;
	}


	const string& getPrintVar(unsigned int i) const
	{
		static const string joker = "?";
		if(i >= mons->printvar.size())
		{
			return joker;
		} else {
			return mons->printvar[i];
		}
	}

	unsigned int getNumRealVars() const
	{
		return mons->realvars;
	}

	unsigned int getNumVars() const
	{
		return mons->vars;
	}

	const vector<Poly<_Scalar> >& getGB() const
	{
		return gB;
	}

	PolyRing* getBaseRing()
	{
		return BaseRing;
	}

	bool isSmooth();

//Monomial arithmetic:


	mBmonomial mongcd(const mBmonomial &lhs, const mBmonomial &rhs);

	mBmonomial monlcm(const mBmonomial &lhs, const mBmonomial &rhs);

	mBmonomial multiply (const mBmonomial &lhs, const mBmonomial &rhs);

	mBmonomial divide (const mBmonomial &dividend, const mBmonomial &divisor);

	pair<mBmonomial,_Scalar> multiply(const pair<mBmonomial,_Scalar> &lhs, const pair<mBmonomial,_Scalar> &rhs);

	pair<mBmonomial,_Scalar> divide(const pair<mBmonomial,_Scalar> &dividend, const pair<mBmonomial,_Scalar> &divisor);

	pair<mBmonomial,_Scalar> pairgcd(const pair<mBmonomial,_Scalar> &lhs, const pair<mBmonomial,_Scalar> &rhs);

	pair<mBmonomial,_Scalar> pairlcm(const pair<mBmonomial,_Scalar> &lhs, const pair<mBmonomial,_Scalar> &rhs);

	pair<unsigned int,_Scalar> taylorDerive(const mBmonomial &monomial, const mBmonomial &diffop);

//Computations:
	void quotientHilbertPoly(Poly<mpq_class> &hilbertPolynomial,const vector<Poly<_Scalar> > &relations);

	static bool groebnerDivide(Poly<_Scalar> f, vector<Poly<_Scalar> > &coeffs, Poly<_Scalar> &rest,const vector<Poly<_Scalar> > &generators);

	void reduce(vector<Poly<_Scalar> > &generators);

	void computeGroebnerBasis(vector<Poly<_Scalar> > &generators);

//Output:
	bool print(std::ostream &os,const mBmonomial monomial, bool longform = false);

};

#endif //SRC_MONOMBASIS_HPP
