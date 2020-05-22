//============================================================================
// Name        : PolyRing.cpp
// Author      : Jonathan Steinbuch
//============================================================================


#include "PolyRing.hpp"

#include "../helpers/Modulus.hpp"
#include "Poly.hpp"

monlookups::monlookups(const unsigned int vars, const vector<string> (&printvar), const unsigned int order, int realvars, const vector<int> &weights, bool invert) : lu(vars), vars(vars), order(order), realvars(realvars)
{
	list.clear();
	if (realvars < 0)
	{
		this->realvars = vars;
	}
	if (invert)
	{
		this->printvar = printvar;
		this->weights = weights;
		if (weights.size() < vars)
		{
			this->weights.resize(vars, 1);
		}
	}
	else
	{
		for (unsigned int i = 0; i < vars; i++)
		{
			this->printvar.push_back(printvar[vars - i - 1]);
			if (weights.size() <= i)
			{
				this->weights.push_back(1);
			}
			else
			{
				this->weights.push_back(weights[vars - i - 1]);
			}
		}
	}

	maxdegree = 0;
}

template<typename _Scalar>
PolyRing<_Scalar>::PolyRing(const unsigned int vars, const vector<string> (&printvar), const unsigned int order, int realvars, const vector<int> &weights, bool invert)
{
	BaseRing = this;
	mons = new monlookups(vars,printvar,order, realvars, weights,invert);

	reSize(1);
}

template<typename _Scalar>
PolyRing<_Scalar>::PolyRing(PolyRing* BaseRing, const vector<Poly<_Scalar> > &relations) : BaseRing(BaseRing){
	mons = BaseRing->mons;

	setRelations(relations);
}


template<typename _Scalar>
PolyRing<_Scalar>::PolyRing(PolyRing* BaseRing, const vector<string> &relations) :
		BaseRing(BaseRing)
{
	mons = BaseRing->mons;

	vector<Poly<_Scalar> > polrels;
	for (auto &strrelation : relations)
	{
		Poly<_Scalar> polrel = Poly<_Scalar>(BaseRing, strrelation);
		polrels.push_back(polrel);
	}

	setRelations(polrels);
}




template<typename _Scalar>
PolyRing<_Scalar>::~PolyRing(){
	if(BaseRing == this)
	{
		delete mons;
	}
}


template<typename _Scalar>
void PolyRing<_Scalar>::reCalcRelations(const unsigned int minDegree)
{
	if(minDegree >= mons->maxdegree) return;
	vector<mBmonomial> nbasis;
	unsigned int nbasisSizeOffset = 0;
	nbasisSizeOffset = basisdegstart[minDegree];

	vector<unsigned int> nbasisdegstart;

	unsigned int degree = minDegree;

	for (unsigned int i = nbasisSizeOffset; i < basis.size(); i++)
	{

		while (degree < mons->maxdegree && i >= basisdegstart[degree + 1])
		{
			degree++;
			nbasisdegstart.push_back(nbasis.size()+nbasisSizeOffset);
		}
		mBmonomial elem = basis[i];

		bool inIdeal  = false;
		for(auto & gbElem : gB)
		{
			mBmonomial quotient = divide(elem, gbElem.leadingMonomial());
			if(quotient != mBmonomial::mon_ERROR)
			{
				inIdeal = true;
				Poly<_Scalar> rest(this);
				lookup[elem.getId()] = 0;
				rest = gbElem;
				_Scalar s = rest.leadingCoefficient();
				rest.popLeadingTerm();
				if(s != 1)
				{
					rest /= s;
				}
				while(!rest.isNull())
				{
					mBmonomial rMonomial = multiply(rest.leadingMonomial(),quotient);
					lookup[elem.getId()] += (lookup[rMonomial.getId()])*(-(rest.leadingCoefficient()));
					rest.popLeadingTerm();
				}
				break;
			}
		}
		if (!inIdeal){
			nbasis.push_back(elem);
		}

		vector<Poly<_Scalar> > coeffs;
		coeffs.assign(gB.size(), Poly<_Scalar>(this));
		Poly<_Scalar> rest(this);
		bool inIdeal2 = groebnerDivide(Poly<_Scalar>(this->BaseRing, elem, 1), coeffs, rest,gB);
		assert((inIdeal2 || rest.leadingMonomial() < elem ) == inIdeal);

	}

	basis.resize(nbasisSizeOffset,0);
	basis.insert(basis.end(),nbasis.begin(),nbasis.end());
	basisdegstart.resize(minDegree+1,0);
	basisdegstart.insert(basisdegstart.end(),nbasisdegstart.begin(),nbasisdegstart.end());
}

template<typename _Scalar>
void PolyRing<_Scalar>::reSizeBasis()
{
	if(mons->maxdegree < basisdegstart.size())
		return;
	for (unsigned int degree = basisdegstart.size(); degree <= mons->maxdegree; degree++)
	{
		basisdegstart.push_back(basis.size());

		for (unsigned int i = mons->monN[degree]; i < mons->monN[degree+1]; i++)
		{
			basis.push_back(i);
			lookup.push_back(Poly<_Scalar>(this, i, 1));
		}
	}
}

template<typename _Scalar>
void PolyRing<_Scalar>::reSize(const unsigned int neededDegree)
{
	unsigned int startdeg = mons->maxdegree + 1;
	if (mons->maxdegree < neededDegree)
	{
		unsigned int newMaxDegree = mons->maxdegree;
		if (mons->maxdegree == 0)
		{
			startdeg = 0;
			newMaxDegree = 1;
			mons->monN.push_back(0);
		}
		while (newMaxDegree <= neededDegree)
		{
			newMaxDegree *= 2;
		}

		for (unsigned int degree = startdeg; degree <= newMaxDegree; degree++)
		{
			vector<unsigned int> mOff;

			kSumTon(mons->list, degree, mons->vars, mons->order, mons->weights, mOff);

			mons->monomOffset.push_back(mOff);
			mons->monN.push_back(mons->list.size());
		}
		mons->maxdegree = newMaxDegree;
	}

	startdeg = basisdegstart.size();
	reSizeBasis();
	if(gB.size() > 0)
	{
		reCalcRelations(startdeg);
	}
}

template<typename _Scalar>
inline void PolyRing<_Scalar>::reSizeToIndex(const unsigned int neededIndex)
{
	while(neededIndex >= mons->list.size())
	{
		reSize(mons->maxdegree*2);
	}
}

template<typename _Scalar>
vector<unsigned int> PolyRing<_Scalar>::getHilbert()
{
	vector<unsigned int> hilbert;
	cout << "Hilbertfunktion:";
	for (unsigned int i = 0; i < basisdegstart.size() - 1; i++)
	{
		hilbert.push_back(basisdegstart[i + 1] - basisdegstart[i]);
		cout << " " << i << ":" << hilbert[i];
	}
	cout << endl << endl;
	return hilbert;
}

template<typename _Scalar>
unsigned int PolyRing<_Scalar>::getDegree(const vector<unsigned int> &x) const
{
	unsigned int deg = 0;
	for (unsigned int j = 0; j < x.size(); j++)
	{
		deg += x[j] * mons->weights[j];
	}
	return deg;
}

template<typename _Scalar>
unsigned int PolyRing<_Scalar>::getDegree(const unsigned int index) const
{
	assert(index < mons->list.size());

	return getDegree(mons->list[index]);
}

template<typename _Scalar>
unsigned int PolyRing<_Scalar>::getDegree(const mBmonomial monomial) const
{
	return getDegree(monomial.getId());
}

template<typename _Scalar>
unsigned int PolyRing<_Scalar>::getIndex(const vector<unsigned int> &x)
{
	return getIndex(getDegree(x), x);
}

template<typename _Scalar>
mBmonomial PolyRing<_Scalar>::getMonomial(const vector<unsigned int> &x)
{
	return mBmonomial(getIndex(x));
}

template<typename _Scalar>
unsigned int PolyRing<_Scalar>::getVar(string var, unsigned int pos)
{
	for (unsigned int i = 0; i < mons->printvar.size(); i++)
	{
		if (var.find(mons->printvar[i], pos) == pos)
		{
			return i;
		}
	}
	return mB_ERROR;
}

template<typename _Scalar>
unsigned int PolyRing<_Scalar>::getIndex(const unsigned int degree, const vector<unsigned int> &x)
{
	reSize(degree);

	unsigned int pos = 0;
	switch (mons->order)
	{
	case degrevlex:
	{
		unsigned int nj = degree;
		pos = mons->monN[nj];
		for (unsigned int j = 0; j < mons->vars - 1; j++)
		{
			if (j < x.size())
			{
				nj -= x[j] * mons->weights[j];
			}
			unsigned int krest = j;
			pos += mons->monomOffset[nj][krest];
		}
		break;
	}
	case deglex:
	{
		unsigned int nj = degree;
		pos = mons->monN[nj + 1] - 1;
		for (unsigned int j = mons->vars - 1; j > 0; --j)
		{
			if (j < x.size())
			{
				nj -= x[j] * mons->weights[j];
			}
			unsigned int krest = mons->vars - j - 1;
			pos -= mons->monomOffset[nj][krest - 1];
		}
		break;
	}
	default:
	{
		return false;
		break;
	}
	}
	return pos;
}

template<typename _Scalar>
mBmonomial PolyRing<_Scalar>::getMonomial(const unsigned int degree, const vector<unsigned int> &x)
{
	reSize(degree);
	return mBmonomial(getIndex(degree, x));
}

/*	unsigned int getIndex(const monom<mons->vars,_Scalar>& monomial)
 {
 return getIndex(monomial.x);
 }*/

template<typename _Scalar>
unsigned int PolyRing<_Scalar>::getBasisIndexfromDegBasisIndex(const unsigned int degree, const unsigned int degBasisIndex) const
{
	assert(degree < mons->maxdegree);

	unsigned int ret = basisdegstart[degree] + degBasisIndex;
	assert(ret < basisdegstart[degree + 1]);
	return ret;
}

template<typename _Scalar>
mBmonomial PolyRing<_Scalar>::getMonomialfromBasis(const unsigned int degree, const unsigned int degBasisIndex) const
{
	return basis[getBasisIndexfromDegBasisIndex(degree, degBasisIndex)];
}

template<typename _Scalar>
unsigned int PolyRing<_Scalar>::getDegBasisIndex(const mBmonomial monomial) const
{
	const unsigned int Index = monomial.getId();
	unsigned int degree = getDegree(Index);
	return getDegBasisIndex(degree, Index);
}

template<typename _Scalar>
unsigned int PolyRing<_Scalar>::getDegBasisIndex(const unsigned int degree, const mBmonomial monomial) const
{
	assert(degree < basisdegstart.size());
	auto b = lower_bound(basis.begin() + basisdegstart[degree], basis.begin() + basisdegstart[degree + 1], monomial.getId());
	//assert(*b == Index);
	return distance(basis.begin(), b);
}

template<typename _Scalar>
inline const Poly<_Scalar>& PolyRing<_Scalar>::getLookup(const mBmonomial monomial)
{
	reSizeToIndex(monomial.getId());
	return lookup[monomial.getId()];
}

template<typename _Scalar>
inline const Poly<_Scalar> PolyRing<_Scalar>::getLookup(const mBmonomial monomial, const _Scalar coefficient)
{
	reSizeToIndex(monomial.getId());
	return lookup[monomial.getId()] * coefficient;
}

template<typename _Scalar>
inline const vector<unsigned int> &PolyRing<_Scalar>::getExponents(const mBmonomial monomial)
{
	reSizeToIndex(monomial.getId());
	return mons->list[monomial.getId()];
}

/*	const mBpolynomial<_Scalar> monomBasis<_Scalar>::getLookup(const monom<mons->vars,_Scalar>& monomial,const _Scalar coefficient) const
 {
 unsigned int index = getIndex(monomial);
 return getLookup(index,coefficient);
 }*/

template<typename _Scalar>
mBmonomial PolyRing<_Scalar>::mongcd(const mBmonomial &lhs, const mBmonomial &rhs)
{
	const vector<unsigned int> &lhsX = getExponents(lhs);
	const vector<unsigned int> &rhsX = getExponents(rhs);

	vector<unsigned int> x(mons->vars);
	unsigned int degree = 0;
	for (unsigned int i = 0; i < mons->vars; i++)
	{
		x[i] = min(lhsX[i], rhsX[i]);
		degree += x[i] * mons->weights[i];
	}
	return getMonomial(degree, x);
}

template<typename _Scalar>
pair<mBmonomial,_Scalar> PolyRing<_Scalar>::pairgcd(const pair<mBmonomial,_Scalar> &lhs, const pair<mBmonomial,_Scalar> &rhs)
{
	mBmonomial mon = mongcd(lhs.first,rhs.first);
	_Scalar scal = gcd(lhs.second,rhs.second);
	return pair<mBmonomial,_Scalar>(mon,scal);
}

template<typename _Scalar>
pair<mBmonomial,_Scalar> PolyRing<_Scalar>::pairlcm(const pair<mBmonomial,_Scalar> &lhs, const pair<mBmonomial,_Scalar> &rhs)
{
	mBmonomial mon = monlcm(lhs.first,rhs.first);
	_Scalar scal = lcm(lhs.second,rhs.second);
	return pair<mBmonomial,_Scalar>(mon,scal);
}

template<typename _Scalar>
mBmonomial PolyRing<_Scalar>::monlcm(const mBmonomial &lhs, const mBmonomial &rhs)
{
	const vector<unsigned int> &lhsX = getExponents(lhs);
	const vector<unsigned int> &rhsX = getExponents(rhs);

	vector<unsigned int> x(mons->vars);
	unsigned int degree = 0;
	for (unsigned int i = 0; i < mons->vars; i++)
	{
		x[i] = max(lhsX[i], rhsX[i]);
		degree += x[i] * mons->weights[i];
	}
	return getMonomial(degree, x);
}

template<typename _Scalar>
pair<mBmonomial,_Scalar> PolyRing<_Scalar>::multiply(const pair<mBmonomial,_Scalar> &lhs, const pair<mBmonomial,_Scalar> &rhs)
{
	mBmonomial mon = multiply(lhs.first,rhs.first);
	_Scalar scal = lhs.second*rhs.second;
	return pair<mBmonomial,_Scalar>(mon,scal);
}

template<typename _Scalar>
mBmonomial PolyRing<_Scalar>::multiply(const mBmonomial &lhs, const mBmonomial &rhs)
{
	const vector<unsigned int> &lhsX = getExponents(lhs);
	const vector<unsigned int> &rhsX = getExponents(rhs);

	vector<unsigned int> x(mons->vars);
	unsigned int degree = 0;
	for (unsigned int i = 0; i < mons->vars; i++)
	{
		x[i] = lhsX[i] + rhsX[i];
		degree += x[i] * mons->weights[i];
	}
	return getMonomial(degree, x);
}

template<typename _Scalar>
pair<mBmonomial,_Scalar> PolyRing<_Scalar>::divide(const pair<mBmonomial,_Scalar> &dividend, const pair<mBmonomial,_Scalar> &divisor)
{
	mBmonomial mon = divide(dividend.first,divisor.first);
	_Scalar scal = dividend.second/divisor.second;
	return pair<mBmonomial,_Scalar>(mon,scal);
}

template<typename _Scalar>
mBmonomial PolyRing<_Scalar>::divide(const mBmonomial &dividend, const mBmonomial &divisor)
{
	const vector<unsigned int> &dividendX = getExponents(dividend);
	const vector<unsigned int> &divisorX = getExponents(divisor);

	vector<unsigned int> x(mons->vars);
	for (unsigned int i = 0; i < mons->vars; i++)
	{

		if (dividendX[i] >= divisorX[i])
			x[i] = dividendX[i] - divisorX[i];
		else
		{
			return mBmonomial(mBmonomial::mon_ERROR);
		}
	}
	return getMonomial(x);
}

template<typename _Scalar>
pair<unsigned int, _Scalar> PolyRing<_Scalar>::taylorDerive(const mBmonomial &monomial, const mBmonomial &diffop)
{ //computes the derivative divided by the factorial of the exponents in the differential operator
	if (diffop.isOne())
		return pair<unsigned int, _Scalar>(0, 0);
	vector<unsigned int> x(mons->vars);
	_Scalar prefactor = 1;
	for (unsigned int i = 0; i < mons->vars; i++)
	{
		unsigned int mi = getExponents(monomial)[i];
		unsigned int di = getExponents(diffop)[i];
		if (mi >= di)
		{
			x[i] = mi - di;
			prefactor *= mons->lu.choose(mi, di);
		}
		else
		{
			return pair<unsigned int, _Scalar>(0, 0);
		}
	}
	return pair<unsigned int, _Scalar>(getIndex(x), prefactor);
}

template<typename _Scalar>
unsigned int PolyRing<_Scalar>::getBasisSize(const int degree)
{
	if (degree < 0)
	{
		return 0;
	}
	else
	{
		reSize(degree + 1);
		return basisdegstart[degree + 1] - basisdegstart[degree];
	}
}

template<typename _Scalar>
bool PolyRing<_Scalar>::isSmooth()
{
	vector<vector<Poly<_Scalar> > > entries;
	entries.push_back(gB);
	DensePolyMatrix<_Scalar> M(BaseRing, 1, gB.size(), entries);
	DensePolyMatrix<_Scalar> Jacobi = M.jacobiMatrix();

	vector<Poly<_Scalar> > jacideal = Jacobi.computeMinors(min(min((unsigned int)2,Jacobi.getCols()),Jacobi.getRows()));
	jacideal.insert(jacideal.end(),gB.begin(),gB.end());
	BaseRing->reduce(jacideal);
	BaseRing->computeGroebnerBasis(jacideal);

	PolyRing<mpq_class> kX(1,
	{ "n" }, degrevlex);
	Poly<mpq_class> hP(&kX);
	BaseRing->quotientHilbertPoly(hP,jacideal);

	if(hP.isNull())
	{
		return true;
	} else {
		return false;
	}
}

template<typename _Scalar>
void PolyRing<_Scalar>::monpoly(const int deg, const unsigned int vars, Poly<mpq_class> &retVal, const Poly<mpq_class> &n)
{
	mpq_class divisor = 1;
	Poly<mpq_class> numerator = Poly<mpq_class>(n.getRing(), 1);
	for (int i = 1; i < (int) vars; i++)
	{
		Poly<mpq_class> d = Poly<mpq_class>(n.getRing(), i - deg);
		numerator *= (n + d);
		divisor *= i;
	}
	//numerator.trim();
	numerator /= divisor;
	retVal = numerator;
}

template<typename _Scalar>
void PolyRing<_Scalar>::quotientHilbertPoly(Poly<mpq_class> &hilbertPolynomial,const vector<Poly<_Scalar> > &relations)
{
	Poly<mpq_class> n = Poly<mpq_class>::variable(hilbertPolynomial.getRing(), 0);

	vector<map<mBmonomial, int> > leads;
	map<mBmonomial, int> firstleads;
	for (auto & elem : relations)
	{
		firstleads.insert(pair<mBmonomial, int>(elem.leadingMonomial(), -1));
	}
	leads.push_back(firstleads);

	vector<int> degrees;
	degrees.push_back(1);

	monpoly(0, mons->vars, hilbertPolynomial, n);

	int factor = -1;
	for (unsigned int iter = 0; leads[iter].size() != 0; iter++)
	{
		for (auto & elem : leads[iter])
		{
			for (unsigned int i = 0; i < leads.size() - 1; i++)
			{
				for (map<mBmonomial, int>::iterator j = leads[i].begin(); j != leads[i].end(); j++)
				{
					mBmonomial quotient = divide(elem.first, (*j).first);
					if (!quotient.isError())
					{
						elem.second -= (*j).second;
					}
				}
			}
			unsigned int d = getDegree(elem.first);
			if (d >= degrees.size())
			{
				degrees.resize(d + 1, 0);
			}
			degrees[d] += elem.second;
		}
		map<mBmonomial, int> nextleads;
		for (map<mBmonomial, int>::iterator i = leads[iter].begin(); i != leads[iter].end(); i++)
		{
			for (map<mBmonomial, int>::iterator j = next(i); j != leads[iter].end(); j++)
			{
				mBmonomial l = monlcm((*i).first, (*j).first);
				nextleads.insert(pair<mBmonomial, int>(l, -1));
			}
		}
		leads.push_back(nextleads);
		factor = -factor;
	}

	hilbertPolynomial = Poly<mpq_class>(n.getRing());

	for (unsigned int i = 0; i < degrees.size(); i++)
	{
		if (degrees[i] != 0)
		{
			Poly<mpq_class> retVal(n.getRing());
			monpoly(i, mons->vars, retVal, n);
			mpq_class factor = degrees[i];
			hilbertPolynomial += factor * retVal;
		}
	}
}

template<typename _Scalar>
bool PolyRing<_Scalar>::groebnerDivide(Poly<_Scalar> f, vector<Poly<_Scalar> > &coeffs, Poly<_Scalar> &rest, const vector<Poly<_Scalar> > &generators )
{
	Poly<_Scalar> fcopy = f;
	assert(rest.isNull());
	for (unsigned int i = 0; i < generators.size(); i++)
	{
		assert(coeffs[i].isNull());
	}
	while (!f.isNull())
	{
		for (unsigned int i = 0; i < generators.size(); i++)
		{
			mBmonomial quotient = f.getRing()->divide(f.leadingMonomial(), generators[i].leadingMonomial());
			if (!quotient.isError()) //not divisible by remindex
			{
				_Scalar g = gcd(f.leadingCoefficient(),generators[i].leadingCoefficient());
				_Scalar s = f.leadingCoefficient() / g;
				if(g != generators[i].leadingCoefficient())
				{
					if(g == -generators[i].leadingCoefficient()){
						s = -s;
					} else {
						_Scalar m = (generators[i].leadingCoefficient()/g);
						f *= m;
						for (unsigned int i = 0; i < generators.size(); i++)
						{
							coeffs[i]*=m;
						}
						fcopy*=m;
					}
				}
				Poly<_Scalar> a(f.getRing(), quotient, s);
				coeffs[i] += a;

				Poly<_Scalar> minus = a * generators[i];

				f -= minus;

				if(f.isNull())
				{
					break;
				}

				i = 0;

			}
		}

		if (!f.isNull())
		{
			rest += f.leadingTerm();
			f.popLeadingTerm();
		}
	}
	rest.trim();

	return rest.isNull();
}

template<typename _Scalar>
Poly<_Scalar> PolyRing<_Scalar>::groebnerS(const Poly<_Scalar> &poly1, const Poly<_Scalar> &poly2){
	const pair<mBmonomial,_Scalar> l1 = poly1.leadingTerm();
	const pair<mBmonomial,_Scalar> l2 = poly2.leadingTerm();

	const pair<mBmonomial,_Scalar> lcm12 = pairlcm(l1,l2);

	const pair<mBmonomial,_Scalar> f1 = divide(lcm12,l1);
	const pair<mBmonomial,_Scalar> f2 = divide(lcm12,l2);

	Poly<_Scalar> retVal = poly1*f1-poly2*f2;
	return retVal;
}

template<typename _Scalar>
void PolyRing<_Scalar>::reduce(vector<Poly<_Scalar> > &generators)
{
	for (unsigned int i = 0; i < generators.size(); i++)
	{
		generators[i].simplify();
	}

	bool hasChanged = true;
	while (hasChanged == true)
	{
		hasChanged = false;
		for (unsigned int i = 0; i < generators.size(); i++)
		{
			if (!generators[i].isNull())
			{
				for (unsigned int j = 0; j < generators.size(); j++)
				{
					if (j != i && !generators[j].isNull())
					{
						bool reduced = generators[i].reduce(generators[j]);
						hasChanged = hasChanged || reduced;
					}
				}
			}
		}
	}

	vector<Poly<_Scalar> > newgens;
	for (unsigned int i = 0; i < generators.size(); i++)
	{
		if (!generators[i].isNull())
		{
			newgens.push_back(generators[i]);
		}
	}
	generators = newgens;
}

template<typename _Scalar>
void PolyRing<_Scalar>::computeGroebnerBasis(vector<Poly<_Scalar> > &generators)
{
	bool hasAdded = true;

	while (hasAdded == true)
	{
		hasAdded = false;
		for (unsigned int i = 0; i < generators.size(); i++)
		{
			for (unsigned int j = i+1; j < generators.size(); j++)
			{
				Poly<_Scalar> S = groebnerS(generators[i],generators[j]);
				Poly<_Scalar> rest(this);
				vector<Poly<_Scalar> > coeffs(generators.size(),this);
				if(!groebnerDivide(S,coeffs,rest,generators))
				{
					generators.push_back(rest);
					hasAdded = true;
				}
			}
		}
	}

	reduce(generators);
}

template<typename _Scalar>
void PolyRing<_Scalar>::setRelations(const vector<Poly<_Scalar> > &generators)
{
	gB = generators;
	BaseRing->computeGroebnerBasis(gB);
	reSize(mons->maxdegree);

	reCalcRelations(0);
}

template<typename _Scalar>
bool PolyRing<_Scalar>::print(std::ostream &os, const mBmonomial monomial, bool longform)
{
	const vector<unsigned int> &exponents = getExponents(monomial);
	bool visible = false;
	for (unsigned int i = 0; i < mons->realvars; i++)
	{
		if (exponents[i] != 0)
		{
			if (longform == true && visible == true && !(opt->outputLatex))
				os << "*";
			visible = true;
			os << mons->printvar[i];
			if (exponents[i] != 1)
			{
				if (longform == true)
					os << "^";
				os << exponents[i];
			}
		}
	}
	return visible;
}

template class PolyRing<mpz_class> ;
template class PolyRing<mpq_class> ;
template class PolyRing<numbermodulo> ;
