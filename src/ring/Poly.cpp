//============================================================================
// Name        : Poly.cpp
// Author      : Jonathan Steinbuch
//============================================================================


#include "Poly.hpp"

#include "../helpers/Modulus.hpp"

	template <typename _Scalar>
	Poly<_Scalar>::Poly(PolyRing<_Scalar>* mB, string::const_iterator& it, const string& input, bool bracket)
	{
		enum inputState { fresh, add, mult, exp};
		poly.clear();

		inputState state = fresh;
		if(bracket)
		{
			(*this) = Poly(mB, it, input, false);
			state = add;
		}

		 while(it != input.end()){
			if(*it == '(')
			{
				it++;
				auto a = Poly(mB, it, input, true);
				if(state == mult || state == exp)
				{
					(*this) *= a;
					assert(this->getRing() == mB);
					state = mult;
				} else {

					if(state == fresh)
					{
						(*this) = a;
						assert(this->getRing() == mB);
						state = mult;
					} else {
						(*this) += a;
					}
				}
			} else if(*it == ')')
			{
				if(bracket)
				{
					it++;
				}
				return;
			} else if(isalpha(*it))
			{
				if(state == mult || state == exp)
				{
					(*this) *= Poly(mB, it, input,false);
					assert(this->getRing() == mB);
					state =  mult;
				} else if (state == add)
				{
					(*this) += Poly(mB, it, input,false);
				}
				else if (state == fresh)
				{
					unsigned int var = mB->getVar(input, it - input.begin());
					assert(var != mB->mB_ERROR);
					it += mB->getPrintVar(var).size();

					(*this) = Poly::variable(mB, var);
					state = exp;
				}
			}
			else if (isdigit(*it))
			{
				if(state == mult)
				{
					(*this) *= Poly(mB,it,input,false);
				}
				else if (state == add)
				{
					(*this) += Poly(mB, it, input,false);
				}
				else
				{
					string numstring;
					while (it != input.end() && isdigit(*it))
					{
						numstring += *it;
						it++;
					}
					if (state == fresh)
					{
						_Scalar num = (_Scalar) std::stoi(numstring);
					   *this = Poly(mB,num);
						assert(this->getRing() == mB);
						state = mult;
					}
					else if (state == exp)
					{
						unsigned int exp = std::stoi(numstring);
						(*this) ^= exp;
						assert(this->getRing() == mB);
						state = mult;
					}
				}
			} else if(*it == '+')
			{
				if(!bracket && (state == mult || state == exp))
				{
					return;
				}
				it++;
				auto a = Poly(mB,it,input,false);
				if(state == fresh)
				{
					(*this) = a;
					assert(this->getRing() == mB);
					state = mult;
				} else {
					(*this) += a;
					state = add;
				}

			}else if(*it == '*')
			{
				it++;
				assert(state != fresh);
				(*this) *= Poly(mB,it,input,false);
				assert(this->getRing() == mB);
				state = mult;
			} else if(*it == '^')
			{
				assert(state == exp);
				it++;
				state = exp;
			}else if(*it == '-')
			{
				if(!bracket && (state == mult || state == exp))
				{
					return;
				}
				it++;
				auto a = Poly(mB,it,input,false);
				if(state == fresh)
				{
					(*this) = ((_Scalar)-1)*a;
					assert(this->getRing() == mB);
					state = mult;
				} else {
					(*this) -= a;
					state = add;
				}
				state = add;
			}else if(*it == ' ')
			{
				it++;
			}  else if(*it != ')'){
				assert(false);
			}
		};

	}

	template <typename _Scalar>
	Poly<_Scalar>::Poly(PolyRing<_Scalar>* mB, const string& input) : Ring(mB)
	{
		auto it = input.cbegin();
		*this = Poly(mB,it,input, true);
	}

	template <typename _Scalar>
	Poly<_Scalar>& Poly<_Scalar>::trim()
	{
		Poly<_Scalar> nthis(Ring);
		for(auto &it : poly){
			if(it.second != 0)
			{
				Poly<_Scalar> lookup = Ring->getLookup(it.first,it.second);
				nthis += lookup;
				//assert(mB->getDegBasisIndex(it.first.id) != UINT_MAX);
			}
		}
		*this = nthis;
		return *this;
	}

	template <typename Scalar>
	void Poly<Scalar>::simplify()
		{
			Scalar factor = 0;
			for (auto &x : poly)
			{
				factor = gcd(x.second, factor);
				if (factor == 1 || factor == -1)
				{
					return;
				}
			}
			if(factor != 0)
			{
				for (auto &x : poly)
				{
					x.second /= factor;
				}
			}
		}

	template <typename _Scalar>
	bool Poly<_Scalar>::reduce(Poly<_Scalar>& reducer)
	{
		bool reduced = false;
		auto it = poly.rbegin();
		while( it != poly.rend())
		{
			mBmonomial curr = (*it).first;
			mBmonomial quotient = Ring->divide(curr,reducer.leadingMonomial());
			if(quotient != Ring->mB_ERROR){
				_Scalar g = gcd((*it).second,reducer.leadingCoefficient());
				_Scalar s = (*it).second / g;
				if(g != reducer.leadingCoefficient())
				{
					if(g == -reducer.leadingCoefficient()){
						s = -s;
					} else {
						_Scalar m = (reducer.leadingCoefficient()/g);
						(*this) *= m;
					}
				}
				Poly<_Scalar> a(Ring, quotient, s);
				*this -= a*reducer;
				reduced = true;
				auto itlow = poly.upper_bound(curr);
				it = make_reverse_iterator(itlow);
				if(it == poly.rend())
				{
					break;
				}
			}
			it++;
		}
		simplify();
		return reduced;
	}

	template <typename _Scalar>
	const unsigned int Poly<_Scalar>::degree() const
	{
		unsigned int deg = 0;
		for(auto &it : poly){
			unsigned int mdeg = Ring->getDegree(it.first);
			if(deg < mdeg)
			{
				deg = mdeg;
			}
		}
		return deg;
	}

	template <typename _Scalar>
	const pair<mBmonomial,_Scalar> Poly<_Scalar>::leadingTerm() const
	{
		if(isNull())
			return pair<mBmonomial,_Scalar>(0,0);
		return (*(--poly.end()));
	}

	template <typename _Scalar>
	const mBmonomial Poly<_Scalar>::leadingMonomial() const
	{
		if(isNull())
			return mBmonomial(0);
		return (*(--poly.end())).first;
	}

	template <typename _Scalar>
	const _Scalar Poly<_Scalar>::leadingCoefficient() const
	{
		if(isNull())
			return 0;
		return (*(--poly.end())).second;
	}

	template <typename _Scalar>
	const pair<mBmonomial,_Scalar> Poly<_Scalar>::trailingTerm() const
	{
		if(isNull())
			return pair<mBmonomial,_Scalar>(0,0);
		return (*(poly.begin()));
	}

	template <typename _Scalar>
	const mBmonomial Poly<_Scalar>::trailingMonomial() const
	{
		if(isNull())
			return mBmonomial(0);
		return (*(poly.begin())).first;
	}

	template <typename _Scalar>
	const _Scalar Poly<_Scalar>::trailingCoefficient() const
	{
		if(isNull())
			return 0;
		return (*(poly.begin())).second;
	}

	template <typename _Scalar>
	bool Poly<_Scalar>::popLeadingTerm()
	{
		if(isNull())
		{
			return false;
		}
		poly.erase(--poly.end());
		return true;
	}

	template <typename _Scalar>
	bool Poly<_Scalar>::isHomogeneous() const
	{
		bool set = false;
		unsigned int deg = 0;
		for(auto &it : poly){
			unsigned int mdeg = Ring->getDegree(it.first);
			if(set && deg != mdeg)
			{
				return false;
			}
			deg = mdeg;
			set = true;
		}
		return true;
	}


	template <typename _Scalar>
	const bool Poly<_Scalar>::isNull() const{
		return poly.size() == 0;
	}

	template <typename _Scalar>
	Poly< _Scalar>& Poly<_Scalar>::operator*=(const _Scalar &rhs)
	{
		auto it = poly.begin();
		while( it != poly.end())
		{
			(*it).second *= rhs;
			if((*it).second == 0)
				it = poly.erase(it);
			else
				it++;
		}
		return *this;
	}

	template <typename _Scalar>
	Poly< _Scalar>& Poly<_Scalar>::operator/=(const _Scalar &rhs)
	{
		assert(rhs != 0);
		auto it = poly.begin();
		while( it != poly.end())
		{
			_Scalar div = (*it).second / rhs;
			assert(div*rhs == (*it).second);
			(*it).second = div;
			if((*it).second == 0)
				it = poly.erase(it);
			else
				it++;
		}
		return *this;
	}

	template <typename _Scalar>
	Poly<_Scalar>& Poly<_Scalar>::monMultEq(const mBmonomial monomial)
	{
		Poly<_Scalar> nthis(Ring);
		auto it = poly.begin();
		for(;it != poly.end(); it++){
			mBmonomial product = Ring->multiply((*it).first,monomial);
			Poly<_Scalar> lookup = Ring->getLookup(product,(*it).second);
			nthis+=lookup;
		}
		*this = nthis;
		return *this;
	}

	template <typename _Scalar>
	Poly<_Scalar> Poly<_Scalar>::monMult( const mBmonomial monomial) const
	{
		Poly<_Scalar> lhs = *this;
		lhs.monMultEq(monomial);
		return lhs;
	}

	template <typename _Scalar>
	Poly<_Scalar>& Poly<_Scalar>::operator*=(const pair<mBmonomial,_Scalar>& rhs)
	{
		return monMultEq(rhs.first)*=rhs.second;

	}

	template <typename _Scalar>
	Poly<_Scalar> Poly<_Scalar>::operator*(const pair<mBmonomial,_Scalar>& rhs) const
	{
		return monMult(rhs.first)*rhs.second;
	}

	template <typename _Scalar>
	Poly<_Scalar> Poly<_Scalar>::operator*(const _Scalar &rhs) const
	{
		Poly<_Scalar> lhs = *this;
		lhs *= rhs;
		return lhs;
	}

	template <typename _Scalar>
	Poly<_Scalar> operator*(const _Scalar &lhs,Poly<_Scalar> rhs)
	{
		rhs *= lhs;
		return rhs;
	}
	template Poly<mpz_class> operator*(const mpz_class &lhs,Poly<mpz_class> rhs);
	template Poly<mpq_class> operator*(const mpq_class &lhs,Poly<mpq_class> rhs);


	template <typename _Scalar>
	Poly<_Scalar>& Poly<_Scalar>::operator*=(Poly<_Scalar> rhs)
	{
		Poly<_Scalar> nthis(Ring);
		auto it = rhs.poly.begin();
		for(;it != rhs.poly.end(); it++){
			nthis += (*this).monMult((*it).first)*(*it).second;
		}
		*this = nthis;
		return *this;
	}

	template <typename _Scalar>
	Poly<_Scalar> Poly<_Scalar>::operator*(const Poly<_Scalar>& rhs) const
	{
		Poly<_Scalar> lhs = *this;
		lhs *= rhs;
		return lhs;
	}

	template <typename _Scalar>
	Poly<_Scalar> Poly<_Scalar>::taylorDerive(mBmonomial monomial) const{
		Poly<_Scalar> ret(Ring);
		for(auto& monom : poly){
			pair<mBmonomial,_Scalar> nmonom = Ring->taylorDerive(monom.first,monomial);
			nmonom.second*=monom.second;
			ret += nmonom;
		}
		return ret;
	}

	template <typename _Scalar>
	Poly<_Scalar>& Poly<_Scalar>::operator+=(const pair<mBmonomial,_Scalar>& rhs)
	{
		if(rhs.second != 0)
		{
			auto ret = poly.insert(pair<mBmonomial,_Scalar>(rhs.first,rhs.second));
			if(!ret.second)
			{
				(*(ret.first)).second += rhs.second;
				if((*(ret.first)).second == 0)
					poly.erase(ret.first);
			}
		}
		return *this;
	}

	template <typename _Scalar>
	Poly<_Scalar>& Poly<_Scalar>::operator+=(const Poly<_Scalar>& rhs)
	{

		auto it = rhs.poly.cbegin();
		for(;it != rhs.poly.cend(); it++)
		{
			*this += (*it);
		}
		return *this;
	}

	template <typename _Scalar>
	Poly<_Scalar> Poly<_Scalar>::operator+(const Poly<_Scalar>&  rhs) const
	{
		Poly<_Scalar> lhs = *this;
		lhs += rhs;
		return lhs;
	}

	template <typename _Scalar>
	Poly<_Scalar>& Poly<_Scalar>::operator-=(const Poly<_Scalar>& rhs)
	{
		*this += ((_Scalar)(-1))*rhs;
		return *this;
	}

	template <typename _Scalar>
	Poly<_Scalar> Poly<_Scalar>::operator-(const Poly<_Scalar>&  rhs) const
	{
		Poly<_Scalar> lhs = *this;
		lhs-=rhs;
		return lhs;
	}

	template <typename _Scalar>
	Poly<_Scalar> Poly<_Scalar>::monAdd(const mBmonomial rhs) const
	{
		Poly<_Scalar> lhs = *this;
		return lhs+=pair<mBmonomial,_Scalar>(rhs,1);
	}

	template <typename _Scalar>
	Poly<_Scalar>& Poly<_Scalar>::operator^=(const unsigned int exponent)
	{
		if(exponent > 1)
		{
			Poly<_Scalar> out = *this;
			for(unsigned int i=1; i < exponent; i++)
			{
				(*this) *= out;
			}
			return *this;
		} else if(exponent == 1){
			return *this;
		} else {
			throw domain_error("^0");
		}
	}

	template <typename _Scalar>
	bool Poly<_Scalar>::operator==(const Poly<_Scalar>& rhs){
		Poly<_Scalar> difference = (*this)-rhs;
		difference.trim();
		return difference.isNull();
	}

	template <typename _Scalar>
	Poly<_Scalar> Poly<_Scalar>::operator^(const unsigned int exponent) const
	{
		Poly<_Scalar> lhs = *this;
		lhs ^= exponent;
		return lhs;
	}

	template <typename _Scalar>
	Poly<_Scalar> Poly<_Scalar>::variable(PolyRing<_Scalar>* mB, const unsigned int i)
	{
		vector<unsigned int> v(mB->getNumVars(),0);
		v[i] = 1;
		return Poly<_Scalar>(mB,mB->getMonomial(v),1);
	}

	template <typename _Scalar>
	Poly<_Scalar> Poly<_Scalar>::generic(PolyRing<_Scalar>* mB, const unsigned int degree)
	{
		Poly<_Scalar> ret(mB);
		for(unsigned int i = mB->getIndexOfBasisDegStart(degree); i < mB->getIndexOfBasisDegStart(degree+1);i++)
		{
			const mBmonomial mx = mB->getBasisElem(i);
			ret += pair<mBmonomial,_Scalar>(mx,rand()%5-2);//*rand()
		}
		return ret;
	}

	template <typename _Scalar>
	template <typename tripScalar>
	vector<rcvTriplet<tripScalar> >& Poly<_Scalar>::appendMatrix(vector<rcvTriplet<tripScalar> >& ret, const int outdegree, const int indegree, unsigned int rows, unsigned int cols, const unsigned int rowoffset, const unsigned int coloffset, bool homogeneous) const
	{
		if (homogeneous)
		{

			if (indegree == outdegree - (int) degree())
			{
				for (unsigned int i = 0; i < cols; i++)
				{

					Poly<_Scalar> product = (*this).monMult(Ring->getBasisElem(Ring->getIndexOfBasisDegStart(indegree) + i));
					for (auto& mon : product.poly)
					{
						unsigned int j = Ring->getDegBasisIndex(outdegree, mon.first) - Ring->getIndexOfBasisDegStart(outdegree);
						ret.push_back(rcvTriplet<tripScalar>(j + rowoffset, i + coloffset, (tripScalar) mon.second));
					}
				}
			}
		}
		else
		{

			for (unsigned int i = 0; i < cols; i++)
			{

				Poly<_Scalar> product = (*this).monMult(Ring->getBasisElem(i));
				for (auto& mon : product.poly)
				{
					unsigned int j = Ring->getDegBasisIndex(mon.first);
					assert(j < rows);
					ret.push_back(rcvTriplet<tripScalar>(j + rowoffset, i + coloffset, (tripScalar) mon.second));
				}
			}
		}
		return ret;
	}

	template <typename _Scalar>
	void Poly<_Scalar>::homogenize(unsigned int variable)
	{
		variable = Ring->getNumVars() - variable-1;
		polyContainer<_Scalar> newpoly;
		unsigned int d = degree();
		vector<unsigned int> vec;
		vec.resize(Ring->getNumVars(),0);
		for(auto & mon : poly)
		{
			unsigned int mond = Ring->getDegree(mon.first);
			if(mond < d)
			{
				vec[variable]= d-mond;
				mBmonomial x = Ring->multiply(mon.first,Ring->getMonomial(vec));
				newpoly.insert(pair<mBmonomial,_Scalar>(x,mon.second));
			} else {
				newpoly.insert(mon);
			}
		}
		poly = newpoly;
	}

	template <typename _Scalar>
	vector<rcvTriplet<_Scalar> > Poly<_Scalar>::getMatrix(const int outdegree,const int indegree, unsigned int rows, unsigned int cols, const unsigned int rowoffset, const unsigned int coloffset){
		vector<rcvTriplet<_Scalar> > ret;
		appendMatrix(ret,outdegree,indegree,rows,cols,rowoffset,coloffset);
		return ret;
	}

	template <typename _Scalar>
	std::ostream &operator<<(std::ostream &os, const Poly<_Scalar>  &m) {
		m.print(os,m.getRing()->opt->longFormPolynomials);
	   return os;
	}

template std::ostream &operator<< <mpz_class>(std::ostream &os, const Poly<mpz_class>  &m);
template std::ostream &operator<< <mpq_class>(std::ostream &os, const Poly<mpq_class>  &m);
template std::ostream &operator<< <numbermodulo>(std::ostream &os, const Poly<numbermodulo>  &m);

template<typename _Scalar>
void Poly<_Scalar>::print(std::ostream &os, bool longform) const
{
	if(poly.size() == 0){
		os << 0;
	} else {
		auto it = poly.cbegin();
		for(;it != poly.cend(); it++)
		{
			_Scalar coeff = (*it).second;
			if(coeff >= 0 && it != poly.cbegin())
			{
				os << "+";
			} else if(coeff < 0){
				os << "-";
				coeff *= -1;
			}
			if(coeff != 1)
			{
				os << coeff;
				if(!(((*it).first).isOne()) && longform && !Ring->opt->outputLatex)
				{
					os << "*";
				}
			}

			bool visible = Ring->print(os,(*it).first, longform);

			if(!visible && coeff == 1)
			{
				os << 1;
			}
		}

	}
}

template class Poly <mpz_class>;
template class Poly <mpq_class>;
template class Poly<numbermodulo> ;
