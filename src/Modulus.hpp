//============================================================================
// Name        : Modulus.hpp
// Author      : Jonathan Steinbuch
//============================================================================

#ifndef SRC_MODULUS_HPP_
#define SRC_MODULUS_HPP_

#include <numeric>

using namespace std;

typedef unsigned long datatype;

//Table automating inverses
struct invTable
{
	void set(datatype modulus)
	{
		invs.resize(modulus,0);
		//preset to 0
		for (datatype a = 0; a < modulus; ++a)
		{
			invs[a] = 0;
		}
		//actually compute inverses
		for (datatype a = 1; a < modulus; ++a)
		{
			for (datatype b = a; b < modulus; ++b)
			{
				if(a*b%modulus == 1)
				{
					invs[a] = b;
					invs[b] = a;
				}
			}
		}
	}

	vector<datatype> invs;
};

//Class implementing, F_p for p a prime.
class numbermodulo {
	datatype data; //the value

	//Attention! These values being static means we can only compute with one modulus at a time
	static datatype modulus; //the modulus
	static invTable lookup; //the inversion lookup table

public:
//Constructors:
	static void setModulus(datatype modulus)
	{
		numbermodulo::modulus = modulus;
		lookup.set(modulus);
	}

	numbermodulo(long input){
		if(input < 0)
		{
			data = -input;
			mod();
			data = modulus - data;
		} else {
			data = input;
			mod();
		}
	}

	numbermodulo(){
		data = 0;
	}

//Arithmetic:
	inline void mod()
	{
		data = data % modulus;
	}

	numbermodulo& operator-=(const numbermodulo& rhs)
	{
		data = (data + modulus) - rhs.data;
		mod();
		return *this;
	}

	numbermodulo operator-() const
	{
		return numbermodulo(modulus-data);
	}

	friend numbermodulo operator-(numbermodulo lhs,const numbermodulo&  rhs)
	{
		lhs-=rhs;
		return lhs;
	}

	numbermodulo& operator+=(const numbermodulo& rhs)
	{
		data += rhs.data;
		mod();
		return *this;
	}

	bool operator==(const numbermodulo& rhs) const
	{
		return data == rhs.data;
	}

	bool operator!=(const numbermodulo& rhs) const
	{
		return !(*this == rhs);
	}

	bool operator<(const numbermodulo& rhs) const
	{
		return data < rhs.data;
	}

	bool operator>(const numbermodulo& rhs) const
	{
		return data < rhs.data;
	}

	bool operator>=(const numbermodulo& rhs) const
	{
		return data >= rhs.data;
	}

	friend inline std::ostream &operator<<(std::ostream &os, const numbermodulo  &m)
	{
		os << m.data;
		return os;
	}

	friend numbermodulo operator+(numbermodulo lhs,const numbermodulo&  rhs)
	{
		lhs+=rhs;
		return lhs;
	}

	numbermodulo& operator*=(const numbermodulo& rhs)
	{
		data *= rhs.data;
		mod();
		return *this;
	}

	friend numbermodulo operator*(numbermodulo lhs,const numbermodulo&  rhs)
	{
		lhs*=rhs;
		return lhs;
	}

	numbermodulo& operator/=(const numbermodulo& rhs)
	{
		assert(lookup.invs[rhs.data] != 0);
		data*=lookup.invs[rhs.data];
		return *this;
	}

	friend numbermodulo operator/(numbermodulo lhs,const numbermodulo&  rhs)
	{
		lhs/=rhs;
		return lhs;
	}

	friend numbermodulo gcd(const numbermodulo& lhs,const numbermodulo&  rhs)
	{
		return 1;
	}

	friend numbermodulo lcm(const numbermodulo& lhs,const numbermodulo&  rhs)
	{
		return 1;
	}
};

#endif /* SRC_MODULUS_HPP_ */
