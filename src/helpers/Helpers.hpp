//============================================================================
// Name        : Helpers.hpp
// Author      : Jonathan Steinbuch
//============================================================================


#ifndef SRC_HELPERS_HELPERS_HPP_
#define SRC_HELPERS_HELPERS_HPP_


#include <iostream>
#include <vector>
#include <string>
#include <assert.h>

#include <cstdint>
#include <limits.h>


using namespace std;

struct program_options
{
	string config_file;
	string input_file;
	string output_file;
	int verbosity;
	bool exteriorPowers;
	bool stopIfUnstable;
	bool linearPowerProgresison;
	bool preAnalyze;
	bool outputM2;
	bool outputLatex;
	bool forceQ;
	bool computeFullKernel;
	bool longFormPolynomials;
};

const program_options noOptions = {"","","",1,false,false,false,false,false,false,false,false};

int dotproduct(const vector<unsigned int> &u, const vector <int> &v);

int partsum(const vector<unsigned int> &u, const vector <int> &v);

int sum(const vector<int> &u);

int max(const vector<int> &input);



class mBmonomial{
private:
	unsigned int id;

public:
	static const unsigned int mon_ERROR = UINT_MAX;

	bool operator>(const mBmonomial& other) const {
	    return id > other.id;
	}

	bool operator<=(const mBmonomial& other) const {
	    return id <= other.id;
	}

	bool operator>=(const mBmonomial& other) const {
	    return id >= other.id;
	}

	bool operator<(const mBmonomial& other) const {
	    return id < other.id;
	}

	bool operator==(const mBmonomial& other) const {
	    return id == other.id;
	}

	bool operator!=(const mBmonomial& other) const {
	    return id != other.id;
	}



	mBmonomial(const unsigned int id) : id(id){

	}

	unsigned int getId() const
	{
		return id;
	}

	bool isError() const
	{
		return id == mon_ERROR;
	}

	bool isOne() const
	{
		return id == 0;
	}
};

class lutable
{
	const unsigned int vars;

public:
	lutable(const unsigned int vars) : vars(vars)
	{
		computeuntildegree(vars+1);
	}

	void computeuntildegree(unsigned int degree)
	{
		unsigned int prevdeg = _choose.size();
		if(degree-prevdeg < prevdeg)
		{
			degree = 2*prevdeg;
		}
		_monN.resize(degree+2);
		_choose.resize(degree+1);
		_monomOffset.resize(degree+1);
		_monN[0] = 0;
		assert(vars < degree);
		for (unsigned int n = 0; n <= degree; ++n)
		{
			unsigned int kstart;
			if(n < prevdeg) {
				kstart = prevdeg;
			} else {
				kstart = 0;
			}
			_choose[n].resize(degree+1);
			for (unsigned int k = kstart; k <= degree; ++k)
			{
				if(n < 1 || k < 1){
					_choose[n][k] = 1;
				} else {
					_choose[n][k] = _choose[n][k-1] + _choose[n-1][k];
				}
			}
			_monN[n+1] = _monN[n] + _choose[n][vars - 1];
		}
		for (unsigned int nj = prevdeg; nj <= degree; ++nj)
		{
			for (unsigned int kminusj = 1; kminusj <= vars; ++kminusj)
			{
				unsigned int sum = 0;
				for (unsigned int i = 0; i < nj; i++)
				{
					sum += _choose[i][kminusj - 1];
				}
				_monomOffset[nj].push_back(sum);
			}
		}
	}

	unsigned int monomOffset(unsigned int degree,unsigned int vars)
	{
		if(_monomOffset.size() <= degree)
		{
			computeuntildegree(degree);
		}
		return _monomOffset[degree][vars];
	}

	unsigned int monN(unsigned int degree)
	{
		if(_monN.size() <= degree)
		{
			computeuntildegree(degree);
		}
		return _monN[degree];
	}

	unsigned int choose(unsigned int n, unsigned int k)
	{
		if (k > n)
		{
			return 0;
		}
		unsigned int nk = n - k;
		if (nk > k)
		{
			if (_choose.size() <= nk)
			{
				computeuntildegree(nk);

			}
		}
		else
		{
			if (_choose.size() <= k)
			{
				computeuntildegree(k);
			}
		}
		return _choose[nk][k];
	}


private:

	vector<vector<unsigned int> > _monomOffset;//[degree + 1][vars];
	vector<unsigned int > _monN;//[degree + 2];//monN[N] contains the number of monomials up to degree N-1
	vector<vector<unsigned int> > _choose;//[degree+1][degree+1]; //choose[n][k] contains n+k choose k, i.e. the number of monomials of degree n in k+1 variables
};

const unsigned int deglex = 0;
const unsigned int degrevlex = 1;
bool nextperm(vector<unsigned int> &perm, const unsigned int n, const unsigned int order);

void kSumTon(vector<vector<unsigned int> >& list, const unsigned int n, const unsigned int k, const unsigned int order);
bool applyWeights(const vector<unsigned int> &perm, vector<unsigned int> &wperm, const vector<int> &weights);
void kSumTon(vector<vector<unsigned int> >& list, const unsigned int n, const unsigned int k, const unsigned int order, const vector<int> &weights, vector<unsigned int> &monomOffset);

void subsets(vector<vector<unsigned int> >& list, const unsigned int n, const unsigned int k, const unsigned int order);

void permutations(vector<vector<unsigned int> > &list, const unsigned int n);
int parity(const vector<unsigned int> &permutation);
void permutationsandparities(vector<vector<unsigned int> > &list, const unsigned int n, vector<int> &paritylist);

#endif //SRC_LUTABLE_HPP
