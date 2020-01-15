//============================================================================
// Name        : Helpers.cpp
// Author      : Jonathan Steinbuch
//============================================================================

#include <limits>

#include <algorithm>

#include "Helpers.hpp"

int dotproduct(const vector<unsigned int> &u, const vector <int> &v){
	assert(u.size()==v.size());
	int sum = 0;
	for(unsigned int i = 0; i < u.size(); i++){
		sum += u[i]*v[i];
	}
	return sum;
}


//sum up only the elements in v that occur in u.
int partsum(const vector<unsigned int> &u, const vector <int> &v){
	int sum = 0;
	for(unsigned int i = 0; i < u.size(); i++){
		assert(u[i] < v.size());
		sum += v[u[i]];
	}
	return sum;
}

//sum up all elements of u
int sum(const vector<int> &u){
	int sum = 0;
	for(unsigned int i = 0; i < u.size(); i++){
		sum += u[i];
	}
	return sum;
}

int max(const vector<int> &input){
		int x = INT_MIN;
		for(auto & elem : input)
		{
			if(x < elem)
				x = elem;
		}
		return x;
	}

//Iterate over all possible ways to sum to n.
//Done in the way that translated to monomials it will represent the chosen monomial order.
bool nextsum(vector<unsigned int> &coeffs, const unsigned int n, const unsigned int order)
{
	unsigned int endpos = coeffs.size() - 1;
	unsigned int i = 0;

	switch (order)
	{
	case degrevlex:
	{
		if (coeffs[endpos] == n)
			return false;
		unsigned int t = n;
		for (i = 0; i < endpos; i++)
		{
			t -= coeffs[i];
			if (coeffs[endpos] == t)
			{
				coeffs[i]--;
				coeffs[i + 1] = coeffs[endpos] + 1;
				if (i + 1 != endpos)
					coeffs[endpos] = 0;
				break;
			}
		}
		break;
	}
	case deglex:
	{
		if (coeffs[endpos] == n)
			return false;
		for (i = 0; i < endpos; i++)
		{
			if (coeffs[i] != 0)
			{
				coeffs[i]--;
				coeffs[i + 1] = coeffs[i + 1] + 1;
				if (i != 0)
				{
					coeffs[0] = coeffs[i];
					coeffs[i] = 0;
				}
				break;
			}
		}
		break;
	}
	default:
	{
		return false;
		break;
	}
	}
	return true;
}

//Return the ways to sum to n with k elements, in the order corresponding to the given monomial order
void kSumTon(vector<vector<unsigned int> >& list, const unsigned int n, const unsigned int k, const unsigned int order)
{
	vector<unsigned int> coeffs;
	coeffs.push_back(n);
	for (unsigned int i = 1; i < k; i++)
	{
		coeffs.push_back(0);
	}
	do
	{
		list.push_back(coeffs);
	} while (nextsum(coeffs, n, order));
}

//Change monomials to properly reflect different weights of the exponents
bool applyWeights(const vector<unsigned int> &perm, vector<unsigned int> &wperm, const vector<int> &weights)
{
	wperm.clear();
	for(unsigned int i = 0; i < perm.size(); i++)
	{
		if(perm[i] % weights[i] != 0)
		{
			return false;
		} else {
			wperm.push_back(perm[i]/weights[i]);
		}
	}
	return true;
}

//Return the ways to sum to n with k elements, with the given weights
void kSumTon(vector<vector<unsigned int> >& list, const unsigned int n, const unsigned int k, const unsigned int order, const vector<int> &weights, vector<unsigned int> &monomOffset)
{
	monomOffset = vector<unsigned int>(k-1,numeric_limits<unsigned int>::max());
	unsigned int start = list.size();

	vector<unsigned int> coeffs;
	unsigned int s = 0;
	coeffs.push_back(n);
	for (unsigned int i = 1; i < k; i++)
	{
		coeffs.push_back(0);
	}
	do
	{
		vector<unsigned int> wcoeffs;
		if(applyWeights(coeffs,wcoeffs,weights))
		{
			for(; s < k-1; s++)
			{
				if(coeffs[s] == 0)
				{
					monomOffset[s] = list.size()-start;
					start = list.size();
				}
				else
					break;
			}

			list.push_back(wcoeffs);
		}
	} while ( nextsum(coeffs, n, order));
}

//Return all k element subsets of the set with n elements
void subsets(vector<vector<unsigned int> >& list, const unsigned int n, const unsigned int k, const unsigned int order)
{
	vector<vector<unsigned int> > list2;
	kSumTon(list2, n - k, k + 1, order);
	for (int i = list2.size() - 1; i >= 0; --i)
	{
		for (unsigned int j = 1; j < k; j++)
		{
			list2[i][j] += list2[i][j - 1] + 1;
		}
		list2[i].pop_back();
		list.push_back(list2[i]);
	}
}

//List all permutations of n elements
void permutations(vector<vector<unsigned int> > &list, const unsigned int n)
{
	vector<unsigned int> perm;
	for(unsigned int i = 0; i < n; i++)
	{
		perm.push_back(i);
	}
	list.push_back(perm);
	while(std::next_permutation(perm.begin(),perm.end()))
	{
		list.push_back(perm);
	}
}

//Return the parity of a given permutation
int parity(const vector<unsigned int> &permutation)
{
	vector<bool> visit(permutation.size(),false);
	unsigned int cycles = 0;
	for(unsigned int i = 0; i < permutation.size(); i++)
	{
		if(!visit[i])
		{
			for(unsigned int j = i; visit[j] != true; j = permutation[j])
			{
				visit[j] = true;
			}
			cycles++;
		}
	}
	if(cycles % 2 == 0){
		return 1;
	} else {
		return -1;
	}
}

//Return a  list of all permutations of n elements together with their parities
void permutationsandparities(vector<vector<unsigned int> > &list, const unsigned int n, vector<int> &paritylist)
{
	list.clear();
	paritylist.clear();

	permutations(list,n);
	for(unsigned int i = 0; i < list.size(); i++){
		int p = parity(list[i]);
		paritylist.push_back(p);
	}
}

