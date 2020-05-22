//============================================================================
// Name        : Helpers.cpp
// Author      : Jonathan Steinbuch
//============================================================================

#include "Helpers.hpp"

#include <limits>

#include <algorithm>


int dotproduct(const vector<unsigned int> &u, const vector <int> &v){
	assert(u.size()==v.size());
	int sum = 0;
	for(unsigned int i = 0; i < u.size(); i++){
		sum += u[i]*v[i];
	}
	return sum;
}

int partsum(const vector<unsigned int> &u, const vector <int> &v){
	int sum = 0;
	for(unsigned int i = 0; i < u.size(); i++){
		assert(u[i] < v.size());
		sum += v[u[i]];
	}
	return sum;
}

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



bool nextperm(vector<unsigned int> &perm, const unsigned int n, const unsigned int order)
{
	unsigned int endpos = perm.size() - 1;
	unsigned int i = 0;

	switch (order)
	{
	case degrevlex:
	{
		if (perm[endpos] == n)
			return false;
		unsigned int t = n;
		for (i = 0; i < endpos; i++)
		{
			t -= perm[i];
			if (perm[endpos] == t)
			{
				perm[i]--;
				perm[i + 1] = perm[endpos] + 1;
				if (i + 1 != endpos)
					perm[endpos] = 0;
				break;
			}
		}
		break;
	}
	case deglex:
	{
		if (perm[endpos] == n)
			return false;
		for (i = 0; i < endpos; i++)
		{
			if (perm[i] != 0)
			{
				perm[i]--;
				perm[i + 1] = perm[i + 1] + 1;
				if (i != 0)
				{
					perm[0] = perm[i];
					perm[i] = 0;
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

void kSumTon(vector<vector<unsigned int> >& list, const unsigned int n, const unsigned int k, const unsigned int order)
{
	vector<unsigned int> perm;
	perm.push_back(n);
	for (unsigned int i = 1; i < k; i++)
	{
		perm.push_back(0);
	}
	do
	{
		list.push_back(perm);
	} while (nextperm(perm, n, order));
}

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

void kSumTon(vector<vector<unsigned int> >& list, const unsigned int n, const unsigned int k, const unsigned int order, const vector<int> &weights, vector<unsigned int> &monomOffset)
{
	monomOffset = vector<unsigned int>(k-1,numeric_limits<unsigned int>::max());
	unsigned int start = list.size();

	vector<unsigned int> perm;
	unsigned int s = 0;
	perm.push_back(n);
	for (unsigned int i = 1; i < k; i++)
	{
		perm.push_back(0);
	}
	do
	{
		vector<unsigned int> wperm;
		if(applyWeights(perm,wperm,weights))
		{
			for(; s < k-1; s++)
			{
				if(perm[s] == 0)
				{
					monomOffset[s] = list.size()-start;
					start = list.size();
				}
				else
					break;
			}

			list.push_back(wperm);
		}
	} while ( nextperm(perm, n, order));
}

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

