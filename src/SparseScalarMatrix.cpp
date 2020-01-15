//============================================================================
// Name        : SparseScalarMatrix.cpp
// Author      : Jonathan Steinbuch
//============================================================================

#include "SparseScalarMatrix.hpp"

#include "SparsePolyMatrix.hpp"
#include "Modulus.hpp"

const mpq_class gcd(const mpq_class &u, const mpq_class &v)
{
	return mpq_class(gcd(u.get_num(), v.get_num()), gcd(u.get_den(), v.get_den()));
}

const mpq_class lcm(const mpq_class &u, const mpq_class &v)
{
	return mpq_class(lcm(u.get_num(), v.get_num()), lcm(u.get_den(), v.get_den()));
}

// add the value to the specified row and column
template<typename Scalar>
inline void SparseScalarMatrix<Scalar>::add(unsigned int r, unsigned int c, const Scalar &value)
{
	if (value == 0)
		return;

	auto ret = rows[r].insert(pair<unsigned int, Scalar>(c, value)); //try to insert
	if (!ret.second)
	{ //if it fails add the value to the existing entry
		(*(ret.first)).second += value;
		if ((*(ret.first)).second == 0)
		{ // if the sum is zero remove the zero entry
			rows[r].erase(ret.first);
			cols[c].erase(r);
		}
	}
	else
	{//if row insertion was successfull update the columns
		cols[c].insert(r);
	}
}

// moves the rows specified in sequence backwards
template<typename Scalar>
void SparseScalarMatrix<Scalar>::swapRows(const vector<unsigned int>& sequence)
{
	if (sequence.size() <= 1)
		return;
	spCol pnf, pf;
	for (unsigned int i = 0; i < sequence.size(); i++)
	{
		for (auto &ra : rows[sequence[i]])
		{
			unsigned int c = ra.first;
			if (i == sequence.size() - 1 && pnf.find(c) != pnf.end())
			{
				continue;
			}
			unsigned int nextrow = (sequence.size() + i - 1) % sequence.size();
			auto ret = cols[c].insert(sequence[nextrow]);
			if (i < sequence.size() - 1)
			{
				if (i == 0)
				{

					if (ret.second == false)
					{
						pf.insert(c);
					}
					else
					{
						pnf.insert(c);
					}
				}
				cols[c].erase(sequence[i]);
			}
			else
			{
				if (pf.find(c) == pf.end())
				{
					cols[c].erase(sequence[i]);
				}
			}
		}
	}

	spRow<Scalar> temp;
	temp = move(rows[sequence[0]]);
	for (unsigned int i = 0; i < sequence.size() - 1; i++)
	{
		rows[sequence[i]] = move(rows[sequence[i + 1]]);
	}
	rows[sequence[sequence.size() - 1]] = move(temp);

	//remember what we did by changing the row Permutation in the opposite order
	rowPerm.swap(sequence);
}

//swap two rows
template<typename Scalar>
void SparseScalarMatrix<Scalar>::swapRows(unsigned int a, unsigned int b)
{
	if (a == b)
		return;

	for (auto &ra : rows[a])
	{
		unsigned int c = ra.first;
		auto ret = cols[c].insert(b);
		if (ret.second)
		{
			cols[c].erase(a);
		}
	}
	for (auto &rb : rows[b])
	{
		unsigned int c = rb.first;
		auto ret = cols[c].insert(a);
		if (ret.second)
		{
			cols[c].erase(b);
		}
	}

	swap(rows[a], rows[b]);

	//remember what we did:
	rowPerm.swap(a, b);
}

//swap two columns
template<typename Scalar>
void SparseScalarMatrix<Scalar>::swapCols(unsigned int a, unsigned int b)
{
	if (a == b)
		return;

	for (auto &ca : cols[a])
	{
		auto ret = rows[ca].insert(pair<unsigned int, Scalar>(b, rows[ca][a]));
		if (ret.second)
		{ //insert succeded
			rows[ca].erase(a);
		}
		else
		{
			Scalar temp = rows[ca][b];
			rows[ca][b] = rows[ca][a];
			rows[ca][a] = temp;
		}

	}
	for (auto &cb : cols[b])
	{
		auto ret = rows[cb].insert(pair<unsigned int, Scalar>(a, rows[cb][b]));
		if (ret.second)
		{
			rows[cb].erase(b);
		} // no else, since if insertion failed the values were already exchanged in the first loop.
	}

	swap(cols[a], cols[b]);

	//remember what we did:
	colPerm.swap(a, b);
}

//Symmetrical swaps of rows and columns at the same time
template<typename Scalar>
void SparseScalarMatrix<Scalar>::symmSwaps(const vector<unsigned int>& newpos)
{
	vector<spRow<Scalar>> oldRows = rows;
	vector<spCol> oldCols = cols;

	for (unsigned int i = 0; i < cols.size(); i++)
	{
		spCol newCol;
		for (auto &x : oldCols[i])
		{
			if (x < newpos.size())
			{
				newCol.insert(newpos[x]);
			}
			else
			{
				newCol.insert(x);
			}
		}

		if (i < newpos.size())
		{
			cols[newpos[i]] = newCol;
		}
		else
		{
			cols[i] = newCol;
		}
	}

	for (unsigned int i = 0; i < rows.size(); i++)
	{
		spRow<Scalar> newRow;
		for (auto &x : oldRows[i])
		{
			if (x.first < newpos.size())
			{
				newRow.insert(pair<unsigned int, Scalar>(newpos[x.first], x.second));
			}
			else
			{
				newRow.insert(pair<unsigned int, Scalar>(x.first, x.second));
			}
		}
		if (i < newpos.size())
		{
			rows[newpos[i]] = newRow;
		}
		else
		{
			rows[i] = newRow;
		}

	}

	rowPerm = Permutation(rowPerm.size(), newpos) * rowPerm;
	colPerm = Permutation(colPerm.size(), newpos) * colPerm;

}

template<typename Scalar>
const Scalar SparseScalarMatrix<Scalar>::randomreadaccess(unsigned int row, unsigned int column) const
{
	auto it = rows[row].find(column);
	if (it != rows[row].end())
	{
		return (*it).second;
	}
	else
	{
		return 0;
	}

}

//add two rows, with scalar multiplication i.e. targetrow = targetrow*factortarget + otherrow*factorother.
template<typename Scalar>
inline void SparseScalarMatrix<Scalar>::addrows(const unsigned int targetrow, const unsigned int otherrow, const Scalar &factortarget, const Scalar &factorother)
{
	assert(factortarget != 0 && factorother != 0);

	if (factortarget != 1)
	{
		for (auto &x : rows[targetrow])
		{
			x.second *= factortarget;
		}
	}

	for (auto &x : rows[otherrow])
	{
		add(targetrow, x.first, x.second * factorother);
	}
}

//alternative but slower and thus useless implementation
template<typename Scalar>
void SparseScalarMatrix<Scalar>::addrowsslow(const unsigned int targetrow, const unsigned int otherrow, const Scalar &factortarget, const Scalar &factorother)
{
	assert(factortarget != 0 && factorother != 0);

	auto targetIt = rows[targetrow].cbegin();
	auto otherIt = rows[otherrow].cbegin();

	auto targetEnd = rows[targetrow].cend();
	auto otherEnd = rows[otherrow].cend();

	vector<pair<unsigned int, Scalar> > temp;

	while (targetIt != targetEnd || otherIt != otherEnd)
	{
		if (targetIt == targetEnd || (otherIt != otherEnd && (*targetIt).first > (*otherIt).first))
		{
			Scalar newVal = (*otherIt).second * factorother;
			temp.push_back(pair<unsigned int, Scalar>((*otherIt).first, newVal));
			auto ret = cols[(*otherIt).first].insert(targetrow);
			assert(ret.second);
			++otherIt;
		}
		else if (otherIt == otherEnd || (*targetIt).first < (*otherIt).first)
		{
			Scalar newVal = (*targetIt).second * factortarget;
			temp.push_back(pair<unsigned int, Scalar>((*targetIt).first, newVal));
			++targetIt;
		}
		else
		{
			Scalar newVal = (*targetIt).second * factortarget + (*otherIt).second * factorother;
			if (newVal != 0)
			{
				temp.push_back(pair<unsigned int, Scalar>((*targetIt).first, newVal));
			}
			else
			{
				size_t rem = cols[(*targetIt).first].erase(targetrow);
				assert(rem == 1);
			}
			++otherIt;
			++targetIt;
		}
	}

	rows[targetrow] = spRow<Scalar>(boost::container::ordered_unique_range_t(), temp.begin(), temp.end());
}

//add columns together
template<typename Scalar>
void SparseScalarMatrix<Scalar>::addcols(unsigned int targetcol, unsigned int othercol, Scalar factortarget, Scalar factorother)
{
	assert(factortarget != 0 && factorother != 0);
	for (auto &x : cols[targetcol])
	{
		auto rit = rows[x].find(targetcol);
		(*rit).second *= factortarget;
	}

	for (auto &x : cols[othercol])
	{
		auto rit = rows[x].find(othercol);
		add(x, targetcol, (*rit).second * factorother);
	}
}

//simplify a row by dividing out common scalar factors.
template<typename Scalar>
inline void SparseScalarMatrix<Scalar>::simplifyrow(unsigned int row)
{
	Scalar factor = 0;
	for (auto &x : rows[row])
	{
		factor = gcd(x.second, factor);
		if (factor == 1 || factor == -1)
		{
			return;
		}
	}
	for (auto &x : rows[row])
	{
		x.second /= factor;
	}
}

//returns the number of nonzero entries in a row
template<typename Scalar>
inline unsigned int SparseScalarMatrix<Scalar>::rowcount(unsigned int row) const
{
	return rows[row].size();
}

//returns the number of nonzero entries in a row in the column range given by range
template<typename Scalar>
inline unsigned int SparseScalarMatrix<Scalar>::rowcount(unsigned int row, const pair<unsigned int, unsigned int>& range) const
{
	unsigned int ret = 0;
	auto rowit = rows[row].lower_bound(range.first);
	while (rowit != rows[row].end() && (*rowit).first < range.second)
	{
		rowit++;
		ret++;
	}
	return ret;
}

//Triangularization of the matrix using the Gauss Algorithm
template<typename Scalar>
void SparseScalarMatrix<Scalar>::reduce(unsigned int &defect, const program_options &opt)
{

	unsigned int progress = 0;
	if (opt.verbosity >= 1)
	{
		cout << "Progress";
	}

	unsigned int maxr = min(cols.size(), rows.size());
	for (unsigned int j = 0; j < maxr; j++)
	{
		unsigned int prow = j;
		unsigned int pcount = UINT_MAX;
		auto colit = cols[j].lower_bound(j);
		while (colit != cols[j].end())
		{
			if (rowcount(*colit) < pcount)
			{
				prow = *colit;
				pcount = rowcount(prow);
			}
			colit++;
		}
		if (pcount == UINT_MAX)
		{
			auto rowit = rows[j].begin();
			if (rowit != rows[j].end())
			{

				swapCols(j, (*rowit).first);
			}
			else
			{
				maxr--;

				swapRows(j, maxr);

				swapCols(j, maxr);
			}
			j--;
			continue;
		}
		else
		{
			if (j != prow)
			{
				swapRows(j, prow);
			}

			simplifyrow(j);

			colit = cols[j].upper_bound(j);
			Scalar pivot = (*(rows[j].find(j))).second;
			while (colit != cols[j].end())
			{
				unsigned int target = *colit;
				Scalar val = (*(rows[target].find(j))).second;

				Scalar g = gcd(pivot, val);

				Scalar factor1 = pivot / g;
				Scalar factor2 = -val / g;

				addrows(target, j, factor1, factor2);

				colit = cols[j].upper_bound(j);
			}
		}

		if (opt.verbosity >= 1)
		{
			while (maxr > 0 && j * 100 > progress * maxr)
			{
				if (progress % 10 == 0)
					cout << ":" << flush;
				else
					cout << "." << flush;
				progress++;
			}
		}
	}
	if (opt.verbosity >= 1)
	{
		cout << "!" << endl;
	}

	defect = cols.size() - maxr;
}

//Tarjan's algorithm for turning the matrix into a block diagonal matrix by only looking at the shape of the nonzeros not the actual values
template<typename Scalar>
void SparseScalarMatrix<Scalar>::tarjan(unsigned int maxr, vector<pair<unsigned int, unsigned int> >& blocks)
{
	assert(maxr <= rows.size() && maxr <= cols.size());
	if (maxr == 0)
	{
		return;
	}
	blocks.clear();
	vector<unsigned int> stack;
	vector<unsigned int> pathsFoundTo(maxr, UINT_MAX);
	vector<unsigned int> previousNode(maxr, UINT_MAX);
	vector<unsigned int> posOnStack(maxr, UINT_MAX);
	vector<bool> final(maxr, false);
	unsigned int topPos = maxr;
	unsigned int bottomPos = 0;

	vector<typename spRow<Scalar>::const_iterator> progressInRow;
	progressInRow.reserve(maxr);

	for (unsigned int i = 0; i < maxr; i++)
	{
		progressInRow.push_back(rows[i].begin());
	}

	unsigned int node = 0;
	unsigned int reached = 0;
	bool pathstart = true;

	while (reached < maxr || stack.size() > 0)
	{
		typename spRow<Scalar>::const_iterator nextit = progressInRow[node];
		if (nextit == rows[node].end() || (*nextit).first >= maxr)
		{ //we reached an end;
			if (pathsFoundTo[node] != UINT_MAX && posOnStack[pathsFoundTo[node]] < posOnStack[node])
			{
				assert(previousNode[node] != UINT_MAX);
				if (posOnStack[pathsFoundTo[node]] < posOnStack[pathsFoundTo[previousNode[node]]])
				{
					assert(pathsFoundTo[node] != UINT_MAX);
					pathsFoundTo[previousNode[node]] = pathsFoundTo[node];
				}

				node = previousNode[node];
			}
			else
			{
				if (posOnStack[node] == UINT_MAX)
				{
					topPos--;
					posOnStack[node] = topPos;
					final[node] = true;
				}
				else
				{
					unsigned int pNode = posOnStack[node];
					assert(pNode < stack.size());
					unsigned int blocksize = stack.size() - pNode;
					if (blocksize > 1)
					{
						assert(topPos - blocksize >= 0);
						blocks.push_back(pair<unsigned int, unsigned int>(topPos - blocksize, topPos));
						if (blocks.size() > 1)
							assert(blocks[blocks.size() - 1].second <= blocks[blocks.size() - 2].first);
					}
					for (unsigned int i = pNode; i < stack.size(); i++)
					{
						topPos--;
						posOnStack[stack[i]] = topPos;
						final[stack[i]] = true;
					}
					stack.resize(pNode);
				}

				if (stack.size() == 0)
				{
					for (; reached < maxr; reached++)
					{
						if (!final[reached] && posOnStack[reached] == UINT_MAX)
						{
							node = reached;
							pathstart = true;
							break;
						}
					}
				}
				else
				{
					assert(stack[stack.size() - 1] < maxr); //fail at power 20
					node = stack[stack.size() - 1];
				}
			}
		}
		else
		{
			unsigned int nextnode = (*nextit).first;
			progressInRow[node]++;
			assert(nextnode < maxr);
			if (posOnStack[nextnode] != UINT_MAX)
			{
				if (!final[nextnode])
				{
					if (posOnStack[nextnode] < posOnStack[pathsFoundTo[node]])
					{
						assert(pathsFoundTo[nextnode] != UINT_MAX);
						pathsFoundTo[node] = pathsFoundTo[nextnode];
					}
				}
			}
			else
			{
				posOnStack[nextnode] = stack.size();
				stack.push_back(nextnode);
				pathsFoundTo[nextnode] = nextnode;
				if (pathstart == false)
				{
					assert(node != nextnode);
					previousNode[nextnode] = node;
				}		//TODO: remove Debug Checks

				pathstart = false;
				node = nextnode;
			}
		}
	}
	assert(bottomPos == topPos);
	symmSwaps(posOnStack);
}

//Check whether the columns and rows still agree on which entries are nonzero
template<typename Scalar>
bool SparseScalarMatrix<Scalar>::checkintegrity()
{
	for (unsigned int i = 0; i < cols.size(); i++)
	{
		for (auto &ca : cols[i])
		{
			if (rows[ca].find(i) == rows[ca].end())
				return false;
		}
	}

	for (unsigned int i = 0; i < rows.size(); i++)
	{
		for (auto &ra : rows[i])
		{
			if (cols[ra.first].find(i) == cols[ra.first].end())
				return false;
		}
	}
	return true;
}

//Move the rows in an order which could be beneficial for performing the Gauss algorithm
template<typename Scalar>
unsigned int SparseScalarMatrix<Scalar>::analyze(const program_options & opt)
{
	unsigned int symbdefect = 0;

	auto thiscopy = *this;

	vector<unsigned int> clookaheadprogress(cols.size(), 0);

	for (unsigned int c = 0; c < cols.size() - symbdefect;)
	{
		vector<unsigned int> sequence;
		sequence.push_back(c);
		while (sequence.size() >= 1)
		{
			unsigned int last = *sequence.rbegin();
			spCol::iterator c_it1 = cols[last].lower_bound(c);
			if (c_it1 == cols[last].end())
			{
				spCol::iterator c_it2 = cols[last].lower_bound(clookaheadprogress[last]);
				if (c_it2 == cols[last].end())
				{
					clookaheadprogress[last] = rows.size();
					sequence.pop_back();
					if (sequence.size() > 0)
						clookaheadprogress[*sequence.rbegin()] = last + 1;
				}
				else
				{
					bool otf = false;
					for (unsigned int i = 0; i < sequence.size(); i++)
					{
						if (*c_it2 == sequence[i])
						{
							otf = true;
							break;
						}
					}
					if (otf)
					{
						clookaheadprogress[last] = *c_it2 + 1;
					}
					else
					{
						sequence.push_back(*c_it2);
					}

				}
			}
			else
			{

				unsigned int prow = *c_it1;
				unsigned int pcount = rows[prow].size();
				while (c_it1 != cols[last].end())
				{
					if (rows[*c_it1].size() < pcount)
					{
						prow = *c_it1;
						pcount = rows[prow].size();
					}
					*c_it1++;
				}
				if (prow != c)
				{
					sequence.push_back(prow);
				}

				swapRows(sequence);
				break;
			}
		}
		if (sequence.size() == 0)
		{
			symbdefect++;
			swapCols(c, cols.size() - symbdefect);
		}
		else
		{
			c++;
		}
	}

	vector<pair<unsigned int, unsigned int> > blocks;
	unsigned int maxr = cols.size() - symbdefect;

	if (rows.size() > maxr && cols.size() > maxr)
	{
		maxr = min(rows.size(), cols.size());
	}

	tarjan(maxr, blocks);

	return symbdefect;
}

template<typename Scalar>
Scalar SparseScalarMatrix<Scalar>::rowproduct(const spRow<Scalar>& v1, const spRow<Scalar>& v2)
{
	Scalar ret = 0;
	auto it1 = v1.begin();
	auto it2 = v2.begin();
	while (it1 != v1.end() && it2 != v2.end())
	{
		if ((*it1).first > (*it2).first)
		{
			it2++;
		}
		else if ((*it1).first < (*it2).first)
		{
			it1++;
		}
		else
		{
			ret += (*it1).second * (*it2).second;
			it1++;
			it2++;
		}
	}
	return ret;
}

template<typename Scalar>
void SparseScalarMatrix<Scalar>::multiply(spRow<Scalar>& v, Scalar factor)
{
	for (auto &x : v)
	{
		x.second *= factor;
	}
}

//return the matrix only consisting of the columns selected
template<typename Scalar>
SparseScalarMatrix<Scalar> SparseScalarMatrix<Scalar>::colSelect(vector<unsigned int> selection)
{
	SparseScalarMatrix<Scalar> out;

	for (auto &i : selection)
	{
		out.cols.push_back(cols[i]);
	}
	for (auto &j : rows)
	{
		spRow<Scalar> row;
		for (unsigned int i = 0; i < selection.size(); i++)
		{
			auto a = j.find(selection[i]);
			if (a != j.end())
			{
				row.insert(pair<unsigned int, Scalar>(i, a->second));
			}
		}
		out.rows.push_back(row);
	}
	out.rowPerm = Permutation(out.rows.size());
	out.colPerm = Permutation(out.cols.size());
	assert(out.checkintegrity());
	return out;
}

//Move a row and a column to the end of the new matrix.
template<typename Scalar>
SparseScalarMatrix<Scalar> SparseScalarMatrix<Scalar>::colUnSelect(unsigned int unselectionCol, unsigned int unselectionRow)
{
	SparseScalarMatrix<Scalar> out;

	for (unsigned int i = 0; i < cols.size(); i++)
	{
		if (i != unselectionCol)
		{
			spCol col;
			for (auto &a : cols[i])
			{
				if (a < unselectionRow)
				{
					col.insert(a);
				}
				else if (a > unselectionRow)
				{
					col.insert(a - 1);
				}
			}
			out.cols.push_back(col);
		}
	}
	for (unsigned int j = 0; j < rows.size(); j++)
	{
		if (j != unselectionRow)
		{
			spRow<Scalar> row;
			for (auto &a : rows[j])
			{
				if (a.first < unselectionCol)
				{
					row.insert(pair<unsigned int, Scalar>(a.first, a.second));
				}
				else if (a.first > unselectionCol)
				{
					row.insert(pair<unsigned int, Scalar>(a.first - 1, a.second));
				}
			}
			out.rows.push_back(row);
		}
	}
	out.rowPerm = Permutation(out.rows.size());
	out.colPerm = Permutation(out.cols.size());
	assert(out.checkintegrity());
	return out;
}

//Experimental modulo operator, taking the image of the columns less than lhscols mnodulo the image of the right hand side columns.
template<typename Scalar>
SparseScalarMatrix<Scalar> SparseScalarMatrix<Scalar>::modulo(unsigned int lhscols)
{
	SparseScalarMatrix<Scalar> copy = *this;

	multimap<unsigned int, unsigned int> remCols;
	for (unsigned int i = 0; i < lhscols; i++)
	{
		remCols.insert(pair<unsigned int, unsigned int>(copy.cols[i].size(), i));
	}

	vector<unsigned int> goodcols;

	for (auto rit = remCols.rbegin(); rit != remCols.rend(); rit++)
	{
		goodcols.push_back((*rit).second);
	}

	SparseScalarMatrix<Scalar> ret = copy.colSelect(goodcols);

	//cout << copy << endl << endl;

	for (unsigned int i = lhscols; i < copy.getCols(); i++)
	{
		goodcols.push_back(i);
	}
	copy = copy.colSelect(goodcols);
	//cout << copy << endl << endl;

	SparseScalarMatrix<Scalar> kernel = copy.kernel();

	//cout << kernel << endl << endl;

	for (unsigned int i = 0; i < ret.getCols();)
	{
		auto x = kernel.rows[i].begin();
		if (x != kernel.rows[i].end())
		{
			auto y = x;
			y++;
			for (; y != kernel.rows[i].end();)
			{
				unsigned int tcol = (*y).first;
				Scalar g = gcd((*x).second, (*y).second);
				Scalar tfactor = (*x).second / g;
				Scalar ofactor = -(*y).second / g;
				y++;
				kernel.addcols(tcol, (*x).first, tfactor, ofactor);
			}
			ret = ret.colUnSelect(i);
			kernel = kernel.colUnSelect((*x).first, i);
			//	cout << kernel << endl;
			i = 0;
		}
		else
		{
			i++;
		}
	}
	return ret;
}

//Compute the kernel of a sparse matrix by way of triangularization
template<typename Scalar>
SparseScalarMatrix<Scalar> SparseScalarMatrix<Scalar>::kernel(const program_options & opt)
{
	SparseScalarMatrix<Scalar> oldthis = *this;
	unsigned int defect = 0;
	if (opt.useTarjan)
	{
		defect = analyze(opt);
	}

	reduce(defect, opt);

	vector<spRow<Scalar>> nrows;

	unsigned int maxcol = cols.size();
	if (!opt.computeFullKernel)
	{
		maxcol = min(cols.size(), cols.size() - defect + 1);
	}

	for (unsigned int i = cols.size() - defect; i < maxcol; i++)
	{
		spRow<Scalar> nrow;
		nrow.insert(pair<unsigned int, Scalar>(i, 1)); //1-entry in the generator
		for (unsigned int j = cols.size() - defect; j > 0; j--)
		{
			Scalar result = rowproduct(rows[j - 1], nrow);
			if (result != 0)
			{
				Scalar diag = (*rows[j - 1].begin()).second;
				Scalar g = gcd(result, diag);
				multiply(nrow, diag / g);

				Scalar nres = rowproduct(rows[j - 1], nrow);
				Scalar mult = -nres / diag;

				nrow.insert(pair<unsigned int, Scalar>(j - 1, mult));

				assert(rowproduct(rows[j - 1], nrow) == 0);
			}
		}
		nrows.push_back(nrow);
	}
	if (!opt.computeFullKernel && opt.verbosity >= 1 && defect > 1)
	{
		cout << "Stopped kernel computation. Actual kernel dimension: " << defect << endl;
	}

	SparseScalarMatrix<Scalar> kernel(nrows, cols.size(), false);
	auto cpm = colPerm.matrix<Scalar>(true);
	kernel = cpm * kernel; //unpermute the rows of the kernel

	return kernel;
}

template<typename Scalar>
SparseScalarMatrix<Scalar>::SparseScalarMatrix(const vector<rcvTriplet<Scalar> > &triplets, unsigned int nrows, unsigned int ncols)
{
	rows = vector<spRow<Scalar> >(nrows);
	cols = vector<spCol>(ncols);
	for (const rcvTriplet<Scalar> &trip : triplets)
	{
		add(trip.row, trip.col, trip.value);
	}
	rowPerm = Permutation(nrows);
	colPerm = Permutation(ncols);
	outrows.clear();
}

template<typename Scalar>
SparseScalarMatrix<Scalar>::SparseScalarMatrix(const SparsePolyMatrix<Scalar> &inputMatrix, const int degree, bool homogeneous)
{
	unsigned int trows, tcols;
	vector<rcvTriplet<Scalar> > tripletList;
	inputMatrix.degreeTriplets(degree, tripletList, trows, tcols, homogeneous);

	rows = vector<spRow<Scalar> >(trows);
	cols = vector<spCol>(tcols);

	for (const rcvTriplet<Scalar> &trip : tripletList)
	{
		add(trip.row, trip.col, trip.value);
	}
	rowPerm = Permutation(trows);
	colPerm = Permutation(tcols);
	outrows.clear();
}

template<typename Scalar>
SparseScalarMatrix<Scalar>::SparseScalarMatrix(const vector<spRow<Scalar> > &rows, const vector<spCol> &cols) :
		rows(rows), cols(cols)
{
	rowPerm = Permutation(rows.size());
	colPerm = Permutation(cols.size());
	assert(checkintegrity());
	outrows.clear();
}

template<typename Scalar>
void SparseScalarMatrix<Scalar>::fillcolumns()
{
	for (unsigned int i = 0; i < rows.size(); i++)
	{
		for (auto &ra : rows[i])
		{
			cols[ra.first].insert(i);
		}
	}
}

template<typename Scalar>
SparseScalarMatrix<Scalar>::SparseScalarMatrix(const vector<spRow<Scalar>>& vecs, unsigned int ncols, bool rowMajor)
{
	if (rowMajor)
	{
		rows = vecs;
		cols = vector<spCol>(ncols);
		fillcolumns();
	}
	else
	{
		rows = vector<spRow<Scalar>>(ncols);
		for (unsigned int i = 0; i < vecs.size(); i++)
		{
			for (auto &ra : vecs[i])
			{
				rows[ra.first].insert(pair<unsigned int, Scalar>(i, ra.second));
			}
		}
		cols = vector<spCol>(vecs.size());
		fillcolumns();
	}

	rowPerm = Permutation(rows.size());
	colPerm = Permutation(cols.size());
	outrows.clear();

	checkintegrity();
}

template<typename Scalar>
SparseScalarMatrix<Scalar>& SparseScalarMatrix<Scalar>::operator*=(const SparseScalarMatrix<Scalar>& rhs)
{
	*this = (*this) * rhs;
	return *this;
}

template<typename Scalar>
bool SparseScalarMatrix<Scalar>::operator==(const SparseScalarMatrix<Scalar>& rhs) const
{
	if (rows.size() != rhs.rows.size() || cols.size() != rhs.cols.size())
		return false;
	for (unsigned int i = 0; i < rows.size(); i++)
	{
		for (auto &x : rows[i])
		{
			auto ret = rhs.rows[i].find(x.first);
			if (ret == rhs.rows[i].end() || (*ret).second != x.second)
				return false;
		}
	}
	return true;
}

template<typename Scalar>
SparseScalarMatrix<Scalar> SparseScalarMatrix<Scalar>::operator*(const SparseScalarMatrix<Scalar>& rhs) const
{
	const SparseScalarMatrix<Scalar> &lhs = *this;
	assert(lhs.cols.size() == rhs.rows.size());
	vector<spRow<Scalar>> newrows;
	newrows.reserve(lhs.rows.size());
	for (unsigned int i = 0; i < lhs.rows.size(); i++)
	{
		spRow<Scalar> row;
		for (auto &x : lhs.rows[i])
		{
			for (auto &y : rhs.rows[x.first])
			{
				Scalar product = x.second * y.second;
				auto ret = row.insert(pair<unsigned int, Scalar>(y.first, product));
				if (ret.second == false)
				{
					(*ret.first).second += product;
				}
			}
		}
		for (auto it = row.begin(); it != row.end();)
		{
			if ((*it).second == 0)
			{
				auto f = (*it).first;
				row.erase(it);
				it = row.lower_bound(f);
			}
			else
			{
				it++;
			}
		}
		newrows.push_back(row);
	}
	return SparseScalarMatrix(newrows, rhs.cols.size());
}

// General functions


//Only print the block sx to ex times sy to ey of the matrix
template<typename Scalar>
void SparseScalarMatrix<Scalar>::printblock(unsigned int sx, unsigned int sy, unsigned int ex, unsigned int ey, std::ostream &os)
{
	for (unsigned int y = sy; y < ey; y++)
	{
		auto rowit = rows[y].lower_bound(sx);
		for (unsigned int x = sx; x < ex; x++)
		{
			if (rowit == rows[y].end() || (*rowit).first > x)
			{
				os << 0;
			}
			else
			{
				os << (*rowit).second;
				rowit++;
			}
			os << "\t";
		}
		os << endl;
	}
	os << endl;
}

template<typename _Scalar>
std::ostream &operator<<(std::ostream &os, SparseScalarMatrix<_Scalar> const &m)
{
	if (m.cols.size() == 0)
	{
		return os;
	}
	if (m.outrows.size() != 0)
	{
		for (auto &x : m.outrows)
		{
			os << x << "\t|";
			for (auto &j : m.rows[x])
			{
				os << j.first << "\t";
			}
			os << endl;
		}
		return os;
	}
	else
	{
		for (unsigned int i = 0; i < m.rows.size(); i++)
		{
			os << i << "\t|";
			for (auto &j : m.rows[i])
			{
				os << j.first << "\t";
			}
			os << endl;
		}
		return os;
	}

	for (unsigned int i = 0; i < m.rows.size(); i++)
	{
		for (unsigned int j = 0; j < m.cols.size(); j++)
		{
			os << m.randomreadaccess(i, j) << "\t";
		}
		os << endl;
	}
	return os;
}
template std::ostream &operator<<(std::ostream &os, SparseScalarMatrix<mpz_class> const &m);
template std::ostream &operator<<(std::ostream &os, SparseScalarMatrix<mpq_class> const &m);
template std::ostream &operator<<(std::ostream &os, SparseScalarMatrix<numbermodulo> const &m);

template class SparseScalarMatrix<mpz_class> ;
template class SparseScalarMatrix<mpq_class> ;
template class SparseScalarMatrix<numbermodulo> ;
