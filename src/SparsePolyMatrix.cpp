//============================================================================
// Name        : SparsePolyMatrix.cpp
// Author      : Jonathan Steinbuch
//============================================================================

#include "SparsePolyMatrix.hpp"

#include "Helpers.hpp"
#include "Modulus.hpp"


template<typename _Scalar>
	int SparsePolyMatrix<_Scalar>::indegree()
	{
		return sum(indegrees);
	}

template<typename _Scalar>
	int SparsePolyMatrix<_Scalar>::outdegree()
	{
		return sum(outdegrees);
	}

//Tries to compute all twists from the given output twists.
//If there are no outdegrees given, the first outdegree is assumed to be 0.
//Will return false if the degree computation is not possible.
template<typename _Scalar>
bool SparsePolyMatrix<_Scalar>::computeDegrees(const vector<int> &outdegs)
{
	vector<bool> inset;
	vector<bool> outset;

	outdegrees.clear();
	indegrees.clear();

	if (outdegs.size() == rows.size())
	{
		outdegrees = outdegs;
		outset.insert(outset.begin(), rows.size(), true);
	}
	else
	{
		if (outdegs.size() == 0)
		{
			outdegrees.push_back(0);
		}
		else if (outdegs.size() == 1)
		{
			outdegrees.push_back(outdegs[0]);
		}
		else
		{
			return false;
		}
		outdegrees.insert(outdegrees.begin() + 1, rows.size() - 1, 0);
		outset.push_back(true);
		outset.insert(outset.begin() + 1, rows.size() - 1, false);
	}

	indegrees.insert(indegrees.begin(), cols.size(), 0);
	inset.insert(inset.begin(), cols.size(), false);

	for (unsigned int k = 0; k < rows.size(); k++)
	{
		for (unsigned int i = 0; i < rows.size(); i++)
		{
			for (auto & entry : rows[i])
			{
				if (!entry.second.isHomogeneous())
					return false;
				if (outset[i])
				{
					int ndeg = -entry.second.degree() + outdegrees[i];
					if (inset[entry.first] && indegrees[entry.first] != ndeg)
					{
						cerr << "Degree computation error!"<< endl;
						return false;

					}
					else
					{
						indegrees[entry.first] = ndeg;
						inset[entry.first] = true;
					}
				}
				else
				{
					if (inset[entry.first])
					{
						outdegrees[i] = indegrees[entry.first] + entry.second.degree();
						outset[i] = true;
					}
				}
			}
		}

		bool set = true;
		for (unsigned int j = 0; j < rows.size(); j++)
		{
			if (!outset[j])
			{
				set = false;
				break;
			}
		}

		for (unsigned int j = 0; j < cols.size(); j++)
		{
			if (!inset[j])
			{
				set = false;
				break;
			}
		}

		if (set == true)
			break;
	}

	for (unsigned int j = 0; j < rows.size(); j++)
	{
		if (!outset[j])
			return false;
	}

	for (unsigned int j = 0; j < cols.size(); j++)
	{
		if (!inset[j])
			return false;
	}

	return true;
}

template <typename _Scalar>
SparsePolyMatrix<_Scalar>::SparsePolyMatrix(const DensePolyMatrix<_Scalar> &old)
{
	Ring = old.Ring;

	indegrees =old.indegrees;
	outdegrees = old.outdegrees;

	inlabels = old.inlabels;
	outlabels = old.outlabels;

	rows.clear();
	cols.clear();

	for (unsigned int c = 0; c < old.cols; c++)
	{
		spPCol col;
		for (unsigned int r = 0; r < old.rows; r++)
		{
			if(!((old.entries[r][c]).isNull()))
			{
				col.insert(r);
			}
		}
		cols.push_back(col);
	}
	for (unsigned int r = 0; r < old.rows; r++)
	{
		spPRow<_Scalar> row;
		for (unsigned int c = 0; c < old.cols; c++)
		{
			if(!old.entries[r][c].isNull())
			{
				row.insert(pair<unsigned int,Poly<_Scalar> >(c,old.entries[r][c]));
			}
		}
		rows.push_back(row);
	}
}

template <typename _Scalar>
SparsePolyMatrix<_Scalar>::SparsePolyMatrix(PolyRing<_Scalar>* mB, const vector<spPRow<_Scalar> > &rows,const vector<spPCol > &cols, const vector<int> &outdegs) :
 Ring(mB) ,rows(rows), cols(cols)
	{
		computeDegrees(outdegs);
	}

template <typename _Scalar>
SparsePolyMatrix<_Scalar>::SparsePolyMatrix(PolyRing<_Scalar>* mB, const vector<spPRow<_Scalar> > &rows,const vector<spPCol > &cols, const vector<int> &outdegs, const vector<int> &indegs) :
 Ring(mB), indegrees(indegs), outdegrees(outdegs),	rows(rows), cols(cols)
	{

	}

template <typename _Scalar>
void SparsePolyMatrix<_Scalar>::computeCols(unsigned int numcols)
{
	cols.clear();
	cols.resize(numcols);
	for(unsigned int r = 0; r < rows.size(); r++){
		for(auto & elem : rows[r]){
			assert(elem.first < numcols);
			cols[elem.first].insert(r);
		}
	}
}

template <typename _Scalar>
SparsePolyMatrix<_Scalar>::SparsePolyMatrix(PolyRing<_Scalar>* mB,const vector<spPRow<_Scalar> > &rows,unsigned int numcols, const vector <int> &outdegs) :
Ring(mB) ,rows(rows)
	{
		computeCols(numcols);
		computeDegrees(outdegs);
	}

template <typename _Scalar>
SparsePolyMatrix<_Scalar>::SparsePolyMatrix(PolyRing<_Scalar>* mB,const vector<spPRow<_Scalar> > &rows,unsigned int numcols,const vector <int> &outdegs,const vector <int> &indegs) :
Ring(mB), indegrees(indegs), outdegrees(outdegs),	rows(rows)
	{
		computeCols(numcols);
	}

//Create polynomial matrix from the scalar representation.
//sx and sy are the dimensions of the output matrix.
//The degree is the applied twist to which the scalar matrix was formed.
//Make sure the input matrix is well-formed.
template<typename _Scalar>
SparsePolyMatrix<_Scalar>::SparsePolyMatrix(PolyRing<_Scalar>* mB, const SparseScalarMatrix<_Scalar>& inmatrix, const unsigned int sx, const unsigned int sy, const vector<int> &outdegs, const int degree, bool homogeneous) :
		Ring(mB)
{
	unsigned int rowoffsets[sx + 1] =
	{ 0 };

	assert(outdegs.size() == sx);

	cols.resize(sy);

	if (homogeneous)
	{

		for (unsigned int i = 0; i < sx; i++)
		{
			int poldeg = degree + outdegs[i];
			unsigned int bSize = mB->getBasisSize(poldeg);
			rowoffsets[i + 1] = rowoffsets[i] + bSize;
			spPRow<_Scalar> row;
			for (unsigned int c = 0; c < sy; c++)
			{
				Poly<_Scalar> poly(mB);
				bool empty = true;
				for (unsigned int j = 0; j < bSize; j++)
				{

					_Scalar val = inmatrix.randomreadaccess(rowoffsets[i] + j, c);
					if (val != 0)
					{
						poly += mB->getLookup(mB->getMonomialfromBasis(poldeg, j), val);
						empty = false;
					}
				}
				if(!empty)
				{
					assert(poly.isHomogeneous());
					row.insert(pair<unsigned int, Poly<_Scalar> >(c,poly));
					cols[c].insert(i);
				}
			}

			rows.push_back(row);
		}
	}
	else
	{

		for (unsigned int i = 0; i < sx; i++)
		{
			unsigned int bSize = mB->getIndexOfBasisDegStart(degree + outdegs[i] + 1);
			;
			rowoffsets[i + 1] = rowoffsets[i] + bSize;
			spPRow<_Scalar> row;
			for (unsigned int c = 0; c < sy; c++)
			{
				Poly<_Scalar> poly(mB);
				bool empty = true;
				for (unsigned int j = 0; j < bSize; j++)
				{

					_Scalar val = inmatrix.randomreadaccess(rowoffsets[i] + j, c);
					if (val != 0)
					{
						poly += Poly<_Scalar>(mB, mB->getBasisElem(j), val);
						empty = false;
					}
				}
				if(!empty)
				{
					row.insert(pair<unsigned int, Poly<_Scalar> >(c,poly));\
					cols[c].insert(i);
				}
			}

			rows.push_back(row);
		}
	}
	assert(rowoffsets[sx] == inmatrix.getRows());
	computeDegrees(outdegs);
}

template <typename _Scalar>
const Poly<_Scalar> SparsePolyMatrix<_Scalar>::randomreadaccess(unsigned int r, unsigned int c) const
	{
		auto it = rows[r].find(c);
		if (it != rows[r].end())
		{
			return (*it).second;
		}
		else
		{
			return Poly<_Scalar>::zero(Ring);
		}

	}

template<typename _Scalar>
std::ostream &operator<<(std::ostream &os, SparsePolyMatrix<_Scalar> const &m)
{
	vector<vector<string> > out;
	unsigned int maxLabel = 4;
	for (unsigned int i = 0; i < m.outlabels.size(); i++)
	{
		unsigned int len = m.outlabels[i].size();
		if (maxLabel < len)
		{
			maxLabel = len;
		}
	}
	vector<unsigned int> maxlen(m.cols.size(), 0);
	for (unsigned int j = 0; j < m.cols.size(); j++)
	{
		if (m.inlabels.size() > j)
		{
			maxlen[j] = m.inlabels[j].size();
		}
		vector<string> col;
		for (unsigned int i = 0; i < m.rows.size(); i++)
		{
			ostringstream stream;
			stream << m.randomreadaccess(i,j);
			string s = stream.str();
			col.push_back(s);
			unsigned int len = s.size();
			if (maxlen[j] < len)
			{
				maxlen[j] = len;
			}
		}
		out.push_back(col);
	}

	os << std::setw(maxLabel) << " " << "| ";
	for (unsigned int j = 0; j < m.cols.size(); j++)
	{
		if (m.inlabels.size() > j)
		{
			os << std::setw(maxlen[j]) << m.inlabels[j] << " ";
		}
		else
		{
			os << std::setw(maxlen[j]) << m.indegrees[j] << " ";
		}
	}
	os << endl;
	os << std::setfill('-') << std::setw(maxLabel) << "-" << "| ";
	for (unsigned int j = 0; j < m.cols.size(); j++)
	{
		os << std::setw(maxlen[j]) << "" << " ";
	}
	os << std::setfill(' ') << endl;

	for (unsigned int i = 0; i < m.rows.size(); i++)
	{
		if (m.outlabels.size() > i)
		{
			os << std::setw(maxLabel) << m.outlabels[i] << "| ";
		}
		else
		{
			os << std::setw(maxLabel) << m.outdegrees[i] << "| ";
		}
		for (unsigned int j = 0; j < m.cols.size(); j++)
		{
			os << std::setw(maxlen[j]) << out[j][i] << " ";
		}
		os << endl;
	}
	return os;
}
template std::ostream &operator<< <mpz_class>(std::ostream &os, const SparsePolyMatrix<mpz_class>  &m);
template std::ostream &operator<< <mpq_class>(std::ostream &os, const SparsePolyMatrix<mpq_class>  &m);
template std::ostream &operator<< <numbermodulo>(std::ostream &os, const SparsePolyMatrix<numbermodulo>  &m);

//Print in such a way that Macaulay2 can read it in automatically
template<typename _Scalar>
void SparsePolyMatrix<_Scalar>::printForMacaulay2(std::ostream &os)
{
	os << "{";
	for (unsigned int i = 0; i < rows.size(); i++)
	{
		os << "{";
		for (unsigned int j = 0; j < cols.size(); j++)
		{
			randomreadaccess(i,j).print(os,true);
			if(j+1 < cols.size())
			{
				os << ", ";
			}
		}
		os << "}";
		if(i+1 < rows.size())
		{
			os << ", ";
		}
	}
	os << "}";
}

//Compute the kernel of the matrix in the twist given by degreeOfInterest via a Scalar representation
template<typename _Scalar>
SparsePolyMatrix<_Scalar> SparsePolyMatrix<_Scalar>::kernel(unsigned int degreeOfInterest, bool homogeneous, const program_options & opt)
{
	SparseScalarMatrix<_Scalar> inMat = this->degreeMatrix(degreeOfInterest, homogeneous);

	//cout << inMat << endl << endl;

	SparseScalarMatrix<_Scalar> kernel = inMat.kernel(opt);

	//cout << kernel << endl;

	SparsePolyMatrix<_Scalar> mBker(Ring, kernel, cols.size(), kernel.getCols(), indegrees, degreeOfInterest, homogeneous);
	mBker.setLabels(
	{ }, this->inlabels);
	return mBker;
}

//Construct the matrix the kernel of which is the symmetric power of the original kernel.
template<typename _Scalar>
SparsePolyMatrix<_Scalar> SparsePolyMatrix<_Scalar>::symMatrix(const unsigned int power) const
{
	vector<vector<unsigned int> > symn;
	vector<vector<unsigned int> > symm;
	kSumTon(symn, power, cols.size(), degrevlex);
	kSumTon(symm, power - 1, cols.size(), degrevlex);
	unsigned int ncols = symn.size();
	unsigned int nrows = symm.size() * rows.size();
	Poly<_Scalar> nullpoly(Ring);
	vector<spPRow<_Scalar> > rowlist;

	vector<int> outdegs, indegs;

	for (unsigned int u = 0; u < symm.size(); u++)
	{
		int d = dotproduct(symm[u], indegrees);
		for (unsigned int t = 0; t < rows.size(); t++)
		{
			outdegs.push_back(d + outdegrees[t]);
			spPRow<_Scalar> row;
			row.clear();
			for (unsigned int v = 0; v < symn.size(); v++)
			{
				int x[cols.size()];
				int pos = -1;
				for (unsigned int c = 0; c < cols.size(); c++)
				{
					x[c] = symn[v][c] - symm[u][c];
					if (x[c] < 0)
					{
						pos = -1;
						break;
					}
					if (x[c] == 1)
					{
						pos = c;
					}
				}
				if((pos != -1))
				{
					row.insert(pair<unsigned int, Poly<_Scalar> >(v,randomreadaccess(t,pos) * ((_Scalar) (symn[v][pos]))));
				}
			}

			assert(rowlist.size() < nrows);

			rowlist.push_back(row);
		}
	}

	for (unsigned int v = 0; v < symn.size(); v++)
	{
		int d = dotproduct(symn[v], indegrees);
		indegs.push_back(d);

	}

	return SparsePolyMatrix<_Scalar>(Ring, rowlist, ncols, outdegs, indegs);
}

//Construct the matrix the kernel of which is the exterior power of the original kernel.
template<typename _Scalar>
SparsePolyMatrix<_Scalar> SparsePolyMatrix<_Scalar>::extMatrix(const unsigned int power) const
{
	vector<vector<unsigned int> > subn;
	vector<vector<unsigned int> > subm;
	subsets(subn, cols.size(), power, degrevlex);
	subsets(subm, cols.size(), power - 1, degrevlex);
	unsigned int ncols = subn.size();
	Poly<_Scalar> nullpoly(Ring);
	vector<spPRow<_Scalar> > list;

	vector<int> outdegs;

	for (unsigned int u = 0; u < subm.size(); u++)
	{
		int d = partsum(subm[u], indegrees);
		for (unsigned int t = 0; t < rows.size(); t++)
		{
			outdegs.push_back(d+outdegrees[t]);
			spPRow<_Scalar> row;
			row.clear();
			for (unsigned int v = 0; v < subn.size(); v++)
			{
				unsigned int vpos = 0;
				unsigned int discrepancy = 0;
				for (unsigned int s = 0; s < subn[v].size(); s++)
				{
					if (s - discrepancy >= subm[u].size() || subn[v][s] != subm[u][s - discrepancy])
					{
						vpos = s;
						discrepancy++;
						if (discrepancy > 1)
							break;
					}
				}
				int sign = (vpos % 2 == 0) ? 1 : -1;
				if(discrepancy == 1)
				{
					row.insert(pair<unsigned int, Poly<_Scalar> >(v,randomreadaccess(t,subn[v][vpos])  * ((_Scalar) sign)));
				}
			}
			list.push_back(row);
		}
	}

	return SparsePolyMatrix<_Scalar>(Ring, list, ncols, outdegs);
}

//Construct the matrix given by taking each entry to the same power.
template<typename _Scalar>
SparsePolyMatrix<_Scalar> SparsePolyMatrix<_Scalar>::powerMatrix(const unsigned int power) const
{
	unsigned int ncols = cols.size();
	vector<spPRow<_Scalar> > list;

	vector<int> outdegs;

	for (unsigned int r = 0; r < rows.size(); r++)
	{
		spPRow<_Scalar> row;
		row.clear();
		for(auto & entry : rows[r])
		{
			row.insert(pair<unsigned int, Poly<_Scalar> >(entry.first,entry.second ^ power));
		}
		list.push_back(row);
		outdegs.push_back(outdegrees[r]*power);
	}

	return SparsePolyMatrix<_Scalar>(Ring,list, ncols, outdegs);
}

//Compute the scalar matrix describing this matrix in the given degree with respect to the monomial basis of the ring.
//The output is in the form of triplets (row, column, value). trows and tcols will give the total dimension of the matrix.
//If all entries are homogeneous use the homogeneous flag to save a lot of computational resources.
template<typename _Scalar>
template<typename tripScalar>
void SparsePolyMatrix<_Scalar>::degreeTriplets(const int degree, vector<rcvTriplet<tripScalar> > &tripletList, unsigned int& trows, unsigned int& tcols, bool homogeneous) const
{
	unsigned int rowoffsets[rows.size() + 1] =
	{ 0 };
	unsigned int coloffsets[cols.size() + 1] =
	{ 0 };
	int indeg[cols.size()];
	unsigned int cs[cols.size()];

	if (homogeneous)
	{

		for (unsigned int j = 0; j < cols.size(); j++)
		{
			indeg[j] = degree + indegrees[j];

			cs[j] = Ring->getBasisSize(indeg[j]);
			coloffsets[j + 1] = coloffsets[j] + cs[j];

		}

		for (unsigned int i = 0; i < rows.size(); i++)
		{
			int poldeg = degree + outdegrees[i];
			unsigned int rs = Ring->getBasisSize(poldeg);
			rowoffsets[i + 1] = rowoffsets[i] + rs;

			for (auto & entry : rows[i])
			{
				entry.second.appendMatrix(tripletList, poldeg, indeg[entry.first], rs, cs[entry.first], rowoffsets[i], coloffsets[entry.first]);
			}
		}
	}
	else
	{
		for (unsigned int j = 0; j < cols.size(); j++)
		{
			indeg[j] = degree + indegrees[j];

			cs[j] = Ring->getIndexOfBasisDegStart(degree - indegrees[j] + 1);
			coloffsets[j + 1] = coloffsets[j] + cs[j];

		}

		for (unsigned int i = 0; i < rows.size(); i++)
		{
			int poldeg = degree + outdegrees[i];
			unsigned int rs = Ring->getIndexOfBasisDegStart(poldeg + 1);
			rowoffsets[i + 1] = rowoffsets[i] + rs;


			for (auto & entry : rows[i])
			{
				entry.second.appendMatrix(tripletList, poldeg, indeg[entry.first], rs, cs[entry.first], rowoffsets[i], coloffsets[entry.first],homogeneous);
			}
		}
	}
	trows = rowoffsets[rows.size()];
	tcols = coloffsets[cols.size()];
}

//Compute the sparse scalar matrix describing this matrix in the given degree with respect to the monomial basis of the ring.
//If all entries are homogeneous use the homogeneous flag to save a lot of computational resources.
template<typename _Scalar>
SparseScalarMatrix<_Scalar> SparsePolyMatrix<_Scalar>::degreeMatrix(const int degree, bool homogeneous)
{
	unsigned int trows, tcols;
	vector<rcvTriplet<_Scalar> > tripletList;
	degreeTriplets<_Scalar>(degree, tripletList, trows, tcols, homogeneous);

	SparseScalarMatrix<_Scalar> M(tripletList, trows, tcols);
	return M;
}

//Matrix multiplication
template <typename Scalar>
SparsePolyMatrix<Scalar> SparsePolyMatrix<Scalar>::operator*(const SparsePolyMatrix<Scalar>& rhs) const
{
	const SparsePolyMatrix<Scalar> &lhs = *this;
	assert(lhs.cols.size() == rhs.rows.size());
	vector<spPRow<Scalar>> newrows;
	newrows.reserve(lhs.rows.size());
	for (unsigned int i = 0; i < lhs.rows.size(); i++)
	{
		spPRow<Scalar> row;
		for (auto &x : lhs.rows[i])
		{
			for (auto &y : rhs.rows[x.first])
			{
				Poly<Scalar> product = x.second * y.second;
				auto ret = row.insert(pair<unsigned int, Poly<Scalar> >(y.first, product));
				if (ret.second == false)
				{
					(*ret.first).second += product;
				}
			}
		}
		for (auto it = row.begin(); it != row.end();)
		{
			if ((*it).second.isNull())
			{
				auto f = (*it).first;
				row.erase(it);
				it = row.lower_bound(f);
			} else {
				it++;
			}
		}
		newrows.push_back(row);
	}
	return SparsePolyMatrix<Scalar>(lhs.Ring,newrows, rhs.cols.size(),lhs.outdegrees);
}


template class SparsePolyMatrix <mpz_class>;
template class SparsePolyMatrix <mpq_class>;
template class SparsePolyMatrix <numbermodulo> ;
