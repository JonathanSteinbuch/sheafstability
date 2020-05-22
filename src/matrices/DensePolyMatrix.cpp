//============================================================================
// Name        : DensePolyMatrix.cpp
// Author      : Jonathan Steinbuch
//============================================================================

#include "DensePolyMatrix.hpp"

#include "../helpers/Helpers.hpp"
#include "../helpers/Modulus.hpp"

template<typename _Scalar>
	int DensePolyMatrix<_Scalar>::indegree()
	{
		return sum(indegrees);
	}

	template<typename _Scalar>
	int DensePolyMatrix<_Scalar>::outdegree()
	{
		return sum(outdegrees);
	}

template<typename _Scalar>
bool DensePolyMatrix<_Scalar>::computeDegrees(const vector<int> &outdegs)
{
	vector<bool> inset;
	vector<bool> outset;

	outdegrees.clear();
	indegrees.clear();

	if (outdegs.size() == rows)
	{
		outdegrees = outdegs;
		outset.insert(outset.begin(), rows, true);
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
		outdegrees.insert(outdegrees.begin() + 1, rows - 1, 0);
		outset.push_back(true);
		outset.insert(outset.begin() + 1, rows - 1, false);
	}

	indegrees.insert(indegrees.begin(), cols, 0);
	inset.insert(inset.begin(), cols, false);

	for (unsigned int k = 0; k < rows; k++)
	{
		for (unsigned int i = 0; i < rows; i++)
		{
			for (unsigned int j = 0; j < cols; j++)
			{
				if (!entries[i][j].isNull())
				{
					if (!entries[i][j].isHomogeneous())
						return false;
					if (outset[i])
					{
						int ndeg = -entries[i][j].degree() + outdegrees[i];
						if (inset[j] && indegrees[j] != ndeg)
						{
							cout << *this << endl;
							return false;

						}
						else
						{
							indegrees[j] = ndeg;
							inset[j] = true;
						}
					}
					else
					{
						if (inset[j])
						{
							outdegrees[i] = indegrees[j] + entries[i][j].degree();
							outset[i] = true;
						}
					}
				}
			}
		}

		bool set = true;
		for (unsigned int j = 0; j < rows; j++)
		{
			if (!outset[j])
			{
				set = false;
				break;
			}
		}

		for (unsigned int j = 0; j < cols; j++)
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

	for (unsigned int j = 0; j < rows; j++)
	{
		if (!outset[j])
			return false;
	}

	for (unsigned int j = 0; j < cols; j++)
	{
		if (!inset[j])
			return false;
	}

	return true;
}

template<typename _Scalar>
DensePolyMatrix<_Scalar>::DensePolyMatrix(PolyRing<_Scalar>* mB,const vector<vector<Poly<_Scalar> > > &entries) :
		Ring(mB), rows(entries.size()), cols(entries[0].size()), entries(entries)
{
	computeDegrees(vector<int>());
}

template<typename _Scalar>
DensePolyMatrix<_Scalar>::DensePolyMatrix(PolyRing<_Scalar>* mB,const vector<vector<string > > &entries) :
		Ring(mB), rows(entries.size()), cols(entries[0].size())
{
	for (auto &strrow : entries)
	{
		vector<Poly<_Scalar> > polrow;
		for (auto &strentry : strrow)
		{
			Poly<_Scalar>  polentry = Poly<_Scalar> (mB, strentry);
			polrow.push_back(polentry);
		}
		this->entries.push_back(polrow);
	}

	computeDegrees(vector<int>());
}


template<typename _Scalar>
DensePolyMatrix<_Scalar>::DensePolyMatrix(PolyRing<_Scalar>* mB, const unsigned int sx, const unsigned int sy,const vector<vector<Poly<_Scalar> > > &entries) :
		Ring(mB), rows(sx), cols(sy), entries(entries)
{
	computeDegrees(vector<int>());
}

template<typename _Scalar>
DensePolyMatrix<_Scalar>::DensePolyMatrix(PolyRing<_Scalar>* mB, const unsigned int sx, const unsigned int sy,const vector<vector<Poly<_Scalar> > > &entries, const vector<int> &outdegs) :
		Ring(mB), rows(sx), cols(sy), entries(entries)
{
	computeDegrees(outdegs);
}

template<typename _Scalar>
DensePolyMatrix<_Scalar>::DensePolyMatrix(PolyRing<_Scalar>* mB, const unsigned int sx, const unsigned int sy,const vector<vector<Poly<_Scalar> > > &entries, const vector<int> &outdegs, const vector<int> &indegs) :
		Ring(mB), rows(sx), cols(sy), indegrees(indegs), outdegrees(outdegs), entries(entries)
{
}

template<typename _Scalar>
DensePolyMatrix<_Scalar>::DensePolyMatrix(PolyRing<_Scalar>* mB, const SparseScalarMatrix<_Scalar>& inmatrix, const unsigned int sx, const unsigned int sy, const vector<int> &outdegs, const int degree, bool homogeneous) :
		Ring(mB), rows(sx), cols(sy)
{
	unsigned int rowoffsets[rows + 1] =
	{ 0 };

	assert(outdegs.size() == rows);

	if (homogeneous)
	{

		for (unsigned int i = 0; i < rows; i++)
		{
			int poldeg = degree + outdegs[i];
			unsigned int bSize = mB->getBasisSize(poldeg);
			rowoffsets[i + 1] = rowoffsets[i] + bSize;
			vector<Poly<_Scalar> > row;
			for (unsigned int c = 0; c < cols; c++)
			{
				Poly<_Scalar> poly(mB);
				for (unsigned int j = 0; j < bSize; j++)
				{

					_Scalar val = inmatrix.randomreadaccess(rowoffsets[i] + j, c);
					if (val != 0)
					{
						poly += mB->getLookup(mB->getMonomialfromBasis(poldeg, j), val);
					}
				}
				assert(poly.isHomogeneous());
				row.push_back(poly);
			}

			entries.push_back(row);
		}
	}
	else
	{
		for (unsigned int i = 0; i < rows; i++)
		{
			unsigned int bSize = mB->getIndexOfBasisDegStart(degree + outdegs[i] + 1);

			rowoffsets[i + 1] = rowoffsets[i] + bSize;
			vector<Poly<_Scalar> > row;
			for (unsigned int c = 0; c < cols; c++)
			{
				Poly<_Scalar> poly(mB);
				for (unsigned int j = 0; j < bSize; j++)
				{

					_Scalar val = inmatrix.randomreadaccess(rowoffsets[i] + j, c);
					if (val != 0)
					{
						poly += Poly<_Scalar>(mB, mB->getBasisElem(j), val);
					}
				}
				//assert(poly.isHomogeneous());
				row.push_back(poly);
			}

			entries.push_back(row);
		}
	}
	assert(rowoffsets[rows] == inmatrix.getRows());
	computeDegrees(outdegs);
}
template<typename _Scalar>
std::ostream &operator<<(std::ostream &os, DensePolyMatrix<_Scalar> const &m)
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
	vector<unsigned int> maxlen(m.cols, 0);
	for (unsigned int j = 0; j < m.cols; j++)
	{
		if (m.inlabels.size() > j)
		{
			maxlen[j] = m.inlabels[j].size();
		}
		vector<string> col;
		for (unsigned int i = 0; i < m.rows; i++)
		{
			ostringstream stream;
			stream << m.entries[i][j];
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
	for (unsigned int j = 0; j < m.cols; j++)
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
	for (unsigned int j = 0; j < m.cols; j++)
	{
		os << std::setw(maxlen[j]) << "" << " ";
	}
	os << std::setfill(' ') << endl;

	for (unsigned int i = 0; i < m.rows; i++)
	{
		if (m.outlabels.size() > i)
		{
			os << std::setw(maxLabel) << m.outlabels[i] << "| ";
		}
		else
		{
			os << std::setw(maxLabel) << m.outdegrees[i] << "| ";
		}
		for (unsigned int j = 0; j < m.cols; j++)
		{
			os << std::setw(maxlen[j]) << out[j][i] << " ";
		}
		os << endl;
	}
	return os;
}

template<typename _Scalar>
void DensePolyMatrix<_Scalar>::printAsOperator(std::ostream &os)
{
	assert(outlabels.size() >= rows);
	for (unsigned int j = 0; j < cols; j++)
	{
		bool otf = false;
		for (unsigned int i = 0; i < rows; i++)
		{
			if (!entries[i][j].isNull())
			{
				if (otf)
				{
					if (entries[i][j].supportSize() > 1 || entries[i][j].leadingCoefficient() > 0)
					{
						os << " +";
					}
					else
					{
						os << " ";
					}
				}
				else
				{
					otf = true;
				}

				if (entries[i][j].supportSize() == 1)
				{
					os << entries[i][j] << outlabels[i];
				}
				else
				{
					os << "(" << entries[i][j] << ")" << outlabels[i];
				}
			}
		}
		os << endl;
	}
	os << endl;
}


template<typename _Scalar>
void DensePolyMatrix<_Scalar>::printForMacaulay2(std::ostream &os)
{
	os << "{";
	for (unsigned int i = 0; i < rows; i++)
	{
		os << "{";
		for (unsigned int j = 0; j < cols; j++)
		{
			entries[i][j].print(os,true);
			if(j+1 < cols)
			{
				os << ", ";
			}
		}
		os << "}";
		if(i+1 < rows)
		{
			os << ", ";
		}
	}
	os << "}";
}

template<typename _Scalar>
void DensePolyMatrix<_Scalar>::printForLatex(std::ostream &os)
{
	for (unsigned int i = 0; i < rows; i++)
	{
		for (unsigned int j = 0; j < cols; j++)
		{
			entries[i][j].print(os,true);
			if(j+1 < cols)
			{
				os << "& ";
			}
		}
		os << "\\\\" << endl;
	}
}


template<typename _Scalar>
DensePolyMatrix<_Scalar> directsum(const vector<DensePolyMatrix<_Scalar> > &summands)
{
	assert(summands.size() != 0);

	PolyRing<_Scalar>* mB = summands[0].Ring;

	vector<vector<Poly<_Scalar> > > list;
	unsigned int sx = 0;
	unsigned int sy = 0;
	for (unsigned int i = 0; i < summands.size(); i++)
	{
		sx += summands[i].getRows();
		sy += summands[i].getCols();
		assert(summands[i].Ring == mB);
	}
	list.reserve(sx);
	unsigned int ypos = 0;
	unsigned int xpos = 0;

	Poly<_Scalar> nullpoly(mB);
	for (unsigned int i = 0; i < summands.size(); i++)
	{
		unsigned int xend = xpos + summands[i].getRows();
		unsigned int yend = ypos + summands[i].getCols();
		for (unsigned int x = xpos; x < xend; x++)
		{
			vector<Poly<_Scalar> > row;
			row.reserve(sy);
			for (unsigned int y = 0; y < sy; y++)
			{
				if (y >= ypos && y < yend)
				{
					row.push_back(summands[i].getEntry(x - xpos, y - ypos));
				}
				else
				{
					row.push_back(nullpoly);
				}
			}
			list.push_back(row);
		}
		xpos = xend;
		ypos = yend;
	}
	return DensePolyMatrix<_Scalar>(mB, xpos, ypos, list);
}
template DensePolyMatrix<mpz_class> directsum<mpz_class>(const vector<DensePolyMatrix<mpz_class> > &summands);
template DensePolyMatrix<mpq_class> directsum<mpq_class>(const vector<DensePolyMatrix<mpq_class> > &summands);
template DensePolyMatrix<numbermodulo> directsum<numbermodulo>(const vector<DensePolyMatrix<numbermodulo> > &summands);

template<typename _Scalar>
DensePolyMatrix<_Scalar> directsum(const vector<const DensePolyMatrix<_Scalar>*> &summands)
{
	assert(summands.size() != 0);

	PolyRing<_Scalar>* mB = summands[0]->getRing();

	vector<vector<Poly<_Scalar> > > list;
	unsigned int sx = 0;
	unsigned int sy = 0;
	for (unsigned int i = 0; i < summands.size(); i++)
	{
		sx += summands[i]->getRows();
		sy += summands[i]->getCols();
		assert(summands[i]->getRing() == mB);
	}
	list.reserve(sx);
	unsigned int ypos = 0;
	unsigned int xpos = 0;

	Poly<_Scalar> nullpoly(mB);
	for (unsigned int i = 0; i < summands.size(); i++)
	{
		unsigned int xend = xpos + summands[i]->getRows();
		unsigned int yend = ypos + summands[i]->getCols();
		for (unsigned int x = xpos; x < xend; x++)
		{
			vector<Poly<_Scalar> > row;
			row.reserve(sy);
			for (unsigned int y = 0; y < sy; y++)
			{
				if (y >= ypos && y < yend)
				{
					row.push_back(summands[i]->getEntry(x - xpos, y - ypos));
				}
				else
				{
					row.push_back(nullpoly);
				}
			}
			list.push_back(row);
		}
		xpos = xend;
		ypos = yend;
	}
	return DensePolyMatrix<_Scalar>(mB, xpos, ypos, list);
}
template DensePolyMatrix<mpz_class> directsum<mpz_class>(const vector<const DensePolyMatrix<mpz_class>* > &summands);
template DensePolyMatrix<mpq_class> directsum<mpq_class>(const vector<const DensePolyMatrix<mpq_class>* > &summands);
template DensePolyMatrix<numbermodulo> directsum<numbermodulo>(const vector<const DensePolyMatrix<numbermodulo>* > &summands);


template<typename _Scalar>
DensePolyMatrix<_Scalar> sum(const vector<const DensePolyMatrix<_Scalar>*> &summands)
{
	assert(summands.size() != 0);

	PolyRing<_Scalar>* mB = summands[0]->Ring;

	vector<vector<Poly<_Scalar> > > list;
	unsigned int sx = summands[0]->getRows();
	unsigned int sy = 0;
	for (unsigned int i = 0; i < summands.size(); i++)
	{
		sy += summands[i]->getCols();
		assert(summands[i]->Ring == mB);
		assert(summands[i]->getRows() == sx);
	}
	list.reserve(sx);
	unsigned int ypos = 0;

	Poly<_Scalar> nullpoly(mB);
	for (unsigned int x = 0; x < sx; x++)
	{
		vector<Poly<_Scalar> > row;
		row.reserve(sy);
		for (unsigned int i = 0; i < summands.size(); i++)
		{
			unsigned int yend = summands[i]->getCols();
			for (unsigned int y = 0; y < yend; y++)
			{
				row.push_back(summands[i]->getEntry(x, y));
			}

			ypos += yend;
		}
		list.push_back(row);
	}
	return DensePolyMatrix<_Scalar>(mB, sx, sy, list);
}
template DensePolyMatrix<mpz_class> sum<mpz_class>(const vector<const DensePolyMatrix<mpz_class>* > &summands);
template DensePolyMatrix<mpq_class> sum<mpq_class>(const vector<const DensePolyMatrix<mpq_class>* > &summands);
template DensePolyMatrix<numbermodulo> sum<numbermodulo>(const vector<const DensePolyMatrix<numbermodulo>* > &summands);

template<typename _Scalar>
DensePolyMatrix<_Scalar> DensePolyMatrix<_Scalar>::ntimes(const unsigned int power) const
{
	vector<const DensePolyMatrix<_Scalar>*> summands;
	summands.reserve(power);
	for (unsigned int u = 0; u < power; u++)
	{
		summands.push_back(this);
	}
	return directsum(summands);
}

template<typename _Scalar>
DensePolyMatrix<_Scalar> DensePolyMatrix<_Scalar>::irrelevantPower(PolyRing<_Scalar> &mB, int degree)
{
	vector<Poly<_Scalar> > row;
	vector<vector<Poly<_Scalar> > > list;
	if (degree > 0)
	{
		for (unsigned int i = mB.getIndexOfBasisDegStart(degree); i < mB.getIndexOfBasisDegStart(degree + 1); i++)
		{
			row.push_back(Poly<_Scalar>(&mB, mB.getBasisElem(i), 1));
		}
		list.push_back(row);
	}
	return DensePolyMatrix<_Scalar>(&mB, 1, row.size(), list);
}

template<typename _Scalar>
DensePolyMatrix<_Scalar> DensePolyMatrix<_Scalar>::kernel(unsigned int degreeOfInterest, bool homogeneous, const  program_options &opt)
{
	SparseScalarMatrix<_Scalar> inMat = this->degreeMatrix(degreeOfInterest, homogeneous);

	//cout << inMat << endl << endl;

	SparseScalarMatrix<_Scalar> kernel = inMat.kernel(opt);

	//cout << kernel << endl;

	DensePolyMatrix<_Scalar> mBker(Ring, kernel, cols, kernel.getCols(), indegrees, degreeOfInterest, homogeneous);
	mBker.setLabels(
	{ }, this->inlabels);
	return mBker;
}

template<typename _Scalar>
DensePolyMatrix<_Scalar> DensePolyMatrix<_Scalar>::symMatrix(const unsigned int power) const
{
	vector<vector<unsigned int> > symn;
	vector<vector<unsigned int> > symm;
	kSumTon(symn, power, cols, degrevlex);
	kSumTon(symm, power - 1, cols, degrevlex);
	unsigned int ncols = symn.size();
	unsigned int nrows = symm.size() * rows;
	Poly<_Scalar> nullpoly(Ring);
	vector<vector<Poly<_Scalar> > > list;

	vector<int> outdegs, indegs;

	for (unsigned int u = 0; u < symm.size(); u++)
	{
		int d = dotproduct(symm[u], indegrees);
		for (unsigned int t = 0; t < rows; t++)
		{
			outdegs.push_back(d + outdegrees[t]);
			vector<Poly<_Scalar> > row;
			row.clear();
			for (unsigned int v = 0; v < symn.size(); v++)
			{
				int x[cols];
				int pos = -1;
				for (unsigned int c = 0; c < cols; c++)
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
				row.push_back((pos != -1) ? entries[t][pos] * ((_Scalar) (symn[v][pos])) : nullpoly);
			/*	if((pos != -1))
				{
					notnull ++;
				}
				all++;*/
			}

			//cout << notnull << " / " << all << " *>= " << sizeof(nullpoly) << endl;
			assert(list.size() < nrows);

			list.push_back(row);
		}
	}

	for (unsigned int v = 0; v < symn.size(); v++)
	{
		int d = dotproduct(symn[v], indegrees);
		indegs.push_back(d);

	}

	return DensePolyMatrix<_Scalar>(Ring, nrows, ncols, list, outdegs, indegs);
}

template<typename _Scalar>
DensePolyMatrix<_Scalar> DensePolyMatrix<_Scalar>::extMatrix(const unsigned int power) const
{
	vector<vector<unsigned int> > subn;
	vector<vector<unsigned int> > subm;
	subsets(subn, cols, power, degrevlex);
	subsets(subm, cols, power - 1, degrevlex);
	unsigned int ncols = subn.size();
	unsigned int nrows = subm.size() * rows;
	Poly<_Scalar> nullpoly(Ring);
	vector<vector<Poly<_Scalar> > > list;

	vector<int> outdegs;

	for (unsigned int u = 0; u < subm.size(); u++)
	{
		int d = partsum(subm[u], indegrees);
		for (unsigned int t = 0; t < rows; t++)
		{
			outdegs.push_back(d+outdegrees[t]);
			vector<Poly<_Scalar> > row;
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
				row.push_back((discrepancy == 1) ? entries[t][subn[v][vpos]] * ((_Scalar) sign) : nullpoly);
			}
			list.push_back(row);
		}
	}

	return DensePolyMatrix<_Scalar>(Ring, nrows, ncols, list, outdegs);
}

template<typename _Scalar>
DensePolyMatrix<_Scalar> DensePolyMatrix<_Scalar>::powerMatrix(const unsigned int power) const
{
	unsigned int ncols = cols;
	unsigned int nrows = rows;
	vector<vector<Poly<_Scalar> > > list;

	vector<int> outdegs;

	for (unsigned int r = 0; r < rows; r++)
	{
		vector<Poly<_Scalar> > row;
		row.clear();
		for (unsigned int c = 0; c < cols; c++)
		{
			row.push_back(entries[r][c] ^ power);
		}
		list.push_back(row);
		outdegs.push_back(outdegrees[r]*power);
	}

	return DensePolyMatrix<_Scalar>(Ring, nrows, ncols, list, outdegs);
}

template<typename _Scalar>
template<typename tripScalar>
vector<rcvTriplet<tripScalar> > DensePolyMatrix<_Scalar>::degreeTriplets(const int degree, unsigned int& trows, unsigned int& tcols, bool homogeneous)
{
	unsigned int rowoffsets[rows + 1] =
	{ 0 };
	unsigned int coloffsets[cols + 1] =
	{ 0 };
	int indeg[cols];
	unsigned int cs[cols];

	vector<rcvTriplet<tripScalar> > tripletList;

	if (homogeneous)
	{

		for (unsigned int j = 0; j < cols; j++)
		{
			indeg[j] = degree + indegrees[j];

			cs[j] = Ring->getBasisSize(indeg[j]);
			coloffsets[j + 1] = coloffsets[j] + cs[j];

		}

		for (unsigned int i = 0; i < rows; i++)
		{
			int poldeg = degree + outdegrees[i];
			unsigned int rs = Ring->getBasisSize(poldeg);
			rowoffsets[i + 1] = rowoffsets[i] + rs;

			for (unsigned int j = 0; j < cols; j++)
			{

				if (!entries[i][j].isNull())
				{
					entries[i][j].appendMatrix(tripletList, poldeg, indeg[j], rs, cs[j], rowoffsets[i], coloffsets[j]);
				}
			}
		}
	}
	else
	{
		for (unsigned int j = 0; j < cols; j++)
		{
			indeg[j] = degree + indegrees[j];

			cs[j] = Ring->getIndexOfBasisDegStart(degree - indegrees[j] + 1);
			coloffsets[j + 1] = coloffsets[j] + cs[j];

		}

		for (unsigned int i = 0; i < rows; i++)
		{
			int poldeg = degree + outdegrees[i];
			unsigned int rs = Ring->getIndexOfBasisDegStart(poldeg + 1);
			rowoffsets[i + 1] = rowoffsets[i] + rs;

			for (unsigned int j = 0; j < cols; j++)
			{

				if (!entries[i][j].isNull())
				{
					entries[i][j].appendMatrix(tripletList, poldeg, indeg[j], rs, cs[j], rowoffsets[i], coloffsets[j], homogeneous);
				}
			}
		}
	}
	trows = rowoffsets[rows];
	tcols = coloffsets[cols];
	return tripletList;
}

/*template<typename _Scalar>
libnormaliz::Matrix<_Scalar> mBmatrix<_Scalar>::nmzdegreeMatrix(const int degree)
{
	unsigned int trows, tcols;
	vector<Eigen::Triplet<_Scalar> > tripletList = degreeTriplets<_Scalar>(degree, trows, tcols);
	cout << "libnormaliz Matrix from Triplets..." << endl;
	vector<vector<_Scalar> > nmzentries;
	nmzentries.reserve(trows);
	vector<_Scalar> row;
	row.reserve(tcols);
	for (unsigned int j = 0; j < tcols; j++)
	{
		row.push_back(0);
	}
	for (unsigned int i = 0; i < trows; i++)
	{
		nmzentries.push_back(row);
	}
	for (unsigned int c = 0; c < tripletList.size(); c++)
	{
		nmzentries[tripletList[c].row()][tripletList[c].col()] = tripletList[c].value();
	}

	return libnormaliz::Matrix<_Scalar>(nmzentries);
}*/

template<typename _Scalar>
SparseScalarMatrix<_Scalar> DensePolyMatrix<_Scalar>::degreeMatrix(const int degree, bool homogeneous)
{
	unsigned int trows, tcols;
	vector<rcvTriplet<_Scalar> > tripletList = degreeTriplets<_Scalar>(degree, trows, tcols, homogeneous);

	SparseScalarMatrix<_Scalar> M(tripletList, trows, tcols);
	return M;
}

template<typename _Scalar>
DensePolyMatrix<_Scalar> DensePolyMatrix<_Scalar>::operator*(const DensePolyMatrix<_Scalar>& rhs) const
{
	const DensePolyMatrix<_Scalar> &lhs = *this;
	vector<vector<Poly<_Scalar> > > newentries;
	assert(rhs.rows == lhs.cols);
	for (unsigned int i = 0; i < lhs.rows; i++)
	{
		vector<Poly<_Scalar> > row;
		for (unsigned int j = 0; j < rhs.cols; j++)
		{
			Poly<_Scalar> poly(lhs.Ring);
			for (unsigned int k = 0; k < rhs.rows; k++)
			{
				poly += lhs.entries[i][k] * rhs.entries[k][j];
			}
			row.push_back(poly.trim());
		}
		newentries.push_back(row);
	}
	DensePolyMatrix<_Scalar> ret(lhs.Ring, lhs.rows, rhs.cols, newentries, lhs.outdegrees);
	return ret;
}

template<typename _Scalar>
DensePolyMatrix<_Scalar> DensePolyMatrix<_Scalar>::operator%(const DensePolyMatrix<_Scalar>& rhs) const
{
	const DensePolyMatrix<_Scalar> &lhs = *this;
	if (lhs.cols == 0)
		return lhs;
	vector<const DensePolyMatrix<_Scalar>*> list;
	list.push_back(&lhs);
	list.push_back(&rhs);
	DensePolyMatrix<_Scalar> sumx = sum(list);
	sumx.computeDegrees(lhs.getOutdegrees());

	//cout << sumx << endl;

	int d = max(lhs.indegrees);
	SparseScalarMatrix<_Scalar> inMat = sumx.degreeMatrix(d);

	/*	auto inMat2 = inMat;
	 dbspMatrix<_Scalar> k = inMat2.kernel();

	 mBmatrix<_Scalar> t = mBmatrix(lhs.mB,k,sumx.cols,k.getCols(),sumx.indegrees,d);

	 cout << t << endl;*/

	//	cout << inMat << endl;
	SparseScalarMatrix<_Scalar> reducedMat = inMat.modulo(lhs.getCols());

	auto m = DensePolyMatrix<_Scalar>(lhs.Ring, reducedMat, sumx.rows, reducedMat.getCols(), sumx.outdegrees, d);
	m.setLabels(
	{ }, lhs.outlabels);
	return m;
}

template<typename _Scalar>
DensePolyMatrix<_Scalar> DensePolyMatrix<_Scalar>::jacobiMatrix()
{
	assert(rows == 1);

	vector<int> jOutdegrees;
	vector<string> jOutlabels;
	vector<string> jInlabels;

	unsigned int snu = Ring->getNumVars() + 1;

	for (unsigned int nu = 1; nu < snu; nu++)
	{
		ostringstream stream;
		stream << "d";
		Ring->print(stream, nu);
		string s = stream.str();
		jInlabels.push_back(s);
	}

	vector<vector<Poly<_Scalar> > > jEntries;
	for (unsigned int c = 0; c < cols; c++)
	{
		ostringstream stream;
		stream << entries[0][c];
		string s = stream.str();
		jOutlabels.push_back(s);

		vector<Poly<_Scalar> > row;
		for (unsigned int nu = 1; nu < snu; nu++)
		{
			row.push_back(entries[0][c].taylorDerive(nu));
		}
		jEntries.push_back(row);
	}

	DensePolyMatrix<_Scalar> M(Ring, cols, snu - 1, jEntries, jOutdegrees);
	M.setLabels(jInlabels, jOutlabels);

	return M;
}

template<typename _Scalar>
DensePolyMatrix<_Scalar> DensePolyMatrix<_Scalar>::jacobiTaylor(PolyRing<_Scalar>* mB, const Poly<_Scalar>& poly, unsigned int n)
{
	PolyRing<mpz_class> rB(mB->getNumRealVars(), mB->getPrintVar(), degrevlex, mB->getNumRealVars(),{ }, true);
	unsigned int snu = rB.getMonN(n + 1);
	unsigned int smu = rB.getMonN(n);
	vector<Poly<_Scalar> > parderivs;
	vector<int> outdegrees;
	vector<string> outlabels;
	vector<string> inlabels;
	for (unsigned int nu = 0; nu < snu; nu++)
	{
		vector<unsigned int> v = rB.getExponents(nu);
		parderivs.push_back(poly.taylorDerive(mB->getMonomial(v)));

		if (nu >= 1)
		{
			ostringstream stream;
			stream << "d";
			rB.print(stream, nu);
			string s = stream.str();
			inlabels.push_back(s);
		}

	}
	for (unsigned int mu = 0; mu < smu; mu++)
	{
		vector<unsigned int> v = rB.getExponents(mu);
		outdegrees.push_back(mB->getDegree(mB->getMonomial(v)));

		ostringstream stream;
		stream << "d";
		rB.print(stream, mu);
		string s = stream.str();
		outlabels.push_back(s);
	}
	Poly<_Scalar> zero = Poly<_Scalar>::zero(mB);
	vector<vector<Poly<_Scalar> > > entries;

	for (unsigned int mu = 0; mu < smu; mu++)
	{
		vector<Poly<_Scalar> > row;
		for (unsigned int nu = 1; nu < snu; nu++)
		{
			mBmonomial diff = rB.divide(nu, mu);
			if (diff.isError())
			{
				row.push_back(zero);
			}
			else
			{
				row.push_back(parderivs[diff.getId()]);
			}
		}
		entries.push_back(row);
	}

	DensePolyMatrix<_Scalar> M(mB, smu, snu-1, entries, outdegrees);
	M.setLabels(inlabels, outlabels);

	return M;
}

template<typename _Scalar>
const vector<Poly<_Scalar> > DensePolyMatrix<_Scalar>::computeMinors(const unsigned int n) const
{
	vector<vector<unsigned int> > colsubs;
	subsets(colsubs, cols, n, degrevlex);

	vector<vector<unsigned int> > rowsubs;
	subsets(rowsubs, rows, n, degrevlex);

	vector<vector<unsigned int> > perms;
	vector<int> parlist;
	permutationsandparities(perms, n, parlist);

	vector<Poly<_Scalar> > minors;

	for (auto & rsub : rowsubs)
	{
		for (const auto & csub : colsubs)
		{
			Poly<_Scalar> minorpol(Ring);

			for (unsigned int p = 0; p < perms.size(); p++)
			{
				Poly<_Scalar> summand(Ring, (_Scalar) parlist[p]);
				for (unsigned int c = 0; c < n; c++)
				{
					summand *= entries[rsub[perms[p][c]]][csub[c]];
					if (summand.isNull())
					{
						break;
					}
				}
				minorpol += summand;
			}
			minors.push_back(minorpol);
		}
	}
	return minors;
}

template class DensePolyMatrix <mpz_class>;
template class DensePolyMatrix <mpq_class>;
template class DensePolyMatrix <numbermodulo>;
