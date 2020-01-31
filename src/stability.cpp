//============================================================================
// Name        : stability.cpp
// Author      : Jonathan Steinbuch
//
//============================================================================

#include <gmpxx.h>
#include <iostream>
#include <boost/array.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>

#include "Modulus.hpp"

#include "PolyRing.hpp"
#include "DensePolyMatrix.hpp"
#include "SparsePolyMatrix.hpp"
#include "Poly.hpp"
#include "Helpers.hpp"

using namespace std;

const float version1 = 0;
const float version2 = 5;
const float version3 = 0;

datatype numbermodulo::modulus;
invTable numbermodulo::lookup;

typedef Poly<mpz_class> Zpoly;
typedef PolyRing<mpz_class> ZRing;
typedef DensePolyMatrix<mpz_class> ZMatrix;
typedef SparsePolyMatrix<mpz_class> spZMatrix;

typedef Poly<numbermodulo> Fpoly;
typedef PolyRing<numbermodulo> FRing;
typedef DensePolyMatrix<numbermodulo> FMatrix;
typedef SparsePolyMatrix<numbermodulo> spFMatrix;

//kernel computation for the semistability check
template <typename Scalar>
int intest(SparsePolyMatrix<Scalar> Matrix, int degreeOfInterest, SparsePolyMatrix<Scalar> &ker, program_options opt)
{
	SparseScalarMatrix<Scalar> inMat(Matrix, degreeOfInterest);

	if (opt.verbosity >= 1)
	{

		cout << "Computing Kernel of " << Matrix.getRows() << "x" << Matrix.getCols() << "-Matrix (internally " << inMat.getRows() << "x" << inMat.getCols() << ")...";
		cout << endl;
	}

	if (opt.verbosity >= 3)
	{
		cout << "Internal Matrix:" << endl;
		cout << inMat;
		cout << endl << endl;
	}

	SparseScalarMatrix<Scalar> kernel = inMat.kernel(opt);

	if (kernel.getCols() > 0)
	{
		SparsePolyMatrix<Scalar> mBker(Matrix.getRing(), kernel, Matrix.getCols(), kernel.getCols(), Matrix.getIndegrees(), degreeOfInterest);

		ker = mBker;
	}

	return kernel.getCols();

}

//Main function implementing the check for semistability
template <typename Scalar>
bool checkSemistability(PolyRing<Scalar> &Ring, DensePolyMatrix<Scalar> &M, unsigned int characteristic, ofstream &output, const program_options opt = noOptions)
{
	bool result = false;

	Ring.setOptions(opt); //make sure the ring has the correct program options

	assert(M.getRows() == 1); //For now, we can only decide semistability of syzygy sheafs

	if (opt.verbosity >= 1) //Write out some information about what we are going to do
	{
		cout << "We compute if the Syzygy sheaf given by M is semistable over the projective variety X, where:" << endl;
		string field;
		if (characteristic == 0)
		{
			cout << "X is given by \\Q[";
		} else
		{
			cout << "X is given by \\F_" << characteristic << "[";
		}

		for (unsigned int i = 0; i < Ring.getNumVars(); i++)
		{
			cout << Poly<Scalar>::variable(&Ring, i);
			if (i != Ring.getNumVars() - 1)
			{
				cout << ",";
			}
		}
		cout << "]/(";
		for (unsigned int i = 0; i < Ring.getGB().size(); i++)
		{
			cout << Ring.getGB()[i];
			if (i != Ring.getGB().size() - 1)
			{
				cout << ",";
			}
		}
		cout << ")" << endl;
		cout << "And M= " << endl << M << endl;
	}

	//Check if we have a curve.
	PolyRing<mpq_class> kX(1,{ "n" }, degrevlex); //Polynomial Ring in One variable for the Hilbert Polynomial
	kX.setOptions(opt);
	Poly<mpq_class> hP(&kX); //Container for the Hilbert Polynomial
	Ring.getBaseRing()->quotientHilbertPoly(hP, Ring.getGB());

	if (opt.verbosity >= 1)
	{
		cout << "The Hilbert polynomial of X is P(n)=" << hP << endl;
	}

	if (hP.degree() != 1)
	{
		cerr << "Not a Curve! Computation aborted." << endl;
		return false;
	}

	//Check if the curve is smooth
	if (!Ring.isSmooth())
	{
		cerr << "Curve not smooth! Computation aborted." << endl;
		return false;
	}
	else if (opt.verbosity >= 1)
	{
		cout << "Curve is smooth: Check." << endl;
	}

	//Compute degree and genus of the curve
	unsigned int reldeg = hP.leadingCoefficient().get_num().get_ui(); //The degree is the leading coefficient of the Hilbert Polynomial
	assert(hP.leadingCoefficient().get_den() == 1); //If the leading coefficient is not an integer, something went wrong
	int genus = 1;
	if (hP.supportSize() > 1) //This check is needed because if hP has only one nonzero monomial the trailing monomial wouldn't capture the constant coefficient
	{
		genus = (1 - hP.trailingCoefficient().get_num().get_si()); //Genus g = 1 - hP(0).
	}

	if (opt.verbosity >= 1)
	{
		cout << "Genus of Curve: " << genus << endl << endl;
	}

	//Compute the image of M and its degree
	auto Mgens = M.getRowEntries(0); //get entries of M as a family of elements of the Ring
	for (auto & p : Mgens)
	{
		p.changeRing(Ring.getBaseRing()); //Replace the ring elements by the polynomials representing them
	}
	Mgens.insert(Mgens.end(), Ring.getGB().begin(), Ring.getGB().end()); //Add the elements of the ideal to the entries
	Ring.getBaseRing()->computeGroebnerBasis(Mgens); //Compute the Gröbner Basis of the resulting family
	Poly<mpq_class> MhP(&kX); //Container for the Hilbert Polynomial of BaseRing/Mgens
	Ring.getBaseRing()->quotientHilbertPoly(MhP, Mgens);
	if (opt.verbosity >= 1)
	{
		cout << "The Hilbert polynomial of M is P(n)=" << MhP << endl;
		cout << endl;
	}
	assert(MhP.degree() == 0); //The Hilbert Polynomial should be a constant
	int idealdegree = -MhP.leadingCoefficient().get_num().get_si(); //The degree of the image is d=-MhP.

	//Compute rank degree and slope of the kernel of M
	int rank = M.getCols() - M.getRows();
	int degree = M.indegree()*reldeg - idealdegree;
	double slope = (double) degree / (double) rank;

	unsigned int threshold = genus - 1 + reldeg;

	SparsePolyMatrix<Scalar> sparseM = SparsePolyMatrix<Scalar>(M);


	//Decide with which variant of the theorem we work
	int maxExt;
	if(opt.exteriorPowers)
	{
		maxExt = rank - 1;
	} else {
		maxExt = 1;
	}

	for (int extPower = 1; extPower <= maxExt; extPower++)
	{
		auto extMatrix = sparseM.extMatrix(extPower);
		double extSlope = slope * extPower;

		if (opt.verbosity >= 1)
		{
			cout << "Exterior Power: " << extPower << ", ";
			cout << extMatrix.getRows() << "x" << extMatrix.getCols() << "-Matrix" << endl;
		}

		int gc; //temporary variable for the gcd
		int n; //The n from the theorem
		if(opt.exteriorPowers)
		{
			gc = __gcd(rank, -extPower*degree);
			n = rank/gc;
		} else {
			gc = __gcd(rank, -degree);
			n = (rank - 1)*rank/gc;
		}

		unsigned int actualThreshold = threshold * n + 1; //The power q from the Theorem

		if (opt.verbosity >= 1)
		{
			cout << "Computing until power " << actualThreshold << endl;
		}

		//window of powers in which to look for destabilizing sheaves:
		unsigned int currmax = actualThreshold;
		unsigned int currmin = 0;

		unsigned int sPower = 1;
		unsigned int oldsPower = 1;

		SparsePolyMatrix<Scalar> sMatrix;

		//Try symmetric powers of M until we get to q, to maybe get lucky and find a destabilizing sheaf earlier
		while (currmin < currmax)
		{
			if (characteristic == 0)
			{
				if (opt.verbosity >= 1)
				{
					cout << "Symmetric Power: " << sPower << endl;
				}
				sMatrix = extMatrix.symMatrix(sPower); //Compute the Matrix to the Symmetric Power
			}
			else
			{
				if (opt.verbosity >= 1)
				{
					cout << "Frobenius Power: " << sPower << endl;
				}
				//Compute the Matrix to the Frobenius Power
				if(sPower > oldsPower)
				{
					sMatrix = sMatrix.powerMatrix(sPower/oldsPower); //Compute Incrementally. It's less work than to compute it from the original matrix every time
				} else {
					sMatrix = extMatrix.powerMatrix(sPower);
				}

				oldsPower = sPower;
			}

			if (opt.verbosity >= 2)
			{
				cout << "Power matrix:" << endl;
				cout << sMatrix << endl;
			}

			if (opt.verbosity >= 1)
			{
				cout << sMatrix.getRows() << "x" << sMatrix.getCols() << "-Matrix" << endl;
				cout << endl;
			}

			double sSlope = extSlope * sPower;
			int twist = ((int) ceil(-sSlope/reldeg) - 1); //The k from the theorem

			if (opt.verbosity >= 1)
			{
				cout << "Slope: " << sSlope << endl;
				cout << "Twist: " << twist << endl;
			}

			clock_t t = clock();

			SparsePolyMatrix<Scalar> ker; //Container for the kernel of the power Matrix

			int kernelRank = intest(sMatrix, twist, ker, opt); //Compute the rank of the kernel of the power matrix

			if (opt.verbosity >= 1)
			{
				t = clock() - t;
				cout << "Clock: " << t << "ticks or " << t / CLOCKS_PER_SEC << "s" << endl;
				cout << "Rank at degree of Interest: " << kernelRank << endl;
			}
			if (kernelRank > 0) //We found a destabilizing section
			{
				if (opt.verbosity >= 1)
				{
					cout << "Not Semistable!" << endl;
				}
				//Write the power and twist where we found the destabilizing section into the output file
				output << 1 << endl;
				output << "{" << "ExteriorPower => " << extPower << ", SymmetricPower => " << sPower << ", Twist => " << twist << "}" << endl;
				if(opt.verbosity >= 2)
				{
					cout << "Kernel:" << endl;
					cout << ker << endl << endl;
				}
				if(opt.verbosity >= 3) //If the verbosity is high even write out the section itself
				{
					if(opt.outputM2)
					{
						ker.printForMacaulay2(output);
					} else {
						output << ker;
					}
				}
				output << endl;
				result = false;
				if(characteristic == 0)
				{
					currmax = sPower - 1; //Adjust the power for which we may maybe not find destabilizing sections
				} else {
					currmax = sPower / characteristic;
				}
				if(opt.stopIfUnstable)
				{
					return result;
				}
			}
			else if (extPower == maxExt && sPower == actualThreshold && kernelRank == 0)
			{
				if (opt.verbosity >= 1)
				{
					cout << "Semistable!" << endl;
				}
				output << 0 << endl;
				result = true;
				return result;
			}
			else
			{
				currmin = sPower;
			}
			if (opt.verbosity >= 1)
			{
				cout << endl;
			}

			if (characteristic == 0)
			{
				if(opt.linearPowerProgresison)
				{
					sPower++;
				} else{
					//some heuristical power progression
					if (currmax == actualThreshold)
					{
						if (sPower * 2 <= actualThreshold)
							sPower *= 2;
						else
							sPower = actualThreshold;
					}
					else
					{
						sPower = (currmax - currmin) / 2 + currmin;
						if (sPower == currmin)
							sPower++;
					}
				}
			}
			else
			{
				sPower *= characteristic;
				if (sPower == 1)
				{
					sPower = 0;
				}
			}

		}
	}
	return result;
}

//Check if we think that the kernels can be computed over \Z instead of \Q
bool ZViable(ZRing &mB,  const vector<string> &relations, vector<Zpoly > &polrels){
	for (auto &strrelation : relations)
	{
		Zpoly polrel = Poly<mpz_class>(&mB, strrelation);
		polrels.push_back(polrel);
	}
	mB.computeGroebnerBasis(polrels);
	for(auto & gen : polrels)
	{
		if(gen.leadingCoefficient() != 1 && gen.leadingCoefficient() != -1)
		{
			return false;
			break;
		}
	}
	return true;
}

//Choose over which ring to perform the computations
bool semistabilityChoice(vector<string> &variables, vector<string> &relations,vector<vector<string>> &matrix, unsigned int characteristic, ofstream &output, const program_options & opt = noOptions)
{
	if(characteristic == 0)
	{
		if(opt.forceQ)
		{
			PolyRing<mpq_class> mB(variables.size(), variables, degrevlex);
			PolyRing<mpq_class> QuotientRing(&mB, relations);
			DensePolyMatrix<mpq_class> M(&QuotientRing, matrix);
			return checkSemistability(*(M.getRing()), M, characteristic, output, opt);
		} else {
			ZRing mB(variables.size(), variables, degrevlex);

			vector<Zpoly > polrels;
			if(ZViable(mB,relations,polrels))
			{
				ZRing QuotientRing(&mB, polrels);

				ZMatrix M(&QuotientRing, matrix);
				return checkSemistability(QuotientRing, M, characteristic, output, opt);
			} else {
				PolyRing<mpq_class> mBq(variables.size(), variables, degrevlex);
				PolyRing<mpq_class> QuotientRing(&mBq, relations);
				DensePolyMatrix<mpq_class> M(&QuotientRing, matrix);
				return checkSemistability(*(M.getRing()), M, characteristic, output, opt);
			}
		}

	} else {
		numbermodulo::setModulus(characteristic);

		PolyRing<numbermodulo> mB(variables.size(), variables, degrevlex);
		PolyRing<numbermodulo> QuotientRing(&mB, relations);
		DensePolyMatrix<numbermodulo> M(&QuotientRing, matrix);
		return checkSemistability(*(M.getRing()), M, characteristic, output, opt);
	}
}

template <typename Scalar>
void powerMatrix( DensePolyMatrix<Scalar> &M, unsigned int characteristic,  unsigned int symPower, unsigned int extPower,ofstream &output,  const program_options &opt = noOptions){
	M.getRing()->setOptions(opt);

	if (opt.verbosity >= 1 )
	{
		if(characteristic == 0)
		{
			cout << "Computing Sym^" << symPower << "(Ext^" << extPower << "(M))." << endl;
		} else {
			cout << "Computing F*^" << symPower << "(Ext^" << extPower << "(M))." << endl;
		}
	}
	if (opt.verbosity >= 2)
	{
		cout << "Input Matrix:" << endl;
		cout << M << endl;
	}
	SparsePolyMatrix<Scalar> Msparse = M;
	SparsePolyMatrix<Scalar> Mext = M.extMatrix(extPower);
	SparsePolyMatrix<Scalar> Msym;
	if(characteristic == 0)
	{
		Msym = Mext.symMatrix(symPower);
	} else {
		Msym = Mext.powerMatrix(symPower);
	}

	if (opt.verbosity >= 2)
	{
		cout << "Power Matrix:" << endl;
		cout << Msym << endl;
	}
	if(opt.outputM2)
	{
		Msym.printForMacaulay2(output);
	} else {
		output << Msym;
	}
	output << endl;
}

void powerChoice(vector<string> &variables, vector<string> &relations,vector<vector<string>> &matrix, unsigned int characteristic, unsigned int symPower, unsigned int extPower,  ofstream &output, const program_options & opt = noOptions)
{
	if(characteristic == 0)
	{
		if(opt.forceQ)
		{
			PolyRing<mpq_class> mB(variables.size(), variables, degrevlex);
			PolyRing<mpq_class> QuotientRing(&mB, relations);
			DensePolyMatrix<mpq_class> M(&QuotientRing, matrix);
			powerMatrix(M, characteristic, symPower, extPower, output, opt);
		} else {
			ZRing mB(variables.size(), variables, degrevlex);

			vector<Zpoly > polrels;
			if(ZViable(mB,relations,polrels))
			{
				ZRing QuotientRing(&mB, polrels);

				ZMatrix M(&QuotientRing, matrix);
				powerMatrix(M,characteristic, symPower, extPower,output, opt);
			} else {
				PolyRing<mpq_class> mB(variables.size(), variables, degrevlex);
				PolyRing<mpq_class> QuotientRing(&mB, relations);
				DensePolyMatrix<mpq_class> M(&QuotientRing, matrix);
				powerMatrix(M,characteristic, symPower, extPower,output, opt);
			}
		}

	} else {
		numbermodulo::setModulus(characteristic);

		PolyRing<numbermodulo> mB(variables.size(), variables, degrevlex);
		PolyRing<numbermodulo> QuotientRing(&mB, relations);
		DensePolyMatrix<numbermodulo> M(&QuotientRing, matrix);
		powerMatrix(M, characteristic, symPower, extPower,output, opt);
	}
}

template <typename Scalar>
void kernelMatrix( DensePolyMatrix<Scalar> &M, unsigned int characteristic,  unsigned int symPower, unsigned int extPower, unsigned int twist,ofstream &output,  const program_options &opt = noOptions){
	M.getRing()->setOptions(opt);
	if (opt.verbosity >= 1 )
	{
		if(characteristic == 0)
		{
			cout << "Computing Kernel of Sym^" << symPower << "(Ext^" << extPower << "(M)) with twist " << twist << endl;
		} else {
			cout << "Computing Kernel of F*^" << symPower << "(Ext^" << extPower << "(M)) with twist " << twist << endl;
		}
	}
	if (opt.verbosity >= 2)
	{
		cout << "Input Matrix:" << endl;
		cout << M << endl;
	}
	SparsePolyMatrix<Scalar> Msparse = M;
	SparsePolyMatrix<Scalar> Mext = M.extMatrix(extPower);
	SparsePolyMatrix<Scalar> Msym;
	if(characteristic == 0)
	{
		Msym = Mext.symMatrix(symPower);
	} else {
		Msym = Mext.powerMatrix(symPower);
	}

	SparsePolyMatrix<Scalar> kernel;
	intest(Msym,twist,kernel,opt);

	if (opt.verbosity >= 2)
	{
		cout << "Kernel of rank " << kernel.getCols() << endl;
		cout << kernel << endl;
	}
	if(opt.outputM2)
	{
		output << kernel.getCols() << endl;
		kernel.printForMacaulay2(output);
	} else {
		output << kernel.getCols() << endl;
		output << kernel;
	}
	output << endl;
}

void kernelChoice(vector<string> &variables, vector<string> &relations,vector<vector<string>> &matrix, unsigned int characteristic, unsigned int symPower, unsigned int extPower,  unsigned int twist, ofstream &output, const program_options & opt = noOptions)
{
	if(characteristic == 0)
	{
		if(opt.forceQ)
		{
			PolyRing<mpq_class> mB(variables.size(), variables, degrevlex);
			PolyRing<mpq_class> QuotientRing(&mB, relations);
			DensePolyMatrix<mpq_class> M(&QuotientRing, matrix);
			kernelMatrix(M, characteristic, symPower, extPower, twist, output, opt);
		} else {
			ZRing mB(variables.size(), variables, degrevlex);

			vector<Zpoly > polrels;
			if(ZViable(mB,relations,polrels))
			{
				ZRing QuotientRing(&mB, polrels);

				ZMatrix M(&QuotientRing, matrix);
				kernelMatrix(M,characteristic, symPower, extPower,twist, output, opt);
			} else {
				PolyRing<mpq_class> mB(variables.size(), variables, degrevlex);
				PolyRing<mpq_class> QuotientRing(&mB, relations);
				DensePolyMatrix<mpq_class> M(&QuotientRing, matrix);
				kernelMatrix(M,characteristic, symPower, extPower,twist, output, opt);
			}
		}

	} else {
		numbermodulo::setModulus(characteristic);

		PolyRing<numbermodulo> mB(variables.size(), variables, degrevlex);
		PolyRing<numbermodulo> QuotientRing(&mB, relations);
		DensePolyMatrix<numbermodulo> M(&QuotientRing, matrix);
		kernelMatrix(M, characteristic, symPower, extPower,twist, output, opt);
	}
}

//Compute Jacobi-Taylor-Matrices and the resulting operators
template <typename Scalar>
void jacobiTaylor(PolyRing<Scalar> ring, unsigned int maxd, unsigned int maxn,ofstream &output, const program_options & opt = noOptions)
{

	if(ring.getGB().size() > 1)
	{
		cout << "We can only compute the Jacobi Taylor Matrix for one relation. Please input only one polynomial." << endl;
		return;
	}
	Poly<Scalar> relation = ring.getGB()[0];

	ring.setOptions(opt);

	if(maxd == 0)
	{

		unsigned int rank = 0;

		for (unsigned int i = 1; i < maxn; i++)
		{
			if (opt.verbosity >= 2)
			{
				cout << "Jacobi-Taylor Matrix for n=" << i << endl;
			}
			output << "Jacobi-Taylor Matrix for n=" << i << endl;

			DensePolyMatrix<Scalar> JT = DensePolyMatrix<Scalar>::jacobiTaylor(&ring, relation, i);

			if (opt.verbosity >= 2)
			{
				cout << JT << endl;
			}
			output << JT << endl;

			clock_t t = clock();

			for (int d = (int) relation.degree(); d < (int) relation.degree() + 1; d++)
			{
				DensePolyMatrix<Scalar> kernel = JT.kernel(d,true,opt);

				if (opt.verbosity >= 3)
				{
					cout << "Kernel(" << d - (int) relation.degree() << "):" << endl;
					kernel.printAsOperator(cout);
				}
				output << "Kernel(" << d - (int) relation.degree() << "):" << endl;
				kernel.printAsOperator(output);

				/*
				//Future improvements would only compute unitary operators

				if (kernel.getCols() > 0)
				{

					vector<DensePolyMatrix<Scalar>> list;
					for (unsigned int i = 0; i < kernel.getRows(); i++)
					{
						list.push_back(DensePolyMatrix<Scalar>::irrelevantPower(ring, d - kernel.getOutdegrees()[i]));
					}
					DensePolyMatrix<Scalar> Mn = directsum(list);

					DensePolyMatrix<Scalar> kernelmod = kernel % Mn;

					if (opt.verbosity >= 3)
					{
						cout << "Kernel(" << d - (int) relation.degree() << "):" << endl;
						kernelmod.printAsOperator(cout);
					}
					output << "Kernel(" << d - (int) relation.degree() << "):" << endl;
					kernelmod.printAsOperator(output);

					if (d == (int) relation.degree() - (int) i)
					{
						rank += kernelmod.getCols();
					}
				}*/
				rank += kernel.getCols();

			}

			if (opt.verbosity >= 2)
			{
				cout << "Rank: " << rank + 1 << ", cum. hilbert function:" << ring.getIndexOfBasisDegStart(i + 1) << endl;

				cout << "Norm: " << rank + 1 << "/" << ring.getMonN(i + 1) << " = " << (float) (rank + 1) / (float) (ring.getMonN(i + 1)) << endl;

				t = clock() - t;
				cout << "Time spent: " << t / CLOCKS_PER_SEC << "s (" << t << " internal ticks)" << endl << endl;
			}
		}
		cout << "Norm: " << rank + 1 << "/" << ring.getMonN(maxn + 1) << " = " << (float) (rank + 1) / (float) (ring.getMonN(maxn + 1)) << endl;
	} else {
		vector<vector<int>> results;

		if(opt.verbosity >= 2)
		{
			cout << "Computing Jacobi-Taylor Table...";
		}

		output << "n\\d";
		for (unsigned int d = relation.degree(); d <= maxd; d++)
		{
			output << "\t" << d - relation.degree();
		}
		output << endl;

		for (unsigned int i = 1; i <= maxn; i++)
		{

			output << i;

			DensePolyMatrix<Scalar> JT = DensePolyMatrix<Scalar>::jacobiTaylor(&ring, relation, i);

			vector<int> row;

			for (unsigned int d = relation.degree(); d <= maxd; d++)
			{
				DensePolyMatrix<Scalar> kernel = JT.kernel(d,true,opt);

				row.push_back(kernel.getCols());

				int col = d - relation.degree();

				if (i == 1)
				{
					output << "\t" << row[col];
				}
				else
				{
					output << "\t" << row[col] - results[i - 2][col];
				}

			}
			output << endl;

			results.push_back(row);
		}

	}
}

void jacobiTaylorChoice(vector<string> &variables, vector<string> &relations, unsigned int characteristic,unsigned int maxd, unsigned int maxn, ofstream &output, const program_options & opt = noOptions)
{
	if(characteristic == 0)
	{
		if(opt.forceQ)
		{
			PolyRing<mpq_class> mB(variables.size(), variables, degrevlex);
			PolyRing<mpq_class> QuotientRing(&mB, relations);
			jacobiTaylor(QuotientRing, maxd, maxn, output, opt);
		} else {
			ZRing mB(variables.size(), variables, degrevlex);

			vector<Zpoly > polrels;
			if(ZViable(mB,relations,polrels))
			{
				ZRing QuotientRing(&mB, polrels);
				jacobiTaylor(QuotientRing, maxd, maxn, output, opt);
			} else {
				PolyRing<mpq_class> mB(variables.size(), variables, degrevlex);
				PolyRing<mpq_class> QuotientRing(&mB, relations);
				jacobiTaylor(QuotientRing, maxd, maxn, output, opt);
			}
		}

	} else {
		numbermodulo::setModulus(characteristic);

		PolyRing<numbermodulo> mB(variables.size(), variables, degrevlex);
		PolyRing<numbermodulo> QuotientRing(&mB, relations);
		jacobiTaylor(QuotientRing, maxd, maxn, output, opt);
	}
}

bool readuint(string line, string keyword, unsigned int &val)
{
	std::size_t pos = line.find(keyword);
	if (pos != std::string::npos)
	{
		string num = line.substr(pos + keyword.length());

		val = std::stoi(num);
		return true;
	} else{
		return false;
	}
}

bool readstrings(string line, string keyword, vector<string> &values)
{
	std::size_t pos = line.find(keyword);
	bool ret = (pos != std::string::npos);
	while (pos != std::string::npos)
	{
		std::size_t p1 = line.find("\"", pos);
		std::size_t p2 = line.find("\"", p1 + 1);
		if (p1 != std::string::npos && p2 != std::string::npos)
		{
			string var = line.substr(p1 + 1, p2 - (p1 + 1));
			values.push_back(var);
		}
		pos = line.find(",", p2);
	}
	return ret;
}

bool readobjects(string line, string keyword, vector<string> &values)
{
	std::size_t pos = line.find(keyword);
	bool ret = (pos != std::string::npos);
	std::size_t p1 = pos + keyword.length();
	while (pos != std::string::npos)
	{

		std::size_t p2 = line.find(",", p1);
		if (p2 == std::string::npos)
		{
			pos = p2;
			p2 = line.length();
		}

		string rel = line.substr(p1, p2 - p1);
		values.push_back(rel);

		p1 = p2 + 1;
	}
	return ret;
}

bool readmatrix(string line, string keyword, vector<vector<string> > &values)
{
	std::size_t pos = line.find(keyword);
	bool ret = (pos != std::string::npos);

	std::size_t y1 = line.find("{", pos);
	while (pos != std::string::npos && y1 != std::string::npos)
	{
		vector<string> row;
		std::size_t x1 = line.find("{", y1 + 1);
		if (x1 == std::string::npos)
		{
			break;
		}
		std::size_t x2p = line.find("}", x1 + 1);
		while (x1 != std::string::npos)
		{
			std::size_t x2 = line.find(",", x1 + 1);
			if (x2 != std::string::npos)
			{
				x2 = min(x2, x2p);
			}
			else
			{
				x2 = x2p;
			}
			string entry = line.substr(x1 + 1, x2 - (x1 + 1));
			row.push_back(entry);
			if (x2 == x2p)
			{
				x1 = std::string::npos;
			}
			else
			{
				x1 = x2;
			}
		}
		values.push_back(row);
		y1 = x2p + 1;
	}
	return ret;
}

bool readsemistability(ifstream &myfile, program_options &opt)
{
	unsigned int characteristic = 0;
	vector<string> variables;
	vector<string> strrelations;
	vector<vector<string> > strmatrix;
	string keywords[] =
	{ "characteristic:", "variables:", "relations:", "matrix:" };
	vector<int> hasRead(4,0);

	string line;
	while (getline(myfile, line))
	{
		//read characteristic
		hasRead[0] |= readuint(line,keywords[0],characteristic);

		//read variables
		hasRead[1] |= readstrings(line,keywords[1],variables);

		//read relations
		hasRead[2] |= readobjects(line,keywords[2],strrelations);

		//read matrix
		hasRead[3] |= readmatrix(line,keywords[3],strmatrix);
	}
	myfile.close();
	for(unsigned int i = 0; i < hasRead.size(); i++)
	{
		if(hasRead[i] != 1)
		{
			cerr << "There has to be exactly one line with " << keywords[i] << "in the input file.";
			return false;
		}
	}

	ofstream outfile(opt.output_file);
	if (!outfile.is_open())
	{
		cerr << "Unable to open output file";
		return false;
	}

	bool ret = semistabilityChoice(variables, strrelations,strmatrix,characteristic,outfile,opt);

	outfile.close();
	return ret;
}


bool readpowers(ifstream &myfile, program_options &opt)
{
	string line;
	unsigned int characteristic = 0;
	vector<string> variables;
	vector<string> strrelations;
	vector<vector<string> > strmatrix;
	unsigned int extPower = 1;
	unsigned int symPower = 1;

	string keywords[] =
	{ "characteristic:", "variables:", "relations:", "matrix:", "exteriorpower:", "spower:" };
	vector<int> hasRead(6, 0);

	while (getline(myfile, line))
	{
		//read characteristic
		hasRead[0] |= readuint(line, keywords[0], characteristic);

		//read variables
		hasRead[1] |= readstrings(line, keywords[1], variables);

		//read relations
		hasRead[2] |= readobjects(line, keywords[2], strrelations);

		//read matrix
		hasRead[3] |= readmatrix(line, keywords[3], strmatrix);

		hasRead[4] |= readuint(line, keywords[4], extPower);

		hasRead[5] |= readuint(line, keywords[5], symPower);
	}
	myfile.close();
	for (unsigned int i = 0; i < hasRead.size(); i++)
	{
		if (hasRead[i] != 1)
		{
			cerr << "There has to be exactly one line with " << keywords[i] << "in the input file.";
			return false;
		}
	}

	ofstream outfile(opt.output_file);
	if (!outfile.is_open())
	{
		cerr << "Unable to open output file";
		return false;
	}

	powerChoice(variables, strrelations,strmatrix,characteristic,symPower,extPower,outfile,opt);


	outfile.close();
	return true;
}

bool readkernel(ifstream &myfile, program_options &opt)
{
	string line;
	unsigned int characteristic = 0;
	vector<string> variables;
	vector<string> strrelations;
	vector<vector<string> > strmatrix;
	unsigned int extPower = 1;
	unsigned int symPower = 1;
	unsigned int twist = 0;

	string keywords[] =
	{ "characteristic:", "variables:", "relations:", "matrix:", "exteriorpower:", "spower:" , "twist:"};
	vector<int> hasRead(7, 0);

	while (getline(myfile, line))
	{
		//read characteristic
		hasRead[0] |= readuint(line, keywords[0], characteristic);

		//read variables
		hasRead[1] |= readstrings(line, keywords[1], variables);

		//read relations
		hasRead[2] |= readobjects(line, keywords[2], strrelations);

		//read matrix
		hasRead[3] |= readmatrix(line, keywords[3], strmatrix);

		hasRead[4] |= readuint(line, keywords[4], extPower);

		hasRead[5] |= readuint(line, keywords[5], symPower);

		hasRead[6] |= readuint(line, keywords[6], twist);
	}
	myfile.close();
	for (unsigned int i = 0; i < hasRead.size(); i++)
	{
		if (hasRead[i] != 1)
		{
			cerr << "There has to be exactly one line with " << keywords[i] << "in the input file.";
			return false;
		}
	}

	ofstream outfile(opt.output_file);
	if (!outfile.is_open())
	{
		cerr << "Unable to open output file";
		return false;
	}

	kernelChoice(variables, strrelations,strmatrix,characteristic,symPower,extPower,twist,outfile,opt);


	outfile.close();

	return true;
}

bool readjacobiTaylor(ifstream &myfile, program_options &opt)
{
	string line;
	unsigned int characteristic = 0;
	vector<string> variables;
	vector<string> strrelations;
	unsigned int maxd = 0;
	unsigned int maxn = 0;

	string keywords[] =
	{ "characteristic:", "variables:", "relations:", "maxn:",  "maxd:"};
	vector<int> hasRead(4, 0);

	while (getline(myfile, line))
	{
		//read characteristic
		hasRead[0] |= readuint(line, keywords[0], characteristic);

		//read variables
		hasRead[1] |= readstrings(line, keywords[1], variables);

		//read relations
		hasRead[2] |= readobjects(line, keywords[2], strrelations);

		//read tabledimensions
		hasRead[3] |= readuint(line, keywords[3], maxn);
		hasRead[4] |= readuint(line, keywords[4], maxd);

	}
	myfile.close();
	for (unsigned int i = 0; i < hasRead.size(); i++)
	{
		if (hasRead[i] != 1)
		{
			cerr << "There has to be exactly one line with " << keywords[i] << "in the input file.";
			return false;
		}
	}

	ofstream outfile(opt.output_file);
	if (!outfile.is_open())
	{
		cerr << "Unable to open output file";
		return false;
	}

	opt.computeFullKernel = true;
	jacobiTaylorChoice(variables, strrelations,characteristic, maxd, maxn,outfile,opt);


	outfile.close();

	return true;
}

//read the first line of the input file and hand over to the appropriate subroutine
bool readandexecutefile(program_options &opt)
{
	string line;
	ifstream myfile(opt.input_file);
	if (myfile.is_open())
	{
		getline(myfile, line);
		if (line.find("semistability") != std::string::npos)
		{ //User wants us to determine semistability
			readsemistability(myfile,opt);

		} else if(line.find("powers") != std::string::npos)
		{
			readpowers(myfile,opt);
		} else if(line.find("kernel") != std::string::npos)
		{
			readkernel(myfile,opt);
		} else if(line.find("jacobiTaylor") != std::string::npos)
		{
			readjacobiTaylor(myfile,opt);
		}
		else
		{
			cerr << "Could not determine task" << endl;
			myfile.close();
			return false;
		}


		return true;
	}

	else
	{
		cerr << "Unable to open input file" << endl;
		return false;
	}
}

template<class T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
	copy(v.begin(), v.end(), ostream_iterator<T>(os, " "));
	return os;
}

//main function, in which we read the command line options and then hand over to reading the input file
int main(int ac, char* av[])
{
	try
	{
		program_options opt;
		po::options_description generic("Generic options");
		generic.add_options()("version,v", "print the version string")("help", "produce a help message")("config,c", po::value<string>(&opt.config_file)->default_value("sheaf_stability.cfg"), "name of a configuration file");

		po::options_description config("Configuration");
		config.add_options()
				("verbosity", po::value<int>(&opt.verbosity)->default_value(1), "Verbosity Level 0-3, Default: 1, O is silent unless there is an error")
				("input-file", po::value<string>(&opt.input_file)->default_value("input.txt"), "input file")
				("output-file", po::value<string>(&opt.output_file)->default_value("sheaf_out.txt"), "output file")
				("exterior-powers,e", po::bool_switch(&opt.exteriorPowers)->default_value(false), "Using exterior powers in the algorithm")
				("stop-unstable,s", po::bool_switch(&opt.stopIfUnstable)->default_value(false), "Stop immediately if a destabilizing section has been found")
				("linear-progression,l", po::bool_switch(&opt.linearPowerProgresison)->default_value(false), "Linear progression of powers instead of fast lookahead")
		 	 	("tarjan,t", po::bool_switch(&opt.useTarjan)->default_value(false), "Use Tarjan's algorithm before the Gauß algorithm")
				("forceQ,q", po::bool_switch(&opt.forceQ)->default_value(false), "Force Computations over Q (in characteristic 0). By Default Z is used where possible")
				("output-for-M2,m", po::bool_switch(&opt.outputM2)->default_value(false), "Write output matrices in a format easily readable by Macaulay2")
				("compute-full-kernel,f", po::bool_switch(&opt.computeFullKernel)->default_value(false), "Always compute the full kernels and not just whether the kernel is nonzero")
				("long-form-polys,p", po::bool_switch(&opt.longFormPolynomials)->default_value(false), "Output Polynomials in a less condensed format");

		po::options_description cmdline_options;
		cmdline_options.add(generic).add(config);

		po::options_description config_file_options;
		config_file_options.add(config);

		po::options_description visible("Allowed options");
		visible.add(generic).add(config);

		po::positional_options_description p;
		p.add("input-file", 1);
		p.add("output-file", 1);

		po::variables_map vm;
		store(po::command_line_parser(ac, av).options(cmdline_options).positional(p).run(), vm);
		notify(vm);

		ifstream ifs(opt.config_file.c_str());
		if (ifs)
		{
			store(parse_config_file(ifs, config_file_options), vm);
			notify(vm);
		}

		if (vm.count("help"))
		{
			cout << visible << "\n";
			return 0;
		}

		if (vm.count("version"))
		{
			cout << "Sheaf stability computer, version " << version1 << "." << version2 << "." << version3 << "\n";
			return 0;
		}

		return readandexecutefile(opt);
	} catch (exception& e)
	{
		cout << e.what() << "\n";
		return 1;
	}
	return 0;
}

