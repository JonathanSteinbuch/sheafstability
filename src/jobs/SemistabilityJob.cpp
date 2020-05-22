//============================================================================
// Name        : SemiStabilityJob.cpp
// Author      : Jonathan Steinbuch
//============================================================================

#include "SemistabilityJob.hpp"

template <typename Scalar>
int SemistabilityJob<Scalar>::intest(SparsePolyMatrix<Scalar> Matrix, int degreeOfInterest, SparsePolyMatrix<Scalar> &ker, program_options opt)
{
	SparseScalarMatrix<Scalar> inMat(Matrix, degreeOfInterest);

	if (opt.verbosity >= 1)
	{

		cout << "Computing Kernel of " << Matrix.getRows() << "x" << Matrix.getCols() << "-Matrix (internally " << inMat.getRows() << "x" << inMat.getCols() << ")...";
		cout << endl;
	}

	if (opt.verbosity >= 4)
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

template <typename Scalar>
bool SemistabilityJob<Scalar>::doJob(PolyRing<Scalar> &Ring, DensePolyMatrix<Scalar> &M,const JobOptionsType &jobOptions, const program_options &opt)
{
	bool result = false;

	Ring.setOptions(opt);

	assert(M.getRows() == 1);

	if (opt.verbosity >= 1)
	{

		cout << "We compute if the Syzygy sheaf given by M is semistable over the projective variety X, where:" << endl;
		string field;
		if (jobOptions.characteristic == 0)
		{
			cout << "X is given by \\Q[";
		} else
		{
			cout << "X is given by \\F_" << jobOptions.characteristic << "[";
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

	PolyRing<mpq_class> kX(1,{ "n" }, degrevlex);
	kX.setOptions(opt);
	Poly<mpq_class> hP(&kX);
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

	if (!Ring.isSmooth())
	{
		cerr << "Curve not smooth! Computation aborted." << endl;
		return false;
	}
	else if (opt.verbosity >= 1)
	{
		cout << "Curve is smooth: Check." << endl;
	}

	unsigned int reldeg = hP.leadingCoefficient().get_num().get_ui();
	assert(hP.leadingCoefficient().get_den() == 1);
	int genus = 1;
	if (hP.supportSize() > 1)
	{
		genus = (1 - hP.trailingCoefficient().get_num().get_si());
	}

	if (opt.verbosity >= 1)
	{
		cout << "Genus of Curve: " << genus << endl << endl;
	}

	auto Mgens = M.getRowEntries(0);
	for (auto & p : Mgens)
	{
		p.changeRing(Ring.getBaseRing());
	}
	Mgens.insert(Mgens.end(), Ring.getGB().begin(), Ring.getGB().end());
	Ring.getBaseRing()->computeGroebnerBasis(Mgens);
	Poly<mpq_class> MhP(&kX);
	Ring.getBaseRing()->quotientHilbertPoly(MhP, Mgens);
	if (opt.verbosity >= 1)
	{
		cout << "The Hilbert polynomial of M is P(n)=" << MhP << endl;
		cout << endl;
	}

	assert(MhP.degree() == 0);
	int idealdegree = -MhP.leadingCoefficient().get_num().get_si();

	int rank = M.getCols() - M.getRows();
	int degree = M.indegree()*reldeg - idealdegree;
	double slope = (double) degree / (double) rank;

	unsigned int threshold = genus - 1 + reldeg;

	SparsePolyMatrix<Scalar> sparseM = SparsePolyMatrix<Scalar>(M);


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

		if (opt.verbosity >= 1)
		{
			cout << "Exterior Power: " << extPower << ", ";
			cout << extMatrix.getRows() << "x" << extMatrix.getCols() << "-Matrix" << endl;
		}

		int gc;
		int n;
		if(opt.exteriorPowers)
		{
			gc = __gcd(rank, -extPower*degree);
			n = rank/gc;
		} else {
			gc = __gcd(rank, -degree);
			n = (rank - 1)*rank/gc;
		}

		unsigned int actualThreshold = threshold * n + 1;
		double extSlope = slope * extPower;

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

		while (currmin < currmax)
		{

			if (jobOptions.characteristic == 0)
			{
				if (opt.verbosity >= 1)
				{
					cout << "Symmetric Power: " << sPower << endl;
				}
				sMatrix = extMatrix.symMatrix(sPower);
			}
			else
			{
				if (opt.verbosity >= 1)
				{
					cout << "Frobenius Power: " << sPower << endl;
				}
				if(sPower > oldsPower)
				{
					sMatrix = sMatrix.powerMatrix(sPower/oldsPower);
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
			int twist = ((int) ceil(-sSlope/reldeg) - 1);

			if (opt.verbosity >= 1)
			{
				cout << "Slope: " << sSlope << endl;
				cout << "Twist: " << twist << endl;
			}

			clock_t t = clock();

			SparsePolyMatrix<Scalar> ker = sMatrix.kernel(twist,true,opt);

			int kernelRank = ker.getCols();

			if (opt.verbosity >= 1)
			{
				t = clock() - t;
				cout << "Clock: " << t << "ticks or " << t / CLOCKS_PER_SEC << "s" << endl;
				cout << "Rank at degree of Interest: " << kernelRank << endl;
			}
			if (kernelRank > 0)
			{
				if (opt.verbosity >= 1)
				{
					cout << "Not Semistable!" << endl;
				}
				*(jobOptions.output) << 1 << endl;
				*(jobOptions.output) << "{" << "ExteriorPower => " << extPower << ", SymmetricPower => " << sPower << ", Twist => " << twist << "}" << endl;
				if(opt.verbosity >= 2)
				{
					cout << "Kernel:" << endl;
					cout << ker << endl << endl;
				}
				if(opt.verbosity >= 3)
				{
					if(opt.outputM2)
					{
						ker.printForMacaulay2(*(jobOptions.output));
					} else if(opt.outputLatex){
						ker.printForLatex(*(jobOptions.output));
					} else {
						*(jobOptions.output) << ker;
					}
				}
				*(jobOptions.output) << endl;
				result = false;
				if(jobOptions.characteristic == 0)
				{
					currmax = sPower - 1;
				} else {
					currmax = sPower / jobOptions.characteristic;
				}
				if(opt.stopIfUnstable)
				{
					return result;
				}
			}
			else if (extPower == rank - 1 && sPower == actualThreshold && kernelRank == 0)
			{
				if (opt.verbosity >= 1)
				{
					cout << "Semistable!" << endl;
				}
				*(jobOptions.output) << 0 << endl;
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

			if (jobOptions.characteristic == 0)
			{
				if(opt.linearPowerProgresison)
				{
					sPower++;
				} else{
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
				sPower *= jobOptions.characteristic;
				if (sPower == 1)
				{
					sPower = 0;
				}
			}

		}
	}
	return result;
}


template class SemistabilityJob<mpz_class>;
template class SemistabilityJob<mpq_class>;
template class SemistabilityJob<numbermodulo>;
