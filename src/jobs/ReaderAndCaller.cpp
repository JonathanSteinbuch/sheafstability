//============================================================================
// Name        : ReaderAndCaller.cpp
// Author      : Jonathan Steinbuch
//============================================================================


#include "ReaderAndCaller.hpp"


datatype numbermodulo::modulus;
invTable numbermodulo::lookup;


template<typename Scalar>
inline bool ReaderAndCaller::startJob(DensePolyMatrix<Scalar> &M, JobOptionsType& jobOptions)
{

	JobType<Scalar>* Job;
	switch(jobOptions.type)
	{
	case 0: //semistability
		Job = new SemistabilityJob<Scalar>;
	break;
	case 1: //powers
		Job = new PowersJob<Scalar>;
	break;
	case 2: //kernel
		Job = new KernelJob<Scalar>;
	break;
	default:
		cerr << "Internal Job Type confusion.";
		return false;
	}

	ofstream outfile(opt.output_file);
	if (!outfile.is_open())
	{
		cerr << "Unable to open output file";
		return false;
	}
	jobOptions.output = &outfile;

	bool retVal = Job->doJob(*(M.getRing()),M,jobOptions,opt);

	outfile.close();

	delete Job;
	return retVal;
}

bool ReaderAndCaller::ZViable(ZRing &mB,  const vector<string> &relations, vector<Zpoly > &polrels){
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

bool ReaderAndCaller::ringChoice(JobOptionsType& jobOptions)
{
	vector<string> &variables = jobOptions.variables;
	vector<string> &relations = jobOptions.strrelations;
	vector<vector<string>> &matrix = jobOptions.strmatrix;

	if(jobOptions.characteristic == 0)
	{
		if(opt.forceQ)
		{
			PolyRing<mpq_class> mB(variables.size(), variables, degrevlex);
			PolyRing<mpq_class> QuotientRing(&mB, relations);
			DensePolyMatrix<mpq_class> M(&QuotientRing, matrix);
			return startJob(M,jobOptions);
		} else {
			ZRing mB(variables.size(), variables, degrevlex);

			vector<Zpoly > polrels;
			if(ZViable(mB,relations,polrels))
			{
				ZRing QuotientRing(&mB, polrels);

				ZMatrix M(&QuotientRing, matrix);
				return startJob(M,jobOptions);
			} else {
				PolyRing<mpq_class> mB(variables.size(), variables, degrevlex);
				PolyRing<mpq_class> QuotientRing(&mB, relations);
				DensePolyMatrix<mpq_class> M(&QuotientRing, matrix);
				return startJob(M,jobOptions);
			}
		}

	} else {
		numbermodulo::setModulus(jobOptions.characteristic);

		PolyRing<numbermodulo> mB(variables.size(), variables, degrevlex);
		PolyRing<numbermodulo> QuotientRing(&mB, relations);
		DensePolyMatrix<numbermodulo> M(&QuotientRing, matrix);
		return startJob(M,jobOptions);
	}
}

bool ReaderAndCaller::readandexecutefile()
{
	string line;
	ifstream inputfile(opt.input_file);

	if (inputfile.is_open())
	{
		getline(inputfile, line);
		unsigned int type;
		for(type = 0; type < DescriptorList.size(); type++)
		{
			if (line.find(DescriptorList[type].name) != std::string::npos)
			{
				break;
			}
		}
		if(type == DescriptorList.size())
		{
			cerr << "Could not determine task" << endl;
			inputfile.close();
			return false;
		}
		bool retVal;
		JobOptionsType jobOptions(type,&inputfile,retVal);

		inputfile.close();

		if(!retVal)
		{
			return false;
		}
		retVal = ringChoice(jobOptions);

		return retVal;
	}

	else
	{
		cerr << "Unable to open input file" << endl;
		return false;
	}
}

ReaderAndCaller::~ReaderAndCaller()
{
	// TODO Auto-generated destructor stub
}

