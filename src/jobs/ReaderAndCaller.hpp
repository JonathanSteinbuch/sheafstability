//============================================================================
// Name        : ReaderAndCaller.hpp
// Author      : Jonathan Steinbuch
//============================================================================

#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "../helpers/Helpers.hpp"

#include "JobType.hpp"

#include "JacobiTaylorJob.hpp"
#include "KernelJob.hpp"
#include "PowersJob.hpp"
#include "SemistabilityJob.hpp"

#ifndef SRC_JOBS_READERANDCALLER_HPP_
#define SRC_JOBS_READERANDCALLER_HPP_

class ReaderAndCaller
{
	const program_options opt;

	bool ZViable(ZRing &mB,  const vector<string> &relations, vector<Zpoly > &polrels);

	bool readandexecutefile();

	template <typename Scalar>
	bool startJob(DensePolyMatrix<Scalar> &M, JobOptionsType& jobOptions);

	bool ringChoice(JobOptionsType& jobOptions);

public:
	ReaderAndCaller(const program_options &opt) : opt(opt){
		readandexecutefile();
	};

	virtual ~ReaderAndCaller();
};


#endif /* SRC_JOBS_READERANDCALLER_HPP_ */
