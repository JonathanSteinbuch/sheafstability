//============================================================================
// Name        : SemiStabilityJob.hpp
// Author      : Jonathan Steinbuch
//============================================================================

#include "JobType.hpp"

#ifndef SRC_SEMISTABILITYJOB_HPP_
#define SRC_SEMISTABILITYJOB_HPP_


template <typename Scalar>
class SemistabilityJob : public JobType<Scalar>
{
	int intest(SparsePolyMatrix<Scalar> Matrix, int degreeOfInterest, SparsePolyMatrix<Scalar> &ker, program_options opt);

public:

	bool doJob(PolyRing<Scalar> &Ring, DensePolyMatrix<Scalar> &M,const JobOptionsType &jobOptions, const program_options &opt = noOptions);
};

#endif /* SRC_SEMISTABILITYJOB_HPP_ */
