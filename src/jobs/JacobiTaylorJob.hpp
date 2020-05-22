/*
 * JacobiTaylorJob.hpp
 *
 *  Created on: May 13, 2020
 *      Author: jonathan
 */

#ifndef SRC_JOBS_JACOBITAYLORJOB_HPP_
#define SRC_JOBS_JACOBITAYLORJOB_HPP_

#include "JobType.hpp"

template <typename Scalar>
class JacobiTaylorJob: public JobType<Scalar>
{

public:

	bool doJob(PolyRing<Scalar> &Ring, DensePolyMatrix<Scalar> &M,const JobOptionsType &jobOptions, const program_options &opt = noOptions);
};

#endif /* SRC_JOBS_JACOBITAYLORJOB_HPP_ */
