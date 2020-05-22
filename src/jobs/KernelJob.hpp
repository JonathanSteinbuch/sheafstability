/*
 * KernelJob.hpp
 *
 *  Created on: May 13, 2020
 *      Author: jonathan
 */

#ifndef SRC_JOBS_KERNELJOBHPP_
#define SRC_JOBS_KERNELJOBHPP_

#include "JobType.hpp"

template <typename Scalar>
class KernelJob: public JobType<Scalar>
{

public:

	bool doJob(PolyRing<Scalar> &Ring, DensePolyMatrix<Scalar> &M,const JobOptionsType &jobOptions, const program_options &opt = noOptions);
};

#endif /* SRC_JOBS_KERNELJOBHPP_ */
