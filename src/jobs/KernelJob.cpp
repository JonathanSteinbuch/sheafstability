/*
 * KernelJob.cpp
 *
 *  Created on: May 13, 2020
 *      Author: jonathan
 */

#include "KernelJob.hpp"


template <typename Scalar>
bool KernelJob<Scalar>::doJob(PolyRing<Scalar> &Ring, DensePolyMatrix<Scalar> &M,const JobOptionsType &jobOptions, const program_options &opt)
{
	M.getRing()->setOptions(opt);
	if (opt.verbosity >= 1 )
	{
		if(jobOptions.characteristic == 0)
		{
			cout << "Computing Kernel of Sym^" << jobOptions.symPower << "(Ext^" << jobOptions.extPower << "(M)) with twist " << jobOptions.twist << endl;
		} else {
			cout << "Computing Kernel of F*^" << jobOptions.symPower << "(Ext^" << jobOptions.extPower << "(M)) with twist " << jobOptions.twist << endl;
		}
	}
	if (opt.verbosity >= 2)
	{
		cout << "Input Matrix:" << endl;
		cout << M << endl;
	}
	SparsePolyMatrix<Scalar> Msparse = M;
	SparsePolyMatrix<Scalar> Mext = M.extMatrix(jobOptions.extPower);
	SparsePolyMatrix<Scalar> Msym;
	if(jobOptions.characteristic == 0)
	{
		Msym = Mext.symMatrix(jobOptions.symPower);
	} else {
		Msym = Mext.powerMatrix(jobOptions.symPower);
	}

	SparsePolyMatrix<Scalar> kernel;
	kernel = Msym.kernel(jobOptions.twist,true,opt);

	if (opt.verbosity >= 2)
	{
		cout << "Kernel of rank " << kernel.getCols() << endl;
		cout << kernel << endl;
	}
	if(opt.outputM2)
	{
		*(jobOptions.output) << kernel.getCols() << endl;
		kernel.printForMacaulay2(*(jobOptions.output));
	} else if(opt.outputLatex){
		*(jobOptions.output) << kernel.getCols() << endl;
		kernel.printForLatex(*(jobOptions.output));
	} else {
		*(jobOptions.output) << kernel.getCols() << endl;
		*(jobOptions.output) << kernel;
	}
	*(jobOptions.output) << endl;
	return (kernel.getCols() > 0);
}

template class KernelJob<mpz_class>;
template class KernelJob<mpq_class>;
template class KernelJob<numbermodulo>;
