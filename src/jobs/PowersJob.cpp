/*
 * PowersJob.cpp
 *
 *  Created on: May 13, 2020
 *      Author: jonathan
 */

#include "PowersJob.hpp"



template <typename Scalar>
bool PowersJob<Scalar>::doJob(PolyRing<Scalar> &Ring, DensePolyMatrix<Scalar> &M,const JobOptionsType &jobOptions, const program_options &opt){
	M.getRing()->setOptions(opt);

	if (opt.verbosity >= 1 )
	{
		if(jobOptions.characteristic == 0)
		{
			cout << "Computing Sym^" << jobOptions.symPower << "(Ext^" << jobOptions.extPower << "(M))." << endl;
		} else {
			cout << "Computing F*^" << jobOptions.symPower << "(Ext^" << jobOptions.extPower << "(M))." << endl;
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

	if (opt.verbosity >= 2)
	{
		cout << "Power Matrix:" << endl;
		cout << Msym << endl;
	}
	if(opt.outputM2)
	{
		Msym.printForMacaulay2(*(jobOptions.output));
	} else if(opt.outputLatex){
		Msym.printForLatex(*(jobOptions.output));
	} else {
		*(jobOptions.output) << Msym;
	}
	*(jobOptions.output) << endl;

	return true;
}

template class PowersJob<mpz_class>;
template class PowersJob<mpq_class>;
template class PowersJob<numbermodulo>;
