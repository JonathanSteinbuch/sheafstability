/*
 * Job.h
 *
 *  Created on: May 13, 2020
 *      Author: jonathan
 */
#include <gmpxx.h>
#include <iostream>
#include <boost/array.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <fstream>
#include <iterator>
#include <vector>

#include "../helpers/Modulus.hpp"

#include "../ring/PolyRing.hpp"
#include "../matrices/DensePolyMatrix.hpp"
#include "../matrices/SparsePolyMatrix.hpp"
#include "../ring/Poly.hpp"
#include "../helpers/Helpers.hpp"

#ifndef SRC_JOB_HPP_
#define SRC_JOB_HPP_

typedef Poly<mpz_class> Zpoly;
typedef PolyRing<mpz_class> ZRing;
typedef DensePolyMatrix<mpz_class> ZMatrix;
typedef SparsePolyMatrix<mpz_class> spZMatrix;

typedef Poly<numbermodulo> Fpoly;
typedef PolyRing<numbermodulo> FRing;
typedef DensePolyMatrix<numbermodulo> FMatrix;
typedef SparsePolyMatrix<numbermodulo> spFMatrix;

struct JobTypeDescriptor{
	string name;
	bool readMatrix;
	bool readPowers;
	bool readTwist;
};

const vector<JobTypeDescriptor> DescriptorList =
{{"semistability",true,false,false},
{"powers",true,true,false},
{"kernel",true,true,true},
};

struct JobOptionsType{
private:
	bool readuint(string line, string keyword, unsigned int &val);

	bool readstrings(string line, string keyword, vector<string> &values);

	bool readobjects(string line, string keyword, vector<string> &values);

	bool readmatrix(string line, string keyword, vector<vector<string> > &values);

	bool readmarker(string line, string keyword);

public:
	unsigned int type;

	unsigned int characteristic = 0;
	ofstream* output = NULL;
	ifstream* input = NULL;

	vector<string> variables;
	vector<string> strrelations;
	vector<vector<string> > strmatrix;

	unsigned int extPower = 1;
	unsigned int symPower = 1;
	unsigned int twist = 0;

	JobOptionsType(unsigned int type,ifstream* input, bool &retVal) ;
};

template <typename Scalar>
class JobType
{
public:
	virtual ~JobType(){};

	virtual bool doJob(PolyRing<Scalar>  &Ring, DensePolyMatrix<Scalar> &M, const JobOptionsType &jobOptions, const program_options &opt = noOptions) {return false;};
};

#endif /* SRC_JOB_HPP_ */
