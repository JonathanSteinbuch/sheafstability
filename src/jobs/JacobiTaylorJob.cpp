/*
 * JacobiTaylorJob.cpp
 *
 *  Created on: May 13, 2020
 *      Author: jonathan
 */

#include "JacobiTaylorJob.hpp"



template <typename Scalar>
void jacobiTaylor()
{

	const unsigned int vars = 3;

	//monomBasis<vars,maxdegree,Scalar> mB({"x","y","z","w"},degrevlex);
	//monomBasis<Scalar> mB(vars, {"t","y","z","w","x"},degrevlex,4);
	//monomBasis<Scalar> mB(vars,{ "x", "y", "w", "z" }, degrevlex, -1,{ 3, 3, 3, 2 });
	PolyRing<Scalar> mB(vars,
	{ "x", "y", "z" }, degrevlex);

	//poly t = Zpoly(&mB,"t");
	Zpoly x = Zpoly(&mB, "x");
	//poly w = poly(&mB, "w");
	Zpoly y = Zpoly(&mB, "y");
	Zpoly z = Zpoly(&mB, "z");

	//poly relation = (x^3)+(y^3)+(z^3)+(w^3);
	//poly relation = x*y-(z^2);
	//poly relation = (x^2)+(y^2)+(z^2);
	//poly relation = (x^2)*z+(y^3);
	//vector<poly>relations = {poly(&mB,"x8x*x+2y^4(z+x)")};
//	poly relation = (y ^ 2) + (z ^ 2) + (w ^ 3);
	Zpoly relation = (x ^ 3) + (y ^ 3) + (z ^ 3);
	//relation.homogenize(0);
	vector<Zpoly> relations =
	{ relation };

	cout << relation << endl << endl;
	//poly relation = x*y;
	PolyRing<Scalar> S(&mB, relations);

	/*monomBasis<mpq_class> kX(1,
	 { "n" }, degrevlex);
	 mBpolynomial<mpq_class> hP(&kX);
	 mB.computeHilbertPoly(hP, mB.gB);*/
	S.getHilbert();

	//cout << "The Hilbert polynomial of X is P(n)=" << hP << endl;

	/*		int deg = -1;
	 for(unsigned int i=0; i < mB.list.size(); i++)
	 {
	 if(deg < (int)mB.getDegree(i))
	 {
	 deg = mB.getDegree(i);
	 cout << "Degree: " << deg << endl;
	 }
	 mB.print(cout,i);
	 cout << " = " << mB.lookup[i] << endl;
	 }
	 */

	unsigned int rank = 0;

	for (unsigned int i = 1; i < 30; i++)
	{

		cout << "Jacobi-Taylor Matrix for n=" << i << endl;

		ZMatrix JT = ZMatrix::jacobiTaylor(&S, relation, i);

		cout << JT << endl;

		clock_t t = clock();

		for (int d = (int) relation.degree(); d < (int) relation.degree() + 3; d++)
		{
			ZMatrix kernel = JT.kernel(d);

			cout << "Kernel(" << d - (int) relation.degree() << "):" << endl;
			//cout << kernel << endl;
			kernel.printAsOperator(cout);
			if (kernel.getCols() > 0)
			{

				vector<ZMatrix> list;
				for (unsigned int i = 0; i < kernel.getRows(); i++)
				{
					list.push_back(ZMatrix::irrelevantPower(S, d - kernel.getOutdegrees()[i]));
				}
				ZMatrix Mn = directsum(list);

				//cout << Mn << endl;

				//Mn.computeDegrees(kernel.getOutdegrees());

				ZMatrix kernelmod = kernel % Mn;

				//cout << kernelmod << endl;
				kernelmod.printAsOperator(cout);
				if (d == (int) relation.degree() - (int) i)
				{
					rank += kernelmod.getCols();
				}
			}
			//rank += kernel.getCols();

			//cout << JT*kernel << endl;

		}

		cout << "Rank: " << rank + 1 << ", cum. hilbert function:" << S.getIndexOfBasisDegStart(i + 1) << endl;

		cout << "Norm: " << rank + 1 << "/" << S.getMonN(i + 1) << " = " << (float) (rank + 1) / (float) (S.getMonN(i + 1)) << endl;

		t = clock() - t;
		cout << "Time spent: " << t / CLOCKS_PER_SEC << "s (" << t << " internal ticks)" << endl << endl;
	}
}

template <typename Scalar>
void jacobiTaylorTable()
{

	const unsigned int vars = 3;

	//monomBasis<vars,maxdegree,Scalar> mB({"x","y","z","w"},degrevlex);
	//monomBasis<Scalar> mB(vars, {"t","y","z","w","x"},degrevlex,4);
	//monomBasis<Scalar> mB(vars,{ "x", "y", "w", "z" }, degrevlex, -1,{ 3, 3, 3, 2 });
	PolyRing<Scalar> mB(vars,
	{ "x", "y", "z" }, degrevlex);

	//poly t = Zpoly(&mB,"t");
	Zpoly x = Zpoly(&mB, "x");
	//poly w = poly(&mB, "w");
	Zpoly y = Zpoly(&mB, "y");
	Zpoly z = Zpoly(&mB, "z");

	//poly relation = (x^3)+(y^3)+(z^3)+(w^3);
	//poly relation = x*y-(z^2);
	//poly relation = (x^2)+(y^2)+(z^2);
	//poly relation = (x^2)*z+(y^3);
	//vector<poly>relations = {poly(&mB,"x8x*x+2y^4(z+x)")};
//	poly relation = (y ^ 2) + (z ^ 2) + (w ^ 3);
	Zpoly relation = (x ^ 4) + (y ^ 4) + (z ^ 4);
	//poly relation = (x ^ 7) + (y ^ 7) + (z ^ 7);
	//relation.homogenize(0);
	vector<Zpoly> relations =
	{ relation };

	cout << relation << endl << endl;
	//poly relation = x*y;
	PolyRing<Scalar> S(&mB, relations);

	/*monomBasis<mpq_class> kX(1,
	 { "n" }, degrevlex);
	 mBpolynomial<mpq_class> hP(&kX);
	 mB.computeHilbertPoly(hP, mB.getGB());
	 //mB.getHilbert();

	 cout << "The Hilbert polynomial of X is P(n)=" << hP << endl;*/

	/*		int deg = -1;
	 for(unsigned int i=0; i < mB.list.size(); i++)
	 {
	 if(deg < (int)mB.getDegree(i))
	 {
	 deg = mB.getDegree(i);
	 cout << "Degree: " << deg << endl;
	 }
	 mB.print(cout,i);
	 cout << " = " << mB.lookup[i] << endl;
	 }
	 */

	vector<vector<int>> results;

	int maxd = 32;
	unsigned int maxn = 10;

	cout << "n\\d";
	for (int d = relation.degree(); d <= maxd; d++)
	{
		cout << "\t" << d - relation.degree();
	}
	cout << endl;

	for (unsigned int i = 1; i <= maxn; i++)
	{

		cout << i;

		ZMatrix JT = ZMatrix::jacobiTaylor(&S, relation, i);

		vector<int> row;

		for (int d = relation.degree(); d <= maxd; d++)
		{
			ZMatrix kernel = JT.kernel(d);

			row.push_back(kernel.getCols());

			int col = d - relation.degree();

			if (i == 1)
			{
				cout << "\t" << row[col];
			}
			else
			{
				cout << "\t" << row[col] - results[i - 2][col];
			}

		}

		cout << endl;

		results.push_back(row);
	}
}

template<typename Scalar>
bool JacobiTaylorJob<Scalar>::doJob(PolyRing<Scalar> &Ring, DensePolyMatrix<Scalar> &M,const JobOptionsType &jobOptions, const program_options &opt)
{
	return true;
}


template class JacobiTaylorJob<mpz_class>;
template class JacobiTaylorJob<mpq_class>;
template class JacobiTaylorJob<numbermodulo>;
