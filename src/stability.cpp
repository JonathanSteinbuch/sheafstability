//============================================================================
// Name        : stability.cpp
// Author      : Jonathan Steinbuch
//============================================================================

#include "jobs/ReaderAndCaller.hpp"

using namespace std;

const float version1 = 0;
const float version2 = 1;
const float version3 = 0;



int main(int ac, char* av[])
{
//	jacobiTaylorTable();
//	return 0;

	try
	{
		program_options opt;
		po::options_description generic("Generic options");
		generic.add_options()("version,v", "print the version string")("help", "produce a help message")("config,c", po::value<string>(&opt.config_file)->default_value("sheaf_stability.cfg"), "name of a configuration file");

		po::options_description config("Configuration");
		config.add_options()
				("verbosity", po::value<int>(&opt.verbosity)->default_value(1), "Verbosity Level 0-3, Default: 1, O is silent unless there is an error")
				("input-file", po::value<string>(&opt.input_file), "input file")
				("output-file", po::value<string>(&opt.output_file)->default_value("sheaf_out.txt"), "output file")
				("exterior-powers,e", po::bool_switch(&opt.exteriorPowers)->default_value(false), "Using exterior powers in the algorithm")
				("stop-unstable,s", po::bool_switch(&opt.stopIfUnstable)->default_value(false), "Stop immediately if a destabilizing section has been found")
				("linear-progression,l", po::bool_switch(&opt.linearPowerProgresison)->default_value(false), "Linear progression of powers instead of fast lookahead")
		 	 	("pre-analyze,a", po::bool_switch(&opt.preAnalyze)->default_value(false), "Use symbolical analysis step before the Gauß algorithm")
				("forceQ,q", po::bool_switch(&opt.forceQ)->default_value(false), "Force Computations over Q (in characteristic 0). By Default Z is used where possible")
				("output-for-M2,m", po::bool_switch(&opt.outputM2)->default_value(false), "Write output matrices in a format easily readable by Macaulay2")
				("output-for-latex,x", po::bool_switch(&opt.outputLatex)->default_value(false), "Write output matrices in a format easily usable in Latex")
				("compute-full-kernel,f", po::bool_switch(&opt.computeFullKernel)->default_value(false), "Always compute the full kernels and not just whether the kernel is nonzero")
				("long-form-polys,p", po::bool_switch(&opt.longFormPolynomials)->default_value(false), "Output Polynomials in a less condensed format");

		po::options_description cmdline_options;
		cmdline_options.add(generic).add(config);

		po::options_description config_file_options;
		config_file_options.add(config);

		po::options_description visible("Allowed options");
		visible.add(generic).add(config);

		po::positional_options_description p;
		p.add("input-file", 1);
		p.add("output-file", 1);

		po::variables_map vm;
		store(po::command_line_parser(ac, av).options(cmdline_options).positional(p).run(), vm);
		notify(vm);

		ifstream ifs(opt.config_file.c_str());
		if (!ifs)
		{
			cout << "can not open config file: " << opt.config_file << "\n";
			return 0;
		}
		else
		{
			store(parse_config_file(ifs, config_file_options), vm);
			notify(vm);
		}

		if (vm.count("help"))
		{
			cout << visible << "\n";
			return 0;
		}

		if (vm.count("version"))
		{
			cout << "Sheaf stability computer, version " << version1 << "." << version2 << "." << version3 << "\n";
			return 0;
		}

		ReaderAndCaller thisRaC(opt);
		return 0;
	} catch (exception& e)
	{
		cout << e.what() << "\n";
		return 1;
	}
	return 0;
}

