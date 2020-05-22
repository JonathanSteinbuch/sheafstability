/*
 * Job.cpp
 *
 *  Created on: May 13, 2020
 *      Author: jonathan
 */

#include "JobType.hpp"

bool JobOptionsType::readuint(string line, string keyword, unsigned int &val)
{
	std::size_t pos = line.find(keyword);
	if (pos != std::string::npos)
	{
		string num = line.substr(pos + keyword.length());

		val = std::stoi(num);
		return true;
	} else{
		return false;
	}
}

bool JobOptionsType::readstrings(string line, string keyword, vector<string> &values)
{
	std::size_t pos = line.find(keyword);
	bool ret = (pos != std::string::npos);
	while (pos != std::string::npos)
	{
		std::size_t p1 = line.find("\"", pos);
		std::size_t p2 = line.find("\"", p1 + 1);
		if (p1 != std::string::npos && p2 != std::string::npos)
		{
			string var = line.substr(p1 + 1, p2 - (p1 + 1));
			values.push_back(var);
		}
		pos = line.find(",", p2);
	}
	return ret;
}

bool JobOptionsType::readobjects(string line, string keyword, vector<string> &values)
{
	std::size_t pos = line.find(keyword);
	bool ret = (pos != std::string::npos);
	std::size_t p1 = pos + keyword.length();
	while (pos != std::string::npos)
	{

		std::size_t p2 = line.find(",", p1);
		if (p2 == std::string::npos)
		{
			pos = p2;
			p2 = line.length();
		}

		string rel = line.substr(p1, p2 - p1);
		values.push_back(rel);

		p1 = p2 + 1;
	}
	return ret;
}

bool JobOptionsType::readmatrix(string line, string keyword, vector<vector<string> > &values)
{
	std::size_t pos = line.find(keyword);
	bool ret = (pos != std::string::npos);

	std::size_t y1 = line.find("{", pos);
	while (pos != std::string::npos && y1 != std::string::npos)
	{
		vector<string> row;
		std::size_t x1 = line.find("{", y1 + 1);
		if (x1 == std::string::npos)
		{
			break;
		}
		std::size_t x2p = line.find("}", x1 + 1);
		while (x1 != std::string::npos)
		{
			std::size_t x2 = line.find(",", x1 + 1);
			if (x2 != std::string::npos)
			{
				x2 = min(x2, x2p);
			}
			else
			{
				x2 = x2p;
			}
			string entry = line.substr(x1 + 1, x2 - (x1 + 1));
			row.push_back(entry);
			if (x2 == x2p)
			{
				x1 = std::string::npos;
			}
			else
			{
				x1 = x2;
			}
		}
		values.push_back(row);
		y1 = x2p + 1;
	}
	return ret;
}

bool JobOptionsType::readmarker(string line, string keyword){
	std::size_t pos = line.find(keyword);
	if (pos != std::string::npos)
	{
		return true;
	} else{
		return false;
	}
}

const string keywords[] = {"characteristic:", "variables:", "relations:", "matrix:", "exteriorpower:", "spower:" , "twist:"};

JobOptionsType::JobOptionsType(unsigned int type, ifstream* input, bool &retVal) :type(type),input(input)
{
		unsigned int readsize = 3+(DescriptorList[type].readMatrix)+2*(DescriptorList[type].readPowers)+(DescriptorList[type].readTwist);

		vector<int> hasRead(readsize,0);

		string line;
		while (getline(*(input), line))
		{

			//read characteristic
			hasRead[0] |= readuint(line,keywords[0],characteristic);

			//read variables
			hasRead[1] |= readstrings(line,keywords[1],variables);

			//read relations
			hasRead[2] |= readobjects(line,keywords[2],strrelations);

			if(DescriptorList[type].readMatrix)
			{
				//read matrix
				hasRead[3] |= readmatrix(line,keywords[3],strmatrix);
			}
			if(DescriptorList[type].readPowers)
			{
				hasRead[4] |= readuint(line, keywords[4], extPower);

				hasRead[5] |= readuint(line, keywords[5], symPower);
			}

			if(DescriptorList[type].readTwist)
			{
				hasRead[6] |= readuint(line, keywords[6], twist);
			}
		}
		for(unsigned int i = 0; i < readsize; i++)
		{
			if(hasRead[i] != 1)
			{
				cerr << "There has to be exactly one line with " << keywords[i] << "in the input file.";
				retVal = false;
				return;
			}
		}
		retVal = true;

}
