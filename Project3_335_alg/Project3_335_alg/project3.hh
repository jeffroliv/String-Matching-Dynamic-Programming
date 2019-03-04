///////////////////////////////////////////////////////////////////////////////
// maxprotein.hh
//
// Compute the set of foods that maximizes protein, within a calorie budget,
// with the greedy method or exhaustive search.
//
///////////////////////////////////////////////////////////////////////////////


#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <queue>
#include <sstream>
#include <string>
#include <vector>
#include <cstdint>

using namespace std;

// Simple structure for a single protein
struct Protein {
	Protein() {
		description = "";
		sequence = "";
	}
	Protein(std::string desc, std::string seq) {
		description = desc;
		sequence = seq;
	}
	std::string		description;
	std::string 	sequence;
};

// Alias for a vector of shared pointers to Protein objects.
typedef std::vector<std::shared_ptr<Protein>> ProteinVector;


// -------------------------------------------------------------------------
// Load all the proteins from a standard FASTA format file with one line
// per sequence (multi-line sequences are not allowed).
// Returns false on I/O error.
// -------------------------------------------------------------------------
bool load_proteins(ProteinVector & proteins, const std::string& path)
{
	//std::cout << "Loading proteins from [" << path << "]" << std::endl;
	proteins.clear();
	std::ifstream ifs(path.c_str());
	if (!ifs.is_open() || !ifs.good()) {
		std::cout << "Failed to open [" << path << "]" << std::endl;
		return false;
	}
	int proteinsLoaded = 0;
	bool have_description = false;
	std::shared_ptr<Protein> newProtein = nullptr;
	while (!ifs.eof()) {
		std::string lineBuffer;
		std::getline(ifs, lineBuffer);
		if (ifs.eof()) {
			break;
		}
		if (lineBuffer.size() == 0) {
			continue;
		}
		if (lineBuffer[0] == '>') {
			newProtein = std::shared_ptr<Protein>(new Protein);
			newProtein->description = lineBuffer.substr(1);
			have_description = true;
		}
		else if (have_description) {
			newProtein->sequence = lineBuffer;
			proteins.push_back(newProtein);
			proteinsLoaded++;
			have_description = false;
		}
	}

	ifs.close();
	//std::cout << "Loaded " << proteinsLoaded << " proteins from [" << path << "]" << std::endl;

	return true;
}

//Find_max function takes in three numbers as arguments and outputs the max of the three numbers
//--------------------------------------------------------------------------
int find_max(int num1, int num2, int num3)
{
	int max = num1;
	if (num2 > max)
	{
		max = num2;
	}
	if (num3 > max)
	{
		max = num3;
	}

	return max;
}

// -------------------------------------------------------------------------
int dynamicprogramming_longest_common_subsequence(const std::string & string1,
	const std::string & string2)
{
	const int n = string1.size();
	const int m = string2.size();

	int max_len = 0;

	if (n == 0 && m == 0)
	{
		return 0;
	}
	if (n == 1 && m == 1)
	{
		return 1;
	}


	int up, left, diag;

	int D[n + 1][m + 1];

	for (int i = 0; i < n; i++)
	{
		D[i][0] = 0;
	}
	for (int j = 0; j < m; j++)
	{
		D[0][j] = 0;
	}

	for (int i = 1; i < n; i++)
	{
		for (int j = 1; j < m; j++)
		{
			up = D[i - 1][j];
			left = D[i][j - 1];
			diag = D[i - 1][j - 1];
			if (string1[i - 1] == string2[j - 1])
			{
				diag = diag + 1;
			}
			D[i][j] = find_max(up, left, diag);
			max_len = D[i][j] + 1;
		}
	}

	return max_len;
}

// -------------------------------------------------------------------------
std::unique_ptr<std::vector<std::string>> generate_all_subsequences(const std::string & sequence)
{
	std::unique_ptr<std::vector<std::string>> R(new std::vector<std::string>);
	int n = pow(2, sequence.size());
	for (uint64_t bits = 0; bits < n - 1; bits++)
	{
		std::string subsequence = "";
		for (int j = 0; j < sequence.size() - 1; j++)
		{
			if (((bits >> j) & 1) == 1)
			{
				subsequence += subsequence[j];
			}
		}
		R->push_back(subsequence);
	}

	return R;
}


// -------------------------------------------------------------------------
int exhaustive_longest_common_subsequence(const std::string & string1,
	const std::string & string2)
{
	if (string1.empty() || string2.empty())
	{
		return 0;
	}
	std::unique_ptr<std::vector<std::string>> all_subseq1(new std::vector<std::string>);
	std::unique_ptr<std::vector<std::string>> all_subseq2(new std::vector<std::string>);
	int best_score = 0;
	all_subseq1 = generate_all_subsequences(string1);
	all_subseq2 = generate_all_subsequences(string2);
	for (auto s1 : *all_subseq1)
	{
		for (auto s2 : *all_subseq2)
		{
			if (s1 == s2 && s1.size() > best_score)
			{
				best_score = s1.size() + 1;
			}
		}
	}
	return best_score;
}


// -------------------------------------------------------------------------
std::shared_ptr<Protein> exhaustive_best_match(ProteinVector & proteins, const std::string & string1)
{
	int best_i = 0;
	int best_score = 0;

	for (int i = 0; i < proteins.size(); i++)
	{
		int score = exhaustive_longest_common_subsequence(proteins.at(i).get()->sequence, string1);
		if (score > best_score)
		{
			best_score = score;
			best_i = i;
		}
	}

	return proteins[best_i];
}

// -------------------------------------------------------------------------
std::shared_ptr<Protein> dynamicprogramming_best_match(ProteinVector & proteins, const std::string & string1)
{
	int best_i = 0;
	int best_score = 0;

	for (int i = 0; i < proteins.size(); i++)
	{
		int score = dynamicprogramming_longest_common_subsequence(proteins.at(i).get()->sequence, string1);
		if (score > best_score)
		{
			best_score = score;
			best_i = i;
		}
	}

	return proteins[best_i];
}