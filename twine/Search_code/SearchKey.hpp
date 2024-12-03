#ifndef H_SEARCHKEY
#define H_SEARCHKEY

#include <cstdint>
#include <vector>
#include <set>
#include <stack>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <bitset>
#include <random>
#include <tuple>
#include <map>

#include <boost/dynamic_bitset.hpp>
//#include <gmpxx.h>

#include "gurobi_c++.h"
#include "Callback.hpp"






uint64_t countSSBTrails(std::vector<uint8_t> & u,
			std::vector<uint8_t> & v);



std::pair<uint64_t,bool> countTrailsFullCipher(std::vector<uint8_t> const & input,
						std::vector<uint8_t> const & output, 
						std::vector<std::vector<uint8_t>> const & keyval);


std::pair<std::vector<std::vector<std::vector<uint8_t>>>, std::vector<std::vector<boost::dynamic_bitset<>>>>
searchkey(int const nbRounds,
						 std::vector<uint8_t> const & input,
						 std::vector<uint8_t> const & output,
						 std::vector<std::vector<boost::dynamic_bitset<>>> const & allTriedKeyPattern = std::vector<std::vector<boost::dynamic_bitset<>>>());


std::tuple<std::vector<uint8_t>, std::vector<uint8_t>, std::vector<std::vector<uint8_t>>>
MiddleSearch_backtrack(uint const rounds, uint const indexInput, uint const indexOutput);


//std::pair<std::vector<std::vector<std::vector<uint8_t>>>, 
//          std::vector<std::vector<boost::dynamic_bitset<>>>>
std::vector<std::vector<std::vector<uint8_t>>>
searchKeyPattern(int const nbRounds,
                 std::vector<uint8_t> const & input,
                 std::vector<uint8_t> const & output);



std::tuple<std::vector<uint8_t>, std::vector<uint8_t>, std::vector<std::vector<uint8_t>>>
imSearch(uint const rMax, std::vector<uint8_t>  const & u, std::vector<uint8_t> const & v); 


std::vector<std::vector<uint8_t>>
searchKeyall(uint rMax, uint rMiddle, std::vector<uint8_t> & u, std::vector<uint8_t> & v, std::vector<std::vector<uint8_t>> & valk, std::vector<std::vector<uint8_t>> & valx); 
//std::tuple<std::vector<std::vector<uint8_t>>, std::vector<std::vector<uint8_t>>, std::vector<std::vector<uint8_t>>>
std::set<std::pair<std::vector<uint8_t>, std::vector<uint8_t>>>
backsearch(uint rMax, uint rMiddle, std::vector<uint8_t> & u, std::vector<uint8_t> & v, std::vector<std::vector<uint8_t>> & valk, std::vector<std::vector<uint8_t>> & valx); 


std::pair<uint64_t, int> counttrail_middle(uint rounds, std::vector<uint8_t> & u, std::vector<uint8_t> & v, std::vector<std::vector<uint8_t>> & k);

std::pair<uint64_t, bool> counttrail(std::vector<uint8_t> & u, std::vector<uint8_t> & v, std::vector<uint8_t> & k);

std::pair<uint64_t, bool> counttrail_last(std::vector<uint8_t> & u, std::vector<uint8_t> & v, std::vector<uint8_t> & k);

std::pair<std::vector<std::vector<uint8_t>>, std::vector<std::vector<uint8_t>>> countfirst(std::vector<uint8_t> & u);
std::pair<std::vector<std::vector<uint8_t>>, std::vector<std::vector<uint8_t>>> countlast(std::vector<uint8_t> & u);




void printParitySetHex(std::vector<uint8_t> const & x);

void printBits(uint64_t const x, uint nbBits);
//print the first #nbBits of x in cout

void printVec(std::vector<uint8_t> const & x);
//print the elements of x, wihtout separators

std::string hex_char_to_bin(char c);
//hex to bin conversion, LSB on the left in the binary representation

std::vector<uint8_t> hexToVec(std::string const & s);
/*
	From an hex string get the corresponding vector of bits
	LSB of each hex char is put on the left
	e.g. 
	s = "1"  returns {1,0,0,0}
	s = "2B" returns {0,1,0,0,1,1,0,1}
*/




#endif
