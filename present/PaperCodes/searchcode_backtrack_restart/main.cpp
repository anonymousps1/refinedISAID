#include <iostream>
#include <cstdint>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <algorithm>
#include <utility>
#include <stdlib.h>     /* atoi */

#include "BCData.hpp"
#include "configBC.hpp"
#include<ctime>

using namespace std;
clock_t start, end1;
typedef unsigned int uint;


string permuteNibbleSkinny(string const & s){
	string tmp(s);
	vector<uint> SR({0,5,10,15,4,9,14,3,8,13,2,7,12,1,6,11});
	for(uint i = 0; i < 16; i++)
		tmp[SR[i]] = s[i];
	return tmp;
}


int main(int argc, char *argv[])
{
	
	
	uint rMax = 12;

	int input, output, mode;


	for(int i = 0; i < argc; i++)
	{
		if(!strcmp("-i", argv[i])) input = atoi(argv[i+1]); 
		if(!strcmp("-o", argv[i])) output = atoi(argv[i+1]); 
		if(!strcmp("-m", argv[i])) mode = atoi(argv[i+1]); 
	}


	auto BCD = genDataSkinny64(rMax);


	

	if(mode == 1)
	{
		start = clock();

		auto iok = BCD.MiddleSearch_backtrack(output, input);

		end1 = clock();

		double endtime = (double)(end1 - start)/CLOCKS_PER_SEC;

		cout << "compute  time:       " << endtime << endl;

		//If we were able to do so, we get non-empty vectors
		if(get<0>(iok).size() > 0){
			cout << "Input : "; BCD.printParitySetHex(get<0>(iok)); cout << endl;
			cout << "Keys : " << endl;
			auto const & allKeys = get<2>(iok);
			for(auto const & key : allKeys){
				BCD.printParitySetHex(key); cout << endl;
			}
		}
		else
		{
			cout << "can not find a solution for input: " << input << " output: " << output << endl;
		}
	}
	else
	{
		start = clock();
		auto iok = BCD.MiddleSearch(output, input);

		end1 = clock();

		double endtime = (double)(end1 - start)/CLOCKS_PER_SEC;

		cout << "compute  time:       " << endtime << endl;

		//If we were able to do so, we get non-empty vectors
		if(get<0>(iok).size() > 0){
			cout << "Input : "; BCD.printParitySetHex(get<0>(iok)); cout << endl;
			cout << "Keys : " << endl;
			auto const & allKeys = get<2>(iok);
			for(auto const & key : allKeys){
				BCD.printParitySetHex(key); cout << endl;
			}
		}
		else
		{
			cout << "can not find a solution for input: " << input << " output: " << output << endl;
		}
	}

	return 0;

}
