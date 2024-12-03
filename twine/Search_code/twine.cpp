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
#include <fstream>

#include "SearchKey.hpp"
//#include "configBC.hpp"

using namespace std;
typedef unsigned int uint;

int main(int argc, char *argv[])
{

    int nbRounds = 8;



    fstream wd;
    wd.open("nosol.txt", ios::out);
    //vector<vector<uint8_t>> input(64, vector<uint8_t>(64, 1));
    //vector<vector<uint8_t>> output(64, vector<uint8_t>(64, 0));
    //for(uint i = 0; i < 64; i++)
    //{
    //        input[i][i] = 0;
    //        output[i][i] = 1;
    //}

    //#pragma omp parallel for num_threads(16)
    for(uint i = 6; i < 7; i++)
    {
    	   //#pragma omp parallel for num_threads(8)
	   for(uint j = 7; j < 8; j+=8)
	   {
    		   fstream fd;
		   fd.open("./"+ to_string(nbRounds) +"result/k"+to_string(i)+"__"+to_string(j), ios::out);
		   auto solution = MiddleSearch_backtrack(nbRounds, i, j);
		   cout << get<2>(solution).size() << endl;
		   if(get<2>(solution).size() > 0)
		   {
			fd << i  << "    " << j << "    has solution" << endl;
			cout << i  << "    " << j << "    has solution" << endl;
			for(uint index = 0; index < uint(nbRounds-1); index++)
			{
				for(uint index1 = 0; index1 < 64; index1++)
				{
					fd << int(get<2>(solution)[index][index1]);
					cout << int(get<2>(solution)[index][index1]);
				}
				fd << endl;
				cout << endl;
			}
		    }
		    else
		    {
			fd << i  << "    " << j << "    no solution!" << endl;
			wd << i  << "    " << j << "    no solution!" << endl;
			cout << i << "    " << j << "   no solution!" << endl;
		    }

		    fd << endl;
		    fd.close();
	   }
    } 
    wd.close();
    return 0;
}
