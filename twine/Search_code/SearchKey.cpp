#include "SearchKey.hpp"

#define MAX_GRBMIPGapAbs 8

#define MAX_NUMBER_SEARCHKEY 500

#define TIME_LIMIT_SEARCHKEY 300

#define MAX_GRB_THREAD 8

#define ENFORCE_NONZERO_LASTKEY true 
#define ENFORCE_NONZERO_FIRSTKEY true 

#define WEIGHT_THRESHOLD_STEP1 40 

#define TIME_LIMIT_ONESTEPDYNAMIC 30			/*Time limit (in seconds) in improvedDynamicSearch to find one iteration in the first step*/
#define TIME_LIMIT_ONESTEPDYNAMIC 30			/*Time limit (in seconds) in improvedDynamicSearch to find one iteration in the first step*/
#define DYNAMIC_STEP1_MIPFOCUS 1 
#define AGRESSIVE_STRATEGY_MINDEGREE false 

#define MAX_NUMBER_SOLCOUNT_FULLCIPHER 100000

#define TIME_LIMIT_COUNTTRAILS_FULLCIPER 60


using namespace std;
using namespace boost;

typedef unsigned int uint;
typedef pair<uint32_t,uint32_t> pairuint32;
typedef tuple<uint32_t, uint32_t, uint64_t> tup3uint;

std::random_device rd; 
std:: mt19937 prng(rd()); 
std::uniform_int_distribution<uint> randintSeedGurobi(0, 1999999999);
std::uniform_int_distribution<uint> randintCilumn;


uint64_t countSSBTrails(vector<uint8_t> & u, vector<uint8_t> &v)
{
	GRBEnv gurobiEnv;
	gurobiEnv.set(GRB_IntParam_OutputFlag, 0);
	gurobiEnv.set(GRB_IntParam_LogToConsole, 0);

	GRBModel m(gurobiEnv, "./modelSSB/SSB.mps");
	m.set(GRB_IntParam_Threads, 1);

	vector<GRBVar> xvar(16);
	vector<GRBVar> zvar(16);

	for(uint i = 0; i < 16; i++)
	{
		xvar[i] = m.getVarByName("x_" + to_string(i));
		zvar[i] = m.getVarByName("z_" + to_string(i));
	}

	for(uint i = 0; i < 16; i++)
	{
		m.addConstr(xvar[i] == u[i]);
		m.addConstr(zvar[i] == v[i]);
	}

	m.set(GRB_IntParam_PoolSearchMode, 2);
	m.set(GRB_IntParam_PoolSolutions, 20000000);

	m.optimize();

	return m.get(GRB_IntAttr_SolCount);

}



pair<uint64_t, bool> countTrailsFullCipher(vector<uint8_t> const & input, vector<uint8_t> const & output, vector<vector<uint8_t>> const & keyval)
{
	GRBEnv gurobiEnv;
	gurobiEnv.set(GRB_IntParam_OutputFlag, 0);

	int round = keyval.size() + 1;

	GRBModel m(gurobiEnv, "./modeltwine_SSB/"+ to_string(round)+"_twine.mps");

	for(uint i = 0; i < 64; i++)
	{
		m.addConstr(m.getVarByName("x_0_"+ to_string(i)) == input[i]);
		m.addConstr(m.getVarByName("v_"+ to_string(round -1) +"_"+ to_string(i)) == output[i]);
	}
	
	for(uint r  = 0; r< keyval.size(); r++)
	{
		for(uint i = 0; i < 64; i++)
		{
			m.addConstr(m.getVarByName("k_" + to_string(r) + "_" +to_string(i)) == keyval[r][i]);
		}
	}

	//Count trails
	m.set(GRB_IntParam_PoolSearchMode, 2);
	m.set(GRB_IntParam_PoolSolutions, 2000000000); //maximum number of solution allowed by Gurobi
	auto seedGurobi = randintSeedGurobi(prng);
	m.set(GRB_IntParam_Seed, seedGurobi);
	//Dummy objective
	GRBLinExpr objExpr;
	for(uint i = 0; i < 64; i++){
		objExpr += m.getVarByName("x_0_"+to_string(i));
	}
	m.setObjective(objExpr);

	cout << "Counting the total number of trails" << endl;
	// callbackCount cb = callbackCount(MAX_NUMBER_SOLCOUNT_FULLCIPHER);
	m.set(GRB_DoubleParam_TimeLimit, 600);
	m.set(GRB_IntParam_Threads, 8);
	// m.setCallback(&cb);	
	m.set(GRB_IntParam_SolutionLimit, 200000);
	m.update();
	m.optimize();

	cout << "status code: " << m.get(GRB_IntAttr_Status) << endl;

	if(m.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
		return make_pair(0, true);

	uint totalNumberTrails = m.get(GRB_IntAttr_SolCount);

	bool fullOptimized = true;
	if(m.get(GRB_IntAttr_Status) == GRB_OPTIMAL){
		cout << endl << totalNumberTrails << " total trails" << endl;
	}
	else if(m.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT){
		cout << endl << "counting solutions aborted, too much time spent" << endl;
		cout << totalNumberTrails << " trails found" << endl;
		fullOptimized = false;
	}
	else if(m.get(GRB_IntAttr_Status) == GRB_SOLUTION_LIMIT){
		cout << endl << "counting solutions aborted, too many solutions" << endl;
		cout << totalNumberTrails << " trails found" << endl;
		fullOptimized = false;
	}
	else{
		cout << endl << "counting solutions aborted with Status " << m.get(GRB_IntAttr_Status) << endl;
		cout << totalNumberTrails << " trails found" << endl;
		fullOptimized = false;
	}

	return make_pair(totalNumberTrails,fullOptimized);
	
}


pair<vector<vector<vector<uint8_t>>>, vector<vector<dynamic_bitset<>>>>
searchkey(int const nbRounds,
						   vector<uint8_t> const & input,
						   vector<uint8_t> const & output,
						   vector<vector<dynamic_bitset<>>> const & allTriedKeyPattern){


	//Read the base model
	GRBEnv gurobiEnv;
	gurobiEnv.set(GRB_IntParam_OutputFlag, 0);
	GRBModel m(gurobiEnv, "./modeltwine_SSB/"+to_string(nbRounds)+"_tw.mps");
	m.set(GRB_IntParam_LazyConstraints, 1);


	//Read variables
	//input/output variables
	uint blockSize = 64;
	vector<GRBVar> inputVar(blockSize);
	vector<GRBVar> outputVar(blockSize);

	for(uint i = 0; i < blockSize; i++){
		inputVar[i] = m.getVarByName("x_0_"+to_string(i));
		outputVar[i] = m.getVarByName("x_"+to_string(nbRounds)+"_"+to_string(i));
	}


	//All key variables (one key is not involved in an SSB)
	uint const nbroundspoffset = nbRounds;
	vector<vector<GRBVar>> allKeyVar(nbroundspoffset, vector<GRBVar>(blockSize));
	for(uint r = 0; r < nbroundspoffset; r++){
		auto & allKeyVar_r = allKeyVar[r];
		for(uint i = 0; i < blockSize; i++){
			allKeyVar_r[i] = m.getVarByName("k_"+to_string(r)+"_"+to_string(i));
		}
	}

	uint p_in[16] = {0, 1, 2, 3, 4, 5, 8, 9, 6, 7, 10, 11, 12, 13, 14, 15};
	uint p_out[16] = {0, 1, 4, 5, 12, 13, 6, 7, 8, 9, 2, 3, 10, 11, 14, 15};

	//Prepare the callback variables
	vector<vector<vector<GRBVar>>> inSSBVar(nbroundspoffset, 
											vector<vector<GRBVar>>(4,
											vector<GRBVar>(16)));
	vector<vector<vector<GRBVar>>> outSSBVar(nbroundspoffset, 
											 vector<vector<GRBVar>>(4,
											 vector<GRBVar>(16)));

	for(uint r = 0; r < nbroundspoffset; r++){
		auto & inSSBVar_r = inSSBVar[r];
		auto & outSSBVar_r = outSSBVar[r];
		for(uint i = 0; i < 4; i++){
			auto & inSSBVar_r_i = inSSBVar_r[i];
			auto & outSSBVar_r_i = outSSBVar_r[i];
			for(uint j = 0; j < 16; j++){
				uint bit = 4*(p_in[4*i+j/4]) + j%4;
				uint bit1 = 4*(p_out[4*i+j/4]) + j%4;
				inSSBVar_r_i[j] = m.getVarByName("x_"+to_string(r)+"_"+to_string(bit));
				outSSBVar_r_i[j] = m.getVarByName("v_"+to_string(r)+"_"+to_string(bit1));
			}
		}
	}


	//Add input/output constraints
	for(uint i = 0; i < blockSize; i++){
		m.addConstr(inputVar[i] == input[i]);
		m.addConstr(outputVar[i] == output[i]);
	}

	//Remove known key pattern
	uint const allKeyVar_size = allKeyVar.size();
	for(auto const & keypattern : allTriedKeyPattern){
		GRBLinExpr cutExpr(0);
		for(uint r = 0; r < allKeyVar_size; r++){
			auto & keypattern_r = keypattern[r];
			auto & allKeyVar_r = allKeyVar[r];
			for(uint i = 0; i < blockSize; i++){
				if(keypattern_r[i] == 0) cutExpr += allKeyVar_r[i];
				else cutExpr += (1 - allKeyVar_r[i]);
			}
		}
		m.addConstr(cutExpr >= 1);
	}


	//Dummy objective
	GRBLinExpr objExpr(0);
	for(uint r = 0; r < allKeyVar_size; r++){
		auto & allKeyVar_r = allKeyVar[r];
		for(uint i = 0; i < blockSize; i++){
			objExpr += allKeyVar_r[i];
		}
	}
	m.setObjective(objExpr, GRB_MAXIMIZE);

	//setup params for enumerating solutions
	m.set(GRB_IntParam_PoolSearchMode, 2);
	m.set(GRB_IntParam_PoolSolutions, 20000000);
	m.set(GRB_DoubleParam_TimeLimit, 300);
	//m.set(GRB_DoubleParam_MIPGapAbs,MAX_GRBMIPGapAbs);
	//auto seedGurobi = randintSeedGurobi(prng);
	//m.set(GRB_IntParam_Seed, seedGurobi);

	//Search for key patterns so that we have an odd number of trail from input to output
	//Initialize Callback
	callbacksearchkey cb(nbRounds,inputVar,outputVar,allKeyVar,inSSBVar,outSSBVar);
	m.setCallback(&cb);
	m.set(GRB_IntParam_Threads, 8);
	
	m.update();
	m.optimize();
	
	uint const nbSol = m.get(GRB_IntAttr_SolCount);
	if(nbSol > 0){
		cout << "Number of optimal solutions : " << nbSol << endl;
		//We found a key pattern
		vector<vector<vector<uint8_t>>> Allkeyval(nbSol, vector<vector<uint8_t>>(nbroundspoffset, vector<uint8_t>(blockSize)));
		for(uint indexSol = 0; indexSol < nbSol; indexSol++){
			m.set(GRB_IntParam_SolutionNumber, indexSol);
			auto & keyval = Allkeyval[indexSol];
			for(uint r = 0; r < nbroundspoffset; r++){
				auto & keyval_r = keyval[r];
				auto & allKeyVar_r = allKeyVar[r];
				for(uint i = 0; i < blockSize; i++)
					keyval_r[i] = uint8_t(round(allKeyVar_r[i].get(GRB_DoubleAttr_Xn)));
			}
		}
		return make_pair(Allkeyval, cb.allSolutions);
	}
	else{
		//No key pattern found
		return make_pair(vector<vector<vector<uint8_t>>>(), vector<vector<dynamic_bitset<>>>());
	}
}

tuple<vector<uint8_t>, vector<uint8_t>, vector<vector<uint8_t>>>
MiddleSearch_backtrack(uint const rounds, uint const indexInput, uint const indexOutput)
{

	GRBEnv gurobiEnv;
	gurobiEnv.set(GRB_IntParam_OutputFlag, 0);

	uint blockSize = 64;
	vector<uint8_t> output(blockSize, 0);
	output[indexOutput] = 1;
	vector<uint8_t> input(blockSize, 1);
	input[indexInput] = 0;

	vector<uint8_t> output_inter(blockSize, 0);
	output_inter[indexOutput - 1] = 1;

	int rMax = (int)(rounds);
	int rMax_minus_1 = (int)(rounds) - 1;
	int rMax_minus_2 = (int)(rounds) - 2;
	//Generate the model


	vector<vector<uint8_t>> valX(rMax, vector<uint8_t>(blockSize, 0));
	vector<vector<uint8_t>> valK(rMax_minus_1, vector<uint8_t>(blockSize, 0));


	int rMiddle = rMax_minus_2;

	GRBModel m(gurobiEnv,"./modeltwine_SSB/"+to_string(rounds)+"_twine.mps");

	m.set(GRB_IntParam_LazyConstraints, 1);


	//If needed, enforce the last key to be non zero
	//GRBLinExpr sumLastKey(0);
	//for(uint i = 0; i < blockSize; i++)
	//{
	//	sumLastKey += m.getVarByName("k_"+to_string(rMax_minus_2)+"_"+to_string(i));
	//	sumLastKey += m.getVarByName("k_"+to_string(rMax_minus_2-1)+"_"+to_string(i));
	//}
	//m.addConstr(sumLastKey >= 1);

	//GRBLinExpr sumFirstKey(0);
	//for(uint i = 0; i < blockSize; i++)
	//{
	//	sumFirstKey += m.getVarByName("k_0_"+to_string(i));
	//	sumFirstKey += m.getVarByName("k_1_"+to_string(i));
	//}
	//m.addConstr(sumFirstKey >= 1);


	//Fix the input
	for (uint i = 0; i < 64; i++)
		m.addConstr(m.getVarByName("x_0_"+to_string(i)) == input[i]);

	//Fix the output
	for(uint i = 0; i < blockSize; i++)
		m.addConstr(m.getVarByName("v_"+to_string(rMax_minus_1)+"_"+to_string(i)) == output[i]);



	//Fix the known keys
	for(int r = rMiddle+1; r < rMax_minus_1; r++){
		auto const & valK_r = valK[r];
		for(uint i = 0; i < blockSize; i++)
			m.addConstr(m.getVarByName("k_"+to_string(r)+"_"+to_string(i)) == valK_r[i]);
	}

	//Fix the known round input
	for(int r = rMiddle+1; r < rMax_minus_1; r++){
		auto const & valX_r = valX[r];
		for(uint i = 0; i < blockSize; i++)
			m.addConstr(m.getVarByName("x_"+to_string(r)+"_"+to_string(i)) == valX_r[i]);
	}

	//Objective variables
	vector<GRBVar> XmidVar(blockSize);
	vector<GRBVar> KmidVar(blockSize);
	for(uint i = 0; i < blockSize; i++){
		XmidVar[i] = m.getVarByName("x_"+to_string(rMiddle)+"_"+to_string(i));
		KmidVar[i] = m.getVarByName("k_"+to_string(rMiddle)+"_"+to_string(i));
	}

	//Set the objective
	GRBLinExpr objExpr(0);
	for(uint i = 0; i < blockSize; i++){
		objExpr	+= XmidVar[i];
		objExpr -= m.getVarByName("v_"+to_string(rMiddle+1)+"_"+to_string(i));
		objExpr += 0.01*KmidVar[i];
	}
	m.setObjective(objExpr, GRB_MAXIMIZE);

	//Seed
	auto seedGurobi = randintSeedGurobi(prng);
	m.set(GRB_IntParam_Seed, seedGurobi);

	uint p_in[16] = {0, 1, 2, 3, 4, 5, 8, 9, 6, 7, 10, 11, 12, 13, 14, 15};
	uint p_out[16] = {0, 1, 4, 5, 12, 13, 6, 7, 8, 9, 2, 3, 10, 11, 14, 15};

	//Prepare the callback variables
	vector<vector<vector<GRBVar>>> inSSBVar(rMax-rMiddle, 
											vector<vector<GRBVar>>(4,
											vector<GRBVar>(16)));
	vector<vector<vector<GRBVar>>> outSSBVar(rMax-rMiddle, 
											 vector<vector<GRBVar>>(4,
											 vector<GRBVar>(16)));

	for(uint r = 0; r < uint(rMax-rMiddle); r++){
		auto & inSSBVar_r = inSSBVar[r];
		auto & outSSBVar_r = outSSBVar[r];
		for(uint i = 0; i < 4; i++){
			auto & inSSBVar_r_i = inSSBVar_r[i];
			auto & outSSBVar_r_i = outSSBVar_r[i];
			for(uint j = 0; j < 16; j++){
				uint bit = 4*(p_in[4*i+j/4]) + j%4;
				uint bit1 = 4*(p_out[4*i+j/4]) + j%4;
				inSSBVar_r_i[j] = m.getVarByName("x_"+to_string(rMiddle+r)+"_"+to_string(bit));
				outSSBVar_r_i[j] = m.getVarByName("v_"+to_string(rMiddle+r)+"_"+to_string(bit1));
			}
		}
	}

	//Time limit
	m.set(GRB_DoubleParam_TimeLimit, 120);


	m.set(GRB_IntParam_PoolSearchMode, 2);
	m.set(GRB_IntParam_PoolSolutions, rMiddle+5);
	m.set(GRB_IntParam_SolutionLimit, rMiddle+5);

	//MIPFocus
	m.set(GRB_IntParam_MIPFocus, DYNAMIC_STEP1_MIPFOCUS);

	callbackDynamic cb(rMiddle, output, valK, XmidVar, KmidVar, inSSBVar, outSSBVar);
	m.setCallback(&cb);	
	//Optimize
	m.optimize();

	cout << m.get(GRB_IntAttr_Status)  << "status code!" << endl;

	if(m.get(GRB_IntAttr_Status) == GRB_INFEASIBLE){
		cout << "Infeasible for rMiddle = " << rMiddle << ", no solution" << endl;
		// if(indexInput != -1) //Infeasible with this input
		// 	return make_tuple(vector<uint8_t>(), vector<uint8_t>(), vector<vector<uint8_t>());
		// else
		// 	break;
		//if(rMiddle == rMax_minus_2)
			//possibleInput.erase(possibleInput.begin()+indexSelectedInput);
		return make_tuple(vector<uint8_t>(), vector<uint8_t>(), vector<vector<uint8_t>>());
	}

	if(m.get(GRB_IntAttr_SolCount) == 0){
		//Didn't found any solution, start over
		cout << "Could not find a solution for rMiddle = " << rMiddle << ", time or solution number or only even trails" << endl;
		return make_tuple(vector<uint8_t>(), vector<uint8_t>(), vector<vector<uint8_t>>());
	}

	uint nbtrails = m.get(GRB_IntAttr_SolCount);
	cout << "the number of solutions is:" << nbtrails << endl;


	stack<tuple<int, vector<uint8_t>, vector<uint8_t>>> tree;

	for(uint i = 0; i < nbtrails; i++)
	{
		vector<uint8_t> x_in(64);
		vector<uint8_t> k_in(64);

		m.set(GRB_IntParam_SolutionNumber, nbtrails-(i+1));
		for(uint j = 0; j < 64; j++)
		{
			x_in[j] = uint8_t(round(XmidVar[j].get(GRB_DoubleAttr_X)));
			k_in[j] = uint8_t(round(KmidVar[j].get(GRB_DoubleAttr_X)));
		}
		tree.push(make_tuple(rMiddle, x_in, k_in));
	}


	while(!tree.empty())
	{

		auto element = tree.top();
		tree.pop();

		for(uint i = 0; i < 64; i++)
		{
			valK[get<0>(element)][i] = get<2>(element)[i];
			valX[get<0>(element)][i] = get<1>(element)[i];
		}

		int hw = 0;

		for(uint i = 0; i < 64; i++)
			hw += int(get<1>(element)[i]);

		cout << hw << "  hw  " << endl;
		if(hw >= 40)  //this argument should be improved
		{

			uint goodrMiddle = get<0>(element);

			cout << "Found a good middle round, starting from rMiddle = " << goodrMiddle << endl;
			cout << "Input : "; printParitySetHex(input); cout << endl;
			cout << "Keys : " << endl;
			for(int r = goodrMiddle; r < rMax_minus_1; r++){
				printParitySetHex(valK[r]); cout << endl;
			}
			cout << endl;
			cout << "Searching for key pattern from input to middle" << endl;

			auto Lsol_allsol = searchkey(goodrMiddle,input,valX[goodrMiddle]);
			auto const & Lsol = Lsol_allsol.first;
			if(Lsol.size() > 0){
				//We found at least one solution for the first half
				for(auto const & firstKeys : Lsol){

					//Complete the current key pattern
					for(uint r = 0; r < goodrMiddle; r++){
						auto & valK_r = valK[r];
						auto const & firstKeys_r = firstKeys[r];
						for(uint i = 0; i < blockSize; i++)
							valK_r[i] = firstKeys_r[i];
					}

					cout << "Input  : "; printParitySetHex(input); cout << endl;
					cout << "Output : "; printParitySetHex(output); cout << endl;
					cout << "Keys   : " << endl;
					for(auto const & keypattern : valK){
						printParitySetHex(keypattern);
						cout << endl;
					}
					cout << "Counting trails..." << endl;
					//Now count the actual number of trails

					auto nbTrail_fullOpti = countTrailsFullCipher(input,output,valK);
					auto nbTrail = nbTrail_fullOpti.first;
					cout << nbTrail << " trails" << endl;

					auto nbTrail_inter = countTrailsFullCipher(input,output_inter,valK);
					auto inter = nbTrail_inter.first;
					cout << inter << " interfering trails" << endl;


					if(nbTrail_fullOpti.second && nbTrail%2 == 1 && nbTrail_inter.second && inter%2 == 0)
						return make_tuple(input,output,valK);
					else if(nbTrail_fullOpti.second && nbTrail%2 == 0)
						cout << "But even..." << endl;
					//else if(nbTrail_inter.second && nbTrail%2 == 1)
					//	cout << "interfering trails is odd" << endl;
					else
						cout << "But was not able to compute everything..." << endl;
					//Else try next key
				}
			}

			continue;
		}

		else
		{
			uint rMiddle = get<0>(element) - 1;
			GRBModel m(gurobiEnv,"./modeltwine_SSB/"+to_string(rounds)+"_twine.mps");

			m.set(GRB_IntParam_LazyConstraints, 1);


			//If needed, enforce the last key to be non zero
			//GRBLinExpr sumLastKey(0);
			//for(uint i = 0; i < blockSize; i++)
			//{
			//	sumLastKey += m.getVarByName("k_"+to_string(rMax_minus_2)+"_"+to_string(i));
			//	sumLastKey += m.getVarByName("k_"+to_string(rMax_minus_2-1)+"_"+to_string(i));
			//}
			//m.addConstr(sumLastKey >= 1);

			//GRBLinExpr sumFirstKey(0);
			//for(uint i = 0; i < blockSize; i++)
			//{
			//	sumFirstKey += m.getVarByName("k_0_"+to_string(i));
			//	sumFirstKey += m.getVarByName("k_1_"+to_string(i));
			//}
			//m.addConstr(sumFirstKey >= 1);


			//Fix the input
			for (uint i = 0; i < 64; i++)
				m.addConstr(m.getVarByName("x_0_"+to_string(i)) == input[i]);

			//Fix the output
			for(uint i = 0; i < blockSize; i++)
				m.addConstr(m.getVarByName("v_"+to_string(rMax_minus_1)+"_"+to_string(i)) == output[i]);


			//Fix the known keys
			for(int r = rMiddle+1; r < rMax_minus_1; r++){
				auto const & valK_r = valK[r];
				for(uint i = 0; i < blockSize; i++)
					m.addConstr(m.getVarByName("k_"+to_string(r)+"_"+to_string(i)) == valK_r[i]);
			}

			//Fix the known round input
			for(int r = rMiddle+1; r < rMax_minus_1; r++){
				auto const & valX_r = valX[r];
				for(uint i = 0; i < blockSize; i++)
					m.addConstr(m.getVarByName("x_"+to_string(r)+"_"+to_string(i)) == valX_r[i]);
			}

			//Objective variables
			vector<GRBVar> XmidVar(blockSize);
			vector<GRBVar> KmidVar(blockSize);
			for(uint i = 0; i < blockSize; i++){
				XmidVar[i] = m.getVarByName("x_"+to_string(rMiddle)+"_"+to_string(i));
				KmidVar[i] = m.getVarByName("k_"+to_string(rMiddle)+"_"+to_string(i));
			}

			//Set the objective
			GRBLinExpr objExpr(0);
			for(uint i = 0; i < blockSize; i++){
				objExpr	+= XmidVar[i];
				objExpr -= m.getVarByName("v_"+to_string(rMiddle+1)+"_"+to_string(i));
				objExpr += 0.01*KmidVar[i];
			}
			m.setObjective(objExpr, GRB_MAXIMIZE);

			//Seed
			auto seedGurobi = randintSeedGurobi(prng);
			m.set(GRB_IntParam_Seed, seedGurobi);

			//Prepare the callback variables
			//SSB variables for the second half
			vector<vector<vector<GRBVar>>> inSSBVar(rMax-rMiddle, 
													vector<vector<GRBVar>>(4,
													vector<GRBVar>(16)));
			vector<vector<vector<GRBVar>>> outSSBVar(rMax-rMiddle, 
													 vector<vector<GRBVar>>(4,
													 vector<GRBVar>(16)));

			for(uint r = 0; r < rMax-rMiddle; r++)
			{
				auto & inSSBVar_r = inSSBVar[r];
				auto & outSSBVar_r = outSSBVar[r];
				for(uint i = 0; i < 4; i++)
				{
					auto & inSSBVar_r_i = inSSBVar_r[i];
					auto & outSSBVar_r_i = outSSBVar_r[i];
					for(uint j = 0; j < 16; j++)
					{
						uint bit = 4*(p_in[4*i+(j/4)]) + (j%4);
						uint bit1 = 4*(p_out[4*i+(j/4)]) + (j%4);

						inSSBVar_r_i[j] = m.getVarByName("x_"+to_string(rMiddle+r)+"_"+to_string(bit));
						outSSBVar_r_i[j] = m.getVarByName("v_"+to_string(rMiddle+r)+"_"+to_string(bit1));
					}
				}
			}
					
				
			


				//Time limit
			m.set(GRB_DoubleParam_TimeLimit, 180);


			//int NB_SOLUTION;
			//if (rMiddle <= 6)
			//	NB_SOLUTION = 1;
			//else
			if(rMiddle >= 6)
			{
				m.set(GRB_IntParam_PoolSearchMode, 2);
				m.set(GRB_IntParam_PoolSolutions, rMiddle);
				m.set(GRB_IntParam_SolutionLimit, 180);
			}

			//MIPFocus
			m.set(GRB_IntParam_MIPFocus, DYNAMIC_STEP1_MIPFOCUS);

			callbackDynamic cb(rMiddle, output, valK, XmidVar, KmidVar, inSSBVar, outSSBVar);
			m.setCallback(&cb);	
			//Optimize
			m.optimize();

			if(m.get(GRB_IntAttr_Status) == GRB_INFEASIBLE){
				cout << "Infeasible for rMiddle = " << rMiddle << ", backtracking" << endl;
				// if(indexInput != -1) //Infeasible with this input
				// 	return make_tuple(vector<uint8_t>(), vector<uint8_t>(), vector<vector<uint8_t>());
				// else
				// 	break;
				//if(rMiddle == rMax_minus_2)
					//possibleInput.erase(possibleInput.begin()+indexSelectedInput);
				continue;
			}

			if(m.get(GRB_IntAttr_SolCount) == 0){
				//Didn't found any solution, start over
				cout << "Could not find a solution for rMiddle = " << rMiddle << ", backtracking" << endl;
				continue;
			}

			uint nbtrails = m.get(GRB_IntAttr_SolCount);
			cout << "the number of trails is:" << nbtrails << endl;



			for(uint i = 0; i < nbtrails; i++)
			{
				vector<uint8_t> x_in(64);
				vector<uint8_t> k_in(64);
				m.set(GRB_IntParam_SolutionNumber, nbtrails- (i+1));
				for(uint j = 0; j < 64; j++)
				{
					x_in[j] = uint8_t(round(XmidVar[j].get(GRB_DoubleAttr_X)));
					k_in[j] = uint8_t(round(KmidVar[j].get(GRB_DoubleAttr_X)));
				}
				tree.push(make_tuple(rMiddle, x_in, k_in));
			}


		}

	}
	//If we get here, then no input can lead to at least one trail

	cout << " all element have been used, need increase time" << endl;
	return make_tuple(vector<uint8_t>(), vector<uint8_t>(), vector<vector<uint8_t>>());
}



pair<uint64_t,int> counttrail_middle(uint rounds,
		    vector<uint8_t> & u,
                    vector<uint8_t> & v,
                    vector<vector<uint8_t>> & k)
{
    GRBEnv gurobiEnv;
    gurobiEnv.set(GRB_IntParam_OutputFlag, 0);
    GRBModel m(gurobiEnv, to_string(rounds)+"_twine.mps");
    m.set(GRB_IntParam_Threads, 8);

    vector<GRBVar> invar(64);
    vector<GRBVar> outvar(64);
    vector<vector<GRBVar>> keyvar(rounds-1, vector<GRBVar>(64));
    for (uint i = 0; i < 64; i++)
    {
        invar[i] = m.getVarByName("x_0_"+to_string(i));
        outvar[i] = m.getVarByName("v_"+to_string(rounds-1)+"_" + to_string(i));

        m.addConstr(invar[i] == u[i]);
        m.addConstr(outvar[i] == v[i]);

    }
    
    vector<GRBVar> in_2(64);
    vector<GRBVar> out_2(64);
    for (uint i = 0; i < 64; i++)
    {
	    in_2[i] = m.getVarByName("x_1_" + to_string(i));
	    out_2[i] = m.getVarByName("v_8_" + to_string(i));
    }


    for(uint r = 0; r < rounds - 1;r++)
    { 
	    for (uint i = 0; i < 64; i++)
	    {
		keyvar[r][i] = m.getVarByName("k_"+to_string(r)+"_"+to_string(i));

		m.addConstr(keyvar[r][i] == k[r][i]);
	    }
    }

    m.set(GRB_IntParam_PoolSearchMode, 2);
    m.set(GRB_IntParam_PoolSolutions , 2000000000);
    m.set(GRB_IntParam_SolutionLimit, 10000000);
    m.set(GRB_DoubleParam_TimeLimit, 1800);
    m.optimize();
    int optimized;

    uint const nbsol = m.get(GRB_IntAttr_SolCount);

    if(nbsol > 0)
    {
       
        vector<vector<uint8_t>> in(nbsol, vector<uint8_t>(64));
        vector<vector<uint8_t>> out(nbsol, vector<uint8_t>(64));

        for(uint indexSol = 0; indexSol < nbsol; indexSol++){
			m.set(GRB_IntParam_SolutionNumber, indexSol);
				for(uint i = 0; i < 64; i++)
				{
					in[indexSol][i] = uint8_t(round(in_2[i].get(GRB_DoubleAttr_Xn)));
					out[indexSol][i] = uint8_t(round(out_2[i].get(GRB_DoubleAttr_Xn)));
				}
		}

	for(uint i = 0; i < in.size(); i++)
	{
		for(uint j = 0; j < 64; j++)
		{
			//cout << int(in[i][j]);
		}
		//cout << endl;
	}


	cout << "************" << endl;
	for(uint i = 0; i < in.size(); i++)
	{
		for(uint j = 0; j < 64; j++)
		{
			//cout << int(out[i][j]);
		}
		//cout << endl;
	}
    }



    uint64_t nbtrail = 0;
    if(m.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        nbtrail = m.get(GRB_IntAttr_SolCount);
	optimized = int(m.get(GRB_IntAttr_Status));
    }
    //else if(m.get(GRB_IntAttr_Status) == GRB_SOLUTION_LIMIT)
    //    optimized = false;
    else
    {
	    cout << "__reason:" << m.get(GRB_IntAttr_Status) << endl;
	    optimized = int(m.get(GRB_IntAttr_Status));
    }
    return make_pair(nbtrail, optimized);
}






pair<uint64_t,bool> counttrail(vector<uint8_t> & u,
                    vector<uint8_t> & v,
                    vector<uint8_t> & k)
{
    GRBEnv gurobiEnv;
    gurobiEnv.set(GRB_IntParam_OutputFlag, 0);
    GRBModel m(gurobiEnv, "first.mps");
    m.set(GRB_IntParam_Threads, 8);

    vector<GRBVar> invar(64);
    vector<GRBVar> outvar(64);
    vector<GRBVar> keyvar(64);
    for (uint i = 0; i < 64; i++)
    {
        invar[i] = m.getVarByName("x_0_"+to_string(i));
        outvar[i] = m.getVarByName("x_1_" + to_string(i));

        m.addConstr(invar[i] == u[i]);
        m.addConstr(outvar[i] == v[i]);

    }
    for (uint i = 0; i < 64; i++)
    {
        keyvar[i] = m.getVarByName("k_0_"+to_string(i));

        m.addConstr(keyvar[i] == k[i]);
    }

    m.set(GRB_IntParam_PoolSearchMode, 2);
    m.set(GRB_IntParam_PoolSolutions , 2000000000);
    m.set(GRB_IntParam_SolutionLimit, 100000);
    m.set(GRB_DoubleParam_TimeLimit, 60);
    m.optimize();
    bool optimized = true;
    uint64_t nbtrail = 0;
    if(m.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
        nbtrail = m.get(GRB_IntAttr_SolCount);
    else if(m.get(GRB_IntAttr_Status) == GRB_SOLUTION_LIMIT)
        optimized = false;
    return make_pair(nbtrail, optimized);
}


pair<uint64_t,bool> counttrail_last(vector<uint8_t> & u,
                    vector<uint8_t> & v,
                    vector<uint8_t> & k)
{
    GRBEnv gurobiEnv;
    gurobiEnv.set(GRB_IntParam_OutputFlag, 0);
    GRBModel m(gurobiEnv, "last.mps");
    m.set(GRB_IntParam_Threads, 8);

    vector<GRBVar> invar(64);
    vector<GRBVar> outvar(64);
    vector<GRBVar> keyvar(64);
    for (uint i = 0; i < 64; i++)
    {
        invar[i] = m.getVarByName("x_0_"+to_string(i));
        outvar[i] = m.getVarByName("x_1_" + to_string(i));

        m.addConstr(invar[i] == u[i]);
        m.addConstr(outvar[i] == v[i]);

    }
    for (uint i = 0; i < 64; i++)
    {
        keyvar[i] = m.getVarByName("k_0_"+to_string(i));

        m.addConstr(keyvar[i] == k[i]);
    }

    m.set(GRB_IntParam_PoolSearchMode, 2);
    m.set(GRB_IntParam_PoolSolutions , 2000000000);
    m.set(GRB_IntParam_SolutionLimit, 100000);
    m.set(GRB_DoubleParam_TimeLimit, 60);
    m.optimize();
    bool optimized = true;
    uint64_t nbtrail = 0;
    if(m.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
        nbtrail = m.get(GRB_IntAttr_SolCount);
    else if(m.get(GRB_IntAttr_Status) == GRB_SOLUTION_LIMIT)
        optimized = false;
    return make_pair(nbtrail, optimized);
}





pair<vector<vector<uint8_t>>, vector<vector<uint8_t>>> countfirst(vector<uint8_t> &u)
{
	GRBEnv gurobiEnv;
	gurobiEnv.set(GRB_IntParam_OutputFlag, 0);
	GRBModel m(gurobiEnv, "first.mps");
	m.set(GRB_IntParam_Threads, 8);

        m.set(GRB_IntParam_LazyConstraints, 1);
	m.set(GRB_IntParam_MIPFocus, 1);

	vector<GRBVar> invar(64);
	vector<GRBVar> outvar(64);
	vector<GRBVar> keyvar(64);

	for(uint i = 0; i < 64; i++)
	{
		invar[i] = m.getVarByName("x_0_"+to_string(i));
		keyvar[i] = m.getVarByName("k_0_"+to_string(i));
		outvar[i] = m.getVarByName("x_1_"+to_string(i));
	}
	for(uint i = 0; i < 64; i++)
	{
		m.addConstr(invar[i] == u[i]);

	}
	GRBLinExpr sumlastkey(0);
	for(uint i = 0; i < 64; i++)
		sumlastkey += m.getVarByName("k_0_" + to_string(i));
	m.addConstr(sumlastkey >= 1);

	GRBLinExpr sumoutput(0);
	for(uint i = 0; i < 64; i++)
		sumoutput += m.getVarByName("x_1_" + to_string(i));
	m.addConstr(sumoutput <= 63);


        GRBLinExpr objExpr(0);
	for(uint i = 0; i < 64; i++)
	{
	    objExpr += invar[i];
	}
        m.setObjective(objExpr, GRB_MAXIMIZE);
	
   	m.set(GRB_IntParam_PoolSearchMode, 2);
       	m.set(GRB_IntParam_PoolSolutions, MAX_NUMBER_SEARCHKEY);
	m.set(GRB_DoubleParam_TimeLimit, TIME_LIMIT_SEARCHKEY);
       	//m.set(GRB_DoubleParam_MIPGapAbs,MAX_GRBMIPGapAbs);

	callbackfirst cb (invar, outvar, keyvar);
	m.setCallback(&cb);
	m.update();
	m.optimize();

	uint const nbso = m.get(GRB_IntAttr_SolCount);
	if(nbso > 0)
		{

			vector<vector<uint8_t>> allso(nbso, vector<uint8_t>(64));
			vector<vector<uint8_t>> allkey(nbso, vector<uint8_t>(64));
			for(uint i = 0; i < nbso; i++)
			{
				m.set(GRB_IntParam_SolutionNumber, i);
				for(uint j = 0; j < 64; j++)
				{
					allso[i][j] = uint8_t(round(outvar[j].get(GRB_DoubleAttr_Xn)));
					allkey[i][j] = uint8_t(round(keyvar[j].get(GRB_DoubleAttr_Xn)));
				}
			}
			return make_pair(allso, allkey);
		}
		else
			return make_pair(vector<vector<uint8_t>>(), vector<vector<uint8_t>>());
}

pair<vector<vector<uint8_t>>, vector<vector<uint8_t>>> countlast(vector<uint8_t> &u)
{
	GRBEnv gurobiEnv;
	gurobiEnv.set(GRB_IntParam_OutputFlag, 0);
	GRBModel m(gurobiEnv, "last.mps");
	m.set(GRB_IntParam_Threads, 8);
        m.set(GRB_IntParam_LazyConstraints, 1);

	vector<GRBVar> invar(64);
	vector<GRBVar> outvar(64);
	vector<GRBVar> keyvar(64);

	for(uint i = 0; i < 64; i++)
	{
		invar[i] = m.getVarByName("x_0_"+to_string(i));
		keyvar[i] = m.getVarByName("k_0_"+to_string(i));
		outvar[i] = m.getVarByName("x_1_"+to_string(i));
	}
	for(uint i = 0; i < 64; i++)
	{
		m.addConstr(outvar[i] == u[i]);

	}
	GRBLinExpr sumlastkey(0);
	for(uint i = 0; i < 64; i++)
		sumlastkey += m.getVarByName("k_0_" + to_string(i));
	m.addConstr(sumlastkey >= 1);

	GRBLinExpr suminput(0);
	for(uint i = 0; i < 64; i++)
		suminput += m.getVarByName("x_0_" + to_string(i));
	m.addConstr(suminput >= 1);


        GRBLinExpr objExpr(0);
	for(uint i = 0; i < 64; i++)
	{
	    objExpr += outvar[i];
	}
        m.setObjective(objExpr, GRB_MAXIMIZE);

   	m.set(GRB_IntParam_PoolSearchMode, 2);
       	m.set(GRB_IntParam_PoolSolutions, MAX_NUMBER_SEARCHKEY);
	m.set(GRB_DoubleParam_TimeLimit, TIME_LIMIT_SEARCHKEY);
       	m.set(GRB_DoubleParam_MIPGapAbs,MAX_GRBMIPGapAbs);

	callbacklast cb (invar, outvar, keyvar);
	m.setCallback(&cb);
	m.update();
	m.optimize();

	uint const nbso = m.get(GRB_IntAttr_SolCount);
	if(nbso > 0)
		{

			vector<vector<uint8_t>> allso(nbso, vector<uint8_t>(64));
			vector<vector<uint8_t>> allkey(nbso, vector<uint8_t>(64));
			for(uint i = 0; i < nbso; i++)
			{
				m.set(GRB_IntParam_SolutionNumber, i);
				for(uint j = 0; j < 64; j++)
				{
					allso[i][j] = uint8_t(round(invar[j].get(GRB_DoubleAttr_Xn)));
					allkey[i][j] = uint8_t(round(keyvar[j].get(GRB_DoubleAttr_Xn)));
				}
			}
			return make_pair(allso, allkey);
		}
		else
			return make_pair(vector<vector<uint8_t>>(), vector<vector<uint8_t>>());
}



void printParitySetHex(vector<uint8_t> const & x) 
{
	/*
		Hex printer for more compact prints
		LSB of x is x[0]
		An hex character thus represent (x[i],...,x[i+3]) with x[i] LSB
		e.g.
		x = {1,0,0,0} is printed as 1
		x = {0,1,0,0,1,1,0,1} is printed as 2B
	*/
	for(uint s = 0; s < 16; s++){
		uint val = 0;
		for(uint i = 0; i < 4; i++)
			val |= (x[4*s + i] << i);
		cout << hex << uppercase << val << dec;
	}
}


void printBits(uint64_t const x, uint nbBits){
	//print the first #nbBits of x in cout
	uint64_t y = x;
	for(uint i = 0; i < nbBits; i++){
		cout << (y & 1);
		y >>= 1;
	}
}

void printVec(std::vector<uint8_t> const & x){
	//print the elements of x, wihtout separators
	for(auto const & xx : x) cout << int(xx);
}

string hex_char_to_bin(char c)
{
	//hex to bin conversion, LSB on the left in the binary representation
    switch(toupper(c))
    {
        case '0': return "0000";
        case '1': return "1000";
        case '2': return "0100";
        case '3': return "1100";
        case '4': return "0010";
        case '5': return "1010";
        case '6': return "0110";
        case '7': return "1110";
        case '8': return "0001";
        case '9': return "1001";
        case 'A': return "0101";
        case 'B': return "1101";
        case 'C': return "0011";
        case 'D': return "1011";
        case 'E': return "0111";
        case 'F': return "1111";
        default : return ""+c;
    }
}

vector<uint8_t> hexToVec(string const & s){
	/*
		From an hex string get the corresponding vector of bits
		LSB of each hex char is put on the left
		e.g. 
		s = "1"  returns {1,0,0,0}
		s = "2B" returns {0,1,0,0,1,1,0,1}
	*/

    string bin;
    for(unsigned i = 0; i != s.size(); ++i)
       bin += hex_char_to_bin(s[i]);
    
    vector<uint8_t> v(bin.size());
    for(uint i = 0; i < v.size(); i++){
    	if(bin[i] == '0') v[i] = 0;
    	else v[i] = 1;
    }

    return v;
}


