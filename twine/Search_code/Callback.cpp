#include "Callback.hpp"

#define MAX_GRBMIPGapAbs 8						/*Accept non-optimal solutions in searchKeyPattern if the gap between the objective and the best bound on the objective is less than MAX_GRBMIPGapAbs*/
#define MAX_NUMBER_SEARCHKEY 10					/*Maximum number of key patterns to search in searchKeyPattern*/
#define TIME_LIMIT_SEARCHKEY 300				/*Time limit (in seconds) in searchKeyPattern*/
#define MAX_GRB_THREAD 8						/*Max number of thread to use by Gurobi. Setting it to 0 let Gurobi decide it itself*/
#define ENFORCE_NONZERO_LASTKEY true			/*If true, enforce the last key to be non-zero. NECESSARY for the min-degree search, sometimes improve the speed for the degree search. Most results on the algebraic degree were obtained with this parameter set to false*/

#define WEIGHT_THRESHOLD_STEP1 40				/*Threshold on the wheight of the input of a round to go to the next step in the dynamic search algorithm*/
#define TIME_LIMIT_ONESTEPDYNAMIC 30			/*Time limit (in seconds) in improvedDynamicSearch to find one iteration in the first step*/
#define DYNAMIC_STEP1_MIPFOCUS 1				/*Hint for Gurobi to focus more on finding solution than proving optimality. Seems to improve the speed, see the Gurobi doc for more details*/

#define AGRESSIVE_STRATEGY_MINDEGREE false		/*Define how a number of trails > 0 should be treated in the min degree search for non-diagonal coefficient in the matrix. false for GIFT was enough, but for Present, setting this to true + using the dedicated last layer is better*/

#define MAX_NUMBER_SOLCOUNT_FULLCIPHER 100000	/*Solution limit when counting the total number of trails over the whole cipher in countTrailsFullCipher*/
#define TIME_LIMIT_COUNTTRAILS_FULLCIPER 60	/*Time limit (in seconds) when counting the total number of trails over the whole cipher in countTrailsFullCipher. Set to GRB_INFINITY for unlimited time*/

/*The values of the last two limits are the one used to prove the min-degree.
For the algebraic degree for PRESENT, higher limits were used but were not saved.
It should be something like 1 000 000 max solutions and unlimited time
*/

/*
	Macros are used in the code to define different limits on the search (e.g. time limits)
	It's a bit dirty but easier to change and than handling a ton of parameters everywhere
	Our results were obtained with the current values
*/


#define TIME_LIMIT_CBSEARCHKEY 120		/*Time limit (in seconds) when counting the actual number of solutions in callbackSearchKeyPattern*/
#define MAX_NUMBER_CBSEARCHKEY 30000	/*Max number of solution allowed when counting the actual number of solutions in callbackSearchKeyPattern*/
#define MAX_GRB_THREAD 8 				/*Max number of thread to use by Gurobi. Setting it to 0 let Gurobi decide it itself*/

#define TIME_LIMIT_CALLBACKDYNAMIC 60	/*Time limit (in seconds) when counting the actual number of solutions in callbackDynamic*/
#define THRESHOLD_NBSOL_DYNAMIC 21		/*Max number of solution allowed when counting the actual number of solutions in callbackDynamic*/
#define THRESHOLD_NBSOL_SSB_DYNAMIC 1001	/*Max number of solution allowed in one SSB when counting the actual number of solutions in callbackDynamic*/




using namespace std;
using namespace boost;

typedef unsigned int uint;
typedef pair<uint32_t,uint32_t> pairuint32;
typedef tuple<uint32_t,uint32_t,uint64_t> tup3uint;


callbackDynamic::callbackDynamic(uint const xrMiddle,
                    			 vector<uint8_t> const & xoutput,
                    			 vector<vector<uint8_t>> const & xkeyVal,
                    			 vector<GRBVar> const & xinputMiddleVar,
                    			 vector<GRBVar> const & xkeyMiddleVar,
                    			 vector<vector<vector<GRBVar>>> const & xinSSBVar,
                    			 vector<vector<vector<GRBVar>>> const & xoutSSBVar) :
	rMiddle(xrMiddle),
	output(xoutput),
	keyVal(xkeyVal),
	inputMiddleVar(xinputMiddleVar),
	keyMiddleVar(xkeyMiddleVar),
	inSSBVar(xinSSBVar),
	outSSBVar(xoutSSBVar)
{}

void callbackDynamic::callback(){
	try {
        if (where == GRB_CB_MIPSOL) { //If the solver found a solution
		    
        	//Get the solution for each SSB
		cout << "entering the callback" << endl;

        	vector<vector<vector<uint8_t>>> valInSSB(inSSBVar.size(),
        											vector<vector<uint8_t>>(4,
        											vector<uint8_t>(16)));
        	vector<vector<vector<uint8_t>>> valOutSSB(outSSBVar.size(),
        											vector<vector<uint8_t>>(4,
        											vector<uint8_t>(16)));
        	for(uint r = 0; r < inSSBVar.size(); r++){
        		auto & valInSSB_r = valInSSB[r];
        		auto & valOutSSB_r = valOutSSB[r];

        		auto & inSSBVar_r = inSSBVar[r];
        		auto & outSSBVar_r = outSSBVar[r];
        		for(uint i = 0; i < 4; i++){
        			auto & valInSSB_r_i = valInSSB_r[i];
	        		auto & valOutSSB_r_i = valOutSSB_r[i];

	        		auto & inSSBVar_r_i = inSSBVar_r[i];
	        		auto & outSSBVar_r_i = outSSBVar_r[i];
        			for(uint j = 0; j < 16; j++){
        				valInSSB_r_i[j] = int(round(getSolution(inSSBVar_r_i[j])));
        				valOutSSB_r_i[j] = int(round(getSolution(outSSBVar_r_i[j])));
        			}
        		}
        	}

        	//Check if this key pattern can lead to an odd number of trails by checking all SSB
        	bool oddSSB = true;
        	cout << "Number of trails in the SSBs : ";
        	for(uint r = 0; r < inSSBVar.size(); r++){
        		auto & valInSSB_r = valInSSB[r];
        		auto & valOutSSB_r = valOutSSB[r];
        		for(uint i = 0; i < 4; i++){
        			auto & valInSSB_r_i = valInSSB_r[i];
        			auto & valOutSSB_r_i = valOutSSB_r[i];
        			uint64_t nbTrail;
        				
				nbTrail = countSSBTrails(valInSSB_r_i, valOutSSB_r_i);
        			cout << nbTrail << " ";
        			//Check if it exceeds the threshold or even
        			if(nbTrail > THRESHOLD_NBSOL_SSB_DYNAMIC || nbTrail%2 == 0){
        				oddSSB = false;
        				auto & inSSBVar_r_i = inSSBVar[r][i];
		        		auto & outSSBVar_r_i = outSSBVar[r][i];
        				GRBLinExpr cutExpr(0);
        				for(uint j = 0; j < 16; j++){
        					if(valInSSB_r_i[j] == 0) cutExpr += inSSBVar_r_i[j];
        					else cutExpr += (1 - inSSBVar_r_i[j]);
        					if(valOutSSB_r_i[j] == 0) cutExpr += outSSBVar_r_i[j];
        					else cutExpr += (1 - outSSBVar_r_i[j]);
        				}
        				addLazy(cutExpr >= 1);
        			}
	        	}
	        }
	        // cout << endl;

	        if(oddSSB){
	        	//All SSB in the second half are odd, count the actual number of trail in the second half
	        	//Generate a model for the second half
			GRBEnv gurobiEnv;
			gurobiEnv.set(GRB_IntParam_OutputFlag, 0);

	        	//Read the model
	        	GRBModel m(gurobiEnv, "./modeltwine_SSB/"+to_string(inSSBVar.size())+"_twine.mps");

	        	//Add the input/output constraints
	        	for(uint i = 0; i < 64; i++){
	        		m.addConstr(m.getVarByName("x_0_"+to_string(i)) == int(round(getSolution(inputMiddleVar[i]))));
	        		m.addConstr(m.getVarByName("v_"+to_string(inSSBVar.size() - 1)+"_"+to_string(i)) == output[i]);
	        	}

	        	//Add the key constraints
	        	//First key is extracted from current solution
	        	vector<uint8_t> keyMiddleVal(64);
	        	for(uint i = 0; i < 64; i++)
	        		keyMiddleVal[i] = int(round(getSolution(keyMiddleVar[i])));
	        	for(uint i = 0; i < 64; i++)
	        		m.addConstr(m.getVarByName("k_0_"+to_string(i)) == keyMiddleVal[i]);
	        	//The remaining keys are already known
	        	for(uint r = rMiddle+1; r < keyVal.size(); r++){
	        		uint const offset = r - rMiddle;
	        		auto const & keyVal_r = keyVal[r];
	        		for(uint i = 0; i < 64; i++)
	        			m.addConstr(m.getVarByName("k_"+to_string(offset)+"_"+to_string(i)) == keyVal_r[i]);
	        	}

	        	//Dummy objective
				GRBLinExpr objExpr(0);
				for(uint i = 0; i < 64; i++){
					objExpr += m.getVarByName("x_0_"+to_string(i));
				}
				m.setObjective(objExpr);

				//setup params for enumerating solutions
	        	cout << "Counting the actual number of trails... in Callback" << endl;
	        	m.set(GRB_IntParam_PoolSearchMode, 2);
	        	m.set(GRB_IntParam_PoolSolutions, 2000000000); //maximum number of solution allowed by Gurobi
	        	uint64_t threshold = 1;
	        	if(THRESHOLD_NBSOL_DYNAMIC > 1)
	        		threshold = (keyVal.size() + 7 - rMiddle)*THRESHOLD_NBSOL_DYNAMIC;
	        	m.set(GRB_IntParam_SolutionLimit, threshold+1); //We stop if we find more solutions than the threshold
	        	m.set(GRB_IntParam_Threads,MAX_GRB_THREAD);
			m.set(GRB_DoubleParam_TimeLimit, TIME_LIMIT_CALLBACKDYNAMIC);

				m.optimize();

				if(m.get(GRB_IntAttr_Status) == GRB_OPTIMAL){
					uint nbTrail = m.get(GRB_IntAttr_SolCount);
					cout << nbTrail << " different trails" << endl;

					//If the actual number of trails is even or larger than the threshold, we need another key
					if(nbTrail > threshold || nbTrail%2 == 0) oddSSB = false;
				}
				else{
					cout << "Too many or no solutions (" << m.get(GRB_IntAttr_SolCount) << ")" << endl;
					oddSSB = false;
				}

				if(!oddSSB){
			        //If we get here, the key pattern is not good, search for another one
			        //Generate a linear expression to remove he current key pattern from the solution pool
			        GRBLinExpr cutExpr(0);
			        for(uint i = 0; i < 64; i++){
			        	if(keyMiddleVal[i] == 0) cutExpr += keyMiddleVar[i];
			        	else cutExpr += (1 - keyMiddleVar[i]);
			        }
		        	addLazy(cutExpr >= 1);
		        }
	        }

		}
	} catch (GRBException e) {
    cout << "Error number: " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
	} catch (...) {
	cout << "Error during callback" << endl;
	}
}



callbacksearchkey::callbacksearchkey(uint const xnbRounds,
						vector<GRBVar> const & xinputVar,
						vector<GRBVar> const & xoutputVar,
						vector<vector<GRBVar>> const & xallKeyVar,
                        vector<vector<vector<GRBVar>>> const & xinSSBVar,
                        vector<vector<vector<GRBVar>>> const & xoutSSBVar) :
	nbRounds(xnbRounds),
	inputVar(xinputVar),
	outputVar(xoutputVar),
	allKeyVar(xallKeyVar),
	inSSBVar(xinSSBVar),
	outSSBVar(xoutSSBVar),
	ctrsolutions(0),
	allSolutions()
{}

void callbacksearchkey::callback(){
	try {
        if (where == GRB_CB_MIPSOL) { //If the solver found a solution

        	//Increment the solution counter
        	ctrsolutions++;
        	//cout << ctrsolutions << " solutions up to now\r" << flush;

        	//Get the solution for each SSB

		uint number = inSSBVar.size();

        	vector<vector<vector<uint8_t>>> valInSSB(number,
        											vector<vector<uint8_t>>(4,
        											vector<uint8_t>(16)));
        	vector<vector<vector<uint8_t>>> valOutSSB(number,
        											vector<vector<uint8_t>>(4,
        											vector<uint8_t>(16)));
        	for(uint r = 0; r < number; r++){
        		auto & valInSSB_r = valInSSB[r];
        		auto & valOutSSB_r = valOutSSB[r];

        		auto & inSSBVar_r = inSSBVar[r];
        		auto & outSSBVar_r = outSSBVar[r];
        		for(uint i = 0; i < 4; i++){
        			auto & valInSSB_r_i = valInSSB_r[i];
	        		auto & valOutSSB_r_i = valOutSSB_r[i];

	        		auto & inSSBVar_r_i = inSSBVar_r[i];
	        		auto & outSSBVar_r_i = outSSBVar_r[i];
        			for(uint j = 0; j < 16; j++){
        				valInSSB_r_i[j] = int(round(getSolution(inSSBVar_r_i[j])));
        				valOutSSB_r_i[j] = int(round(getSolution(outSSBVar_r_i[j])));
        			}
        		}
        	}

        	//Get the full key solution
        	uint const allKeyVar_size = allKeyVar.size();
        	vector<dynamic_bitset<>> keySolution(allKeyVar_size, 
        										 dynamic_bitset<>(64,0));
			for(uint r = 0; r < allKeyVar_size; r++){
        		auto & allKeyVar_r = allKeyVar[r];
        		auto & keySolution_r = keySolution[r];
        		for(uint i = 0; i < 64; i++){
        			keySolution_r[i] = uint(round(getSolution(allKeyVar_r[i])));
        		}
        	}
        	//save the solution
	        allSolutions.emplace_back(keySolution);

        	//Check if this key pattern can lead to an odd number of trails by checking all SSB
        	bool oddSSB = true;
        	for(uint r = 0; r < number; r++){
        		auto & valInSSB_r = valInSSB[r];
        		auto & valOutSSB_r = valOutSSB[r];
        		for(uint i = 0; i < 4; i++){
        			auto & valInSSB_r_i = valInSSB_r[i];
        			auto & valOutSSB_r_i = valOutSSB_r[i];
        			uint64_t nbTrail;

        				nbTrail = countSSBTrails(valInSSB_r_i, valOutSSB_r_i);

        			if(nbTrail%2 == 0){ //Even number of trails in this SSB, remove it
					cout << "remove invalid SSB value" << endl;
        				oddSSB = false;
        				auto & inSSBVar_r_i = inSSBVar[r][i];
		        		auto & outSSBVar_r_i = outSSBVar[r][i];
        				GRBLinExpr cutExpr(0);
        				for(uint j = 0; j < 16; j++){
        					if(valInSSB_r_i[j] == 0) cutExpr += inSSBVar_r_i[j];
        					else cutExpr += (1 - inSSBVar_r_i[j]);
        					if(valOutSSB_r_i[j] == 0) cutExpr += outSSBVar_r_i[j];
        					else cutExpr += (1 - outSSBVar_r_i[j]);
        				}
        				addLazy(cutExpr >= 1);
        			}
	        	}
	        }


	        if(oddSSB){
	        	//All SSB trails are odd, count the actual number of trails
	        	//Read the base model
				GRBEnv gurobiEnv;
				gurobiEnv.set(GRB_IntParam_OutputFlag, 0);
				GRBModel m(gurobiEnv, "./modeltwine_SSB/"+ to_string(nbRounds) +"_tw.mps");
				//Add input/output constraints
				for(uint i = 0; i < 64; i++){
					m.addConstr(m.getVarByName(inputVar[i].get(GRB_StringAttr_VarName)) == int(round(getSolution(inputVar[i]))));
					m.addConstr(m.getVarByName(outputVar[i].get(GRB_StringAttr_VarName)) == int(round(getSolution(outputVar[i]))));

				}
				//Add key constraints
				for(uint r = 0; r < allKeyVar_size; r++){
	        		auto & allKeyVar_r = allKeyVar[r];
	        		auto & keySolution_r = keySolution[r];
	        		for(uint i = 0; i < 64; i++){
	        			m.addConstr(m.getVarByName(allKeyVar_r[i].get(GRB_StringAttr_VarName)) == uint(keySolution_r[i]));
	        		}
	        	}

	        	//setup params for enumerating solutions
	        	cout << endl << "counting the actual number of trails..." << endl;
	        	m.set(GRB_IntParam_PoolSearchMode, 2);
				m.set(GRB_IntParam_PoolSolutions, 2000000000); //maximum number of solution allowed by Gurobi

				//Dummy objective
				GRBLinExpr objExpr(0);
				for(uint i = 0; i < 64; i++){
					objExpr += m.getVarByName("x_0_"+to_string(i));
				}
				m.setObjective(objExpr);
				
				//Set the counting callback and optimize
				// callbackCount cb = callbackCount(MAX_NUMBER_CBSEARCHKEY);
				// m.setCallback(&cb);
				m.set(GRB_DoubleParam_TimeLimit, TIME_LIMIT_CBSEARCHKEY);
				m.set(GRB_IntParam_SolutionLimit, MAX_NUMBER_CBSEARCHKEY);
				m.set(GRB_IntParam_Threads,MAX_GRB_THREAD);
				
				m.optimize();

				if(m.get(GRB_IntAttr_Status) == GRB_OPTIMAL){
					uint nbTrail = m.get(GRB_IntAttr_SolCount);
					cout << endl << nbTrail << " different trails" << endl;

					//If the actual number of trails is even, we need another key
					if(nbTrail%2 == 0) oddSSB = false;
					else{
						uint weightSol = 0;
						cout << endl << "odd solution " << oddSSB << endl;
						for(auto const & ksol : keySolution){
							for(uint i = 0; i < 64; i++){
								cout << ksol[i];
								if(ksol[i] == 1) weightSol++;
							}
							cout << endl;
						}
						cout << "weight " << weightSol << endl;
					}
				}
				else{
					cout << endl << "Too many solutions, aborted" << endl;
					oddSSB = false;
				}

				if(!oddSSB){
			        //If we get here, the key pattern is not good, search for another one
			        //Generate a linear expression to remove he current key pattern from the solution pool
				cout << "delete invalid key value" << endl;
			        GRBLinExpr cutExpr(0);
			        for(uint r = 0; r < allKeyVar_size; r++){
		        		auto & allKeyVar_r = allKeyVar[r];
		        		auto & keySolution_r = keySolution[r];
		        		for(uint i = 0; i < 64; i++){
		        			if(keySolution_r[i] == 0) cutExpr += allKeyVar_r[i];
		        			else cutExpr += (1 - allKeyVar_r[i]);
		        		}
		        	}
		        	addLazy(cutExpr >= 1);
		        }
	        }
	        //If number of trails is odd, the solution is not removed and thus the solver finishes
        }
      } catch (GRBException e) {
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
      } catch (...) {
        cout << "Error during callback" << endl;
      }
}





















callbackimSearch::callbackimSearch(vector<vector<GRBVar>> const & xinvar,
				vector<vector<GRBVar>> const & xoutvar,
				vector<vector<GRBVar>> const & xkeyvar):
	invar(xinvar),
	outvar(xoutvar),
	keyvar(xkeyvar),
	allsolution(),
	count(0)
{}

void callbackimSearch::callback()
{

	
	try{
	
		if(where == GRB_CB_MIPSOL)
		{
			count++;

			cout << endl << "start#############################" << endl;
			uint blockSize = 64;
			
			GRBEnv gurobiEnv;
			gurobiEnv.set(GRB_IntParam_OutputFlag, 0);

			GRBModel m(gurobiEnv, "./modeltwine_SSB/2_twine.mps");


			vector<uint8_t> keysolution (blockSize);
			vector<uint8_t> insolution (blockSize);
			vector<uint8_t> outsolution (blockSize);
			for(uint i = 0; i < blockSize; i++)
			{
				keysolution[i] = uint(round(getSolution(keyvar[0][i])));
				insolution[i] = uint(round(getSolution(invar[0][i])));
				outsolution[i] = uint(round(getSolution(outvar[0][i])));
			}
			cout << "shuru x: " << endl;
			for(uint r =0; r < invar.size(); r++)
			{
				for(uint i = 0; i < blockSize; i++)
					cout << uint(round(getSolution(invar[r][i])));
				cout << endl;
			}
			cout << endl;
			cout << "shuchu v: " << endl;
			for(uint r =0; r < invar.size(); r++)
			{
				for(uint i = 0; i < blockSize; i++)
					cout << uint(round(getSolution(outvar[r][i])));
				cout << endl;
			}
			cout << endl;
			
			cout << "key k : " << endl;
			vector<vector<uint8_t>> allkey(int(invar.size()), vector<uint8_t>(blockSize));
			for(uint r =0; r < invar.size(); r++)
			{
				for(uint i = 0; i < blockSize; i++)
				{
					allkey[r][i] = uint(round(getSolution(keyvar[r][i])));
					cout << int(allkey[r][i]);
				}
				cout << endl;
			}
			cout << endl;
			cout << "count :" << count << endl;
			//getchar();


			
			for(uint i = 0; i < blockSize; i++)
			{
				m.addConstr(m.getVarByName("x_0_"+to_string(i)) == uint(round(getSolution(invar[0][i]))));
				m.addConstr(m.getVarByName("v_1_"+to_string(i)) == uint(round(getSolution(outvar[0][i]))));
				m.addConstr(m.getVarByName("k_0_"+to_string(i)) == uint(round(getSolution(keyvar[0][i]))));
			}

			m.set(GRB_IntParam_PoolSearchMode, 2);
			m.set(GRB_IntParam_PoolSolutions, 200000000);
			m.set(GRB_IntParam_SolutionLimit,100000);
			m.set(GRB_DoubleParam_TimeLimit, 300);

			GRBLinExpr objExpr(0);
			for(uint i = 0; i < blockSize; i++)
			{
				objExpr += m.getVarByName("x_0_"+to_string(i));
			}
			m.setObjective(objExpr);

			m.set(GRB_IntParam_Threads,8);
			m.update();
			m.optimize();

			bool oddtra = true;

			if(m.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
			{
				uint nbTrail = m.get(GRB_IntAttr_SolCount);
				cout << endl << nbTrail << "trails " << endl;

				if(nbTrail % 2== 0)
				{
					oddtra = false;
					GRBLinExpr cutExpr(0);
					for(uint i = 0; i < blockSize; i++)
					{
						if(keysolution[i] == 0) cutExpr += keyvar[0][i];
						else cutExpr += (1-keyvar[0][i]);
						if(insolution[i] == 0) cutExpr += invar[0][i];
						else cutExpr += (1-invar[0][i]);
						if(outsolution[i] == 0) cutExpr += outvar[0][i];
						else cutExpr += (1-outvar[0][i]);
					}
					addLazy( cutExpr >= 1);
				}
				else
				{
					cout << "odd" << endl;

					for(uint i = 0; i < blockSize; i++)
					{
						
						cout << int(keysolution[i]);
						//cout << int(allkey[0][i]);
					}
					cout << endl;
						//abort();
					//allsolution.insert(make_pair(keysolution, insolution));
					//GRBLinExpr cutExpr(0);
					//	for(uint i = 0; i < blockSize; i++)
					//	{
					//		if(allkey[0][i] == 0) cutExpr += keyvar[0][i];
					//		else cutExpr += (1-keyvar[0][i]);
					//		if(insolution[i] == 0) cutExpr += invar[0][i];
					//		else cutExpr += (1-invar[0][i]);
					//		if(outsolution[i] == 0) cutExpr += outvar[0][i];
					//		else cutExpr += (1-outvar[0][i]);
					//	}
					//addLazy( cutExpr >= 1);
					/*
					//abort();
					for(uint i = 0; i < blockSize; i++)
					{
						
						//cout << int(keysolution[i]);
						cout << int(allkey[0][i]);
					}
					cout << endl;
					*/
				}
			}
			else
			{
				oddtra = false;
				cout << "reason:  " << m.get(GRB_IntAttr_Status) << endl;
			}
			
			
			if(oddtra)
			{
				auto lastround = invar.size();
				printf("%ld round \n", lastround);

			
				
				GRBModel m(gurobiEnv, "./modeltwine_SSB/"+to_string(lastround + 1) +"_twine.mps");

				
			        for(uint i = 0; i < blockSize; i++)
				{
					m.addConstr(m.getVarByName("x_0_"+to_string(i)) == int(round(getSolution(invar[0][i]))));
					m.addConstr(m.getVarByName("v_"+to_string(lastround)+"_"+to_string(i)) == int(round(getSolution(outvar[lastround-1][i]))));
				}

			
				for(uint r = 0; r < lastround; r++)
				{
					for(uint i = 0; i < blockSize; i++)
					{
						m.addConstr(m.getVarByName("k_"+to_string(r)+"_"+to_string(i)) == int(round(getSolution(keyvar[r][i]))));
					}
				}

				m.set(GRB_IntParam_PoolSearchMode, 2);
				m.set(GRB_IntParam_PoolSolutions, 200000000);
				m.set(GRB_IntParam_SolutionLimit, 100000);
				m.set(GRB_DoubleParam_TimeLimit, 300);

				GRBLinExpr objExpr(0);
				for(uint i = 0; i < blockSize; i++)
				{
					objExpr += m.getVarByName("x_0_"+to_string(i));
				}
				m.setObjective(objExpr);

				m.set(GRB_IntParam_Threads,8);
				m.update();
				m.optimize();
				if(m.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
				{
					uint lastnbTrail = m.get(GRB_IntAttr_SolCount);
					cout << endl << lastnbTrail << "trails " << endl;

					if(lastnbTrail % 2== 0)
						oddtra = false;
					else
					{
						cout << "Ok,last round is odd" << endl;
						//abort();
						allsolution.insert(make_pair(keysolution, insolution));
						GRBLinExpr cutExpr(0);
							for(uint i = 0; i < blockSize; i++)
							{
								if(allkey[0][i] == 0) cutExpr += keyvar[0][i];
								else cutExpr += (1-keyvar[0][i]);
								if(insolution[i] == 0) cutExpr += invar[0][i];
								else cutExpr += (1-invar[0][i]);
								if(outsolution[i] == 0) cutExpr += outvar[0][i];
								else cutExpr += (1-outvar[0][i]);
							}
						addLazy( cutExpr >= 1);
					}
				
				}
				else
					oddtra = false;
			}
			
			if(!oddtra)
			{
				GRBLinExpr cutExpr(0);
					for(uint i = 0; i < blockSize; i++)
					{
						if(allkey[0][i] == 0) cutExpr += keyvar[0][i];
						else cutExpr += (1-keyvar[0][i]);
						if(insolution[i] == 0) cutExpr += invar[0][i];
						else cutExpr += (1-invar[0][i]);
						if(outsolution[i] == 0) cutExpr += outvar[0][i];
						else cutExpr += (1-outvar[0][i]);
					}
				addLazy( cutExpr >= 1);
			}

			
		}
	}catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl; 
	}catch(...) {
	cout << "Error during callback" << endl;
	}
}

callbackfirst::callbackfirst(vector<GRBVar> const & xinvar,
                             vector<GRBVar> const & xoutvar,
                             vector<GRBVar> const & xkeyvar) :
    invar(xinvar),
    outvar(xoutvar),
    keyvar(xkeyvar)
{}

void callbackfirst::callback()
{
    try{
        if(where == GRB_CB_MIPSOL)
        {
            uint blockSize = 64;
            //ctrsolution++;

            
            GRBEnv gurobiEnv;
            gurobiEnv.set(GRB_IntParam_OutputFlag, 0);

            GRBModel m(gurobiEnv, "first.mps");
	    m.set(GRB_IntParam_Threads,16);

            for(uint i = 0; i < 64; i++)
            {
                    m.addConstr(m.getVarByName("x_0_"+to_string(i)) == uint(round(getSolution(invar[i]))));
                    m.addConstr(m.getVarByName("x_1_"+to_string(i)) == uint(round(getSolution(outvar[i]))));
                    m.addConstr(m.getVarByName("k_0_"+to_string(i)) == uint(round(getSolution(keyvar[i]))));
            }
                

    
            cout << endl << "counting the actual number of trails..." << endl;
	    
	    vector<uint8_t> keysolution (64);
	    vector<uint8_t> insolution (64);
            vector<uint8_t> outsolution (64);
	    for(uint i = 0; i < blockSize; i++)
	    {
    		keysolution[i] = uint(round(getSolution(keyvar[i])));
		insolution[i] = uint(round(getSolution(invar[i])));
		outsolution[i] = uint(round(getSolution(outvar[i])));
	    }
   
            GRBLinExpr objExpr(0);
            for(uint i = 0; i < blockSize; i++)
	    {
			objExpr += m.getVarByName("x_0_"+to_string(i));
	    }
	    m.setObjective(objExpr);

	    for(uint i = 0; i < 64; i++)
	    {
		    cout << int(insolution[i]);
	    }
	    cout << endl;
				
	    for(uint i = 0; i < 64; i++)
	    {
		    cout << int(keysolution[i]);
	    }
	    cout << endl;

	    for(uint i = 0; i < 64; i++)
	    {
		    cout << int(outsolution[i]);
	    }
	    cout << endl;


	    m.set(GRB_IntParam_PoolSearchMode, 2);
      	    m.set(GRB_IntParam_PoolSolutions, 200000000);
	    m.set(GRB_IntParam_SolutionLimit,10);
	    m.set(GRB_DoubleParam_TimeLimit, 60);


	    m.update();			
	    m.optimize();
            
            bool oddSSB = true;
            if(m.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
            {
                uint nbTrail = m.get(GRB_IntAttr_SolCount);
                cout << endl << nbTrail << "trails" << endl;

                if(nbTrail % 2 == 0) oddSSB = false;
                else
                {
                    cout << endl << "odd solution" << oddSSB << endl;
		}
	    }
			
	    else
            {
		        cout << m.get(GRB_IntAttr_Status) << endl;
			cout << endl << "Too many solutions, aborted" << endl;
			oddSSB = false;
     	    }


            if(!oddSSB)
            {
                GRBLinExpr cutExpr(0);
                for(uint i = 0; i < 64; i++)
                {
                    if(keysolution[i] == 0) cutExpr += keyvar[i];
                    else cutExpr += (1-keyvar[i]);
                }
                for(uint i = 0; i < 64; i++)
                {
                    if(outsolution[i] == 0) cutExpr += outvar[i];
                    else cutExpr += (1-outvar[i]);
                }
                
                addLazy( cutExpr >= 1);
            }
        }
    } catch (GRBException e) {
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Error during callback" << endl;
    }
}



callbacklast::callbacklast(vector<GRBVar> const & xinvar,
                             vector<GRBVar> const & xoutvar,
                             vector<GRBVar> const & xkeyvar) :
    invar(xinvar),
    outvar(xoutvar),
    keyvar(xkeyvar)
{}

void callbacklast::callback()
{
    try{
        if(where == GRB_CB_MIPSOL)
        {
            uint blockSize = 64;
            //ctrsolution++;

            
            GRBEnv gurobiEnv;
            gurobiEnv.set(GRB_IntParam_OutputFlag, 0);

            GRBModel m(gurobiEnv, "last.mps");
	    m.set(GRB_IntParam_Threads,16);

            for(uint i = 0; i < 64; i++)
            {
                    m.addConstr(m.getVarByName("x_0_"+to_string(i)) == uint(round(getSolution(invar[i]))));
                    m.addConstr(m.getVarByName("x_1_"+to_string(i)) == uint(round(getSolution(outvar[i]))));
                    m.addConstr(m.getVarByName("k_0_"+to_string(i)) == uint(round(getSolution(keyvar[i]))));
            }
                

            cout << endl << "counting the actual number of trails..." << endl;
	    
	    vector<uint8_t> keysolution (64);
	    vector<uint8_t> insolution (64);
            vector<uint8_t> outsolution (64);
	    for(uint i = 0; i < blockSize; i++)
	    {
    		keysolution[i] = uint(round(getSolution(keyvar[i])));
		insolution[i] = uint(round(getSolution(invar[i])));
		outsolution[i] = uint(round(getSolution(outvar[i])));
	    }
            
            GRBLinExpr objExpr(0);
            for(uint i = 0; i < blockSize; i++)
	    {
			objExpr += m.getVarByName("x_1_"+to_string(i));
	    }
	    m.setObjective(objExpr);

	    for(uint i = 0; i < 64; i++)
	    {
		    cout << int(insolution[i]);
	    }
	    cout << endl;
				
	    for(uint i = 0; i < 64; i++)
	    {
		    cout << int(keysolution[i]);
	    }
	    cout << endl;

	    for(uint i = 0; i < 64; i++)
	    {
		    cout << int(outsolution[i]);
	    }
	    cout << endl;

	    
	    m.set(GRB_IntParam_PoolSearchMode, 2);
      	    m.set(GRB_IntParam_PoolSolutions, 200000000);
	    m.set(GRB_IntParam_SolutionLimit,10);
	    m.set(GRB_DoubleParam_TimeLimit, 60);


	    m.update();			
	    m.optimize();
            
            bool oddSSB = true;
            if(m.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
            {
                uint nbTrail = m.get(GRB_IntAttr_SolCount);
                cout << endl << nbTrail << "trails" << endl;

                if(nbTrail % 2 == 0) oddSSB = false;
                else
                {
                    cout << endl << "odd solution" << oddSSB << endl;
		}
	    }
			
	    else
            {
		        cout << "reason: " << m.get(GRB_IntAttr_Status) << endl;
			cout << endl << "Too many solutions, aborted" << endl;
			oddSSB = false;
     	    }


            if(!oddSSB)
            {
                GRBLinExpr cutExpr(0);
                for(uint i = 0; i < 64; i++)
                {
                    if(keysolution[i] == 0) cutExpr += keyvar[i];
                    else cutExpr += (1-keyvar[i]);
                }
                for(uint i = 0; i < 64; i++)
                {
                    if(outsolution[i] == 0) cutExpr += outvar[i];
                    else cutExpr += (1-outvar[i]);
                }
                
                addLazy( cutExpr >= 1);
            }
        }
    } catch (GRBException e) {
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Error during callback" << endl;
    }
}


callbacksearchKeyall::callbacksearchKeyall(uint xrMax,
			     vector<GRBVar> const & xinvar,
                             vector<GRBVar> const & xoutvar,
                             vector<vector<GRBVar>> const & xkeyvar) :
    rMax(xrMax),
    invar(xinvar),
    outvar(xoutvar),
    keyvar(xkeyvar)
{}

void callbacksearchKeyall::callback()
{
    try{
        if(where == GRB_CB_MIPSOL)
        {
            uint blockSize = 64;
            //ctrsolution++;

            
            GRBEnv gurobiEnv;
            gurobiEnv.set(GRB_IntParam_OutputFlag, 0);

            GRBModel m(gurobiEnv, "./modeltwine_SSB/"+to_string(rMax)+"_twine.mps");
	    m.set(GRB_IntParam_Threads,16);

            for(uint i = 0; i < 64; i++)
            {
                    m.addConstr(m.getVarByName("x_0_"+to_string(i)) == uint(round(getSolution(invar[i]))));
                    m.addConstr(m.getVarByName("v_"+to_string(rMax-1)+"_"+to_string(i)) == uint(round(getSolution(outvar[i]))));
            }
	    for(uint r = 0; r < rMax-1; r++)
	    {
		    for(uint i = 0; i < blockSize; i++)
		    {
			    m.addConstr(m.getVarByName("k_"+to_string(r)+"_"+to_string(i)) == uint(round(getSolution(keyvar[r][i]))));
		    }
	    }
                
            
            cout << endl << "counting the actual number of trails..." << endl;
	    
	    vector<vector<uint8_t>> keysolution(rMax - 1, vector<uint8_t> (64));
	    vector<uint8_t> insolution (64);
            vector<uint8_t> outsolution (64);
	    for(uint i = 0; i < blockSize; i++)
	    {
		insolution[i] = uint(round(getSolution(invar[i])));
		outsolution[i] = uint(round(getSolution(outvar[i])));
	    }
	    for(uint r = 0; r < rMax - 1; r++)
	    {
		    for(uint i = 0; i < blockSize; i++)
		    {
			    keysolution[r][i] = uint(round(getSolution(keyvar[r][i])));
		    }
	    }
            
	   

            GRBLinExpr objExpr(0);
            for(uint i = 0; i < blockSize; i++)
	    {
			objExpr += m.getVarByName("x_0_"+to_string(i));
	    }
	    m.setObjective(objExpr);

	    
	    m.set(GRB_IntParam_PoolSearchMode, 2);
      	    m.set(GRB_IntParam_PoolSolutions, 200000000);
	    m.set(GRB_IntParam_SolutionLimit,200000);
	    m.set(GRB_DoubleParam_TimeLimit, 500);


	    m.update();			
	    m.optimize();
            
            bool oddSSB = true;
            if(m.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
            {
                uint nbTrail = m.get(GRB_IntAttr_SolCount);
                cout << endl << nbTrail << "trails" << endl;

                if(nbTrail % 2 == 0) oddSSB = false;
                else
                {
                    cout << endl << "odd solution" << oddSSB << endl;
		}
	    }
			
	    else
            {
		        cout << "no solution reason: " << m.get(GRB_IntAttr_Status) << endl;
			cout << endl << "Too many solutions, aborted" << endl;
			oddSSB = false;
     	    }

	    //uint indexin;
	    //for(uint i = 0; i < 64; i++)
	    //{
	    //    if(insolution[i] == 0)
	    //		    indexin = i;
	    //}

	    //cout << "input index:i  " << indexin << endl;
	/*
	    if(oddSSB)
	    {
		    vector<uint8_t> insolution16(64, 1);
		    vector<uint8_t> insolution32(64, 1);
		    insolution16[indexin - 16] = 0;
	    	    cout << "input index:i  " << indexin << endl;

		    insolution32[indexin - 32] = 0;
		    auto result_16 = counttrail_middle(rMax, insolution16, outsolution, keysolution);
		    if(result_16.second == true && result_16.first % 2 == 1)
			    oddSSB = false;
		    auto result_32 = counttrail_middle(rMax, insolution32, outsolution, keysolution);
		    if(result_32.second == true && result_32.first % 2 == 1)
			    oddSSB = false;
	    }
	*/






	    
            if(!oddSSB)
            {
                GRBLinExpr cutExpr(0);
		for(uint r = 0; r < rMax -1; r++)
		{
			for(uint i = 0; i < blockSize; i++)
			{
				if(keysolution[r][i] == 0) cutExpr += keyvar[r][i];
				else cutExpr += (1 - keyvar[r][i]);
			}
		}
                addLazy( cutExpr >= 1);
            }
        }
    } catch (GRBException e) {
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Error during callback" << endl;
    }
}

/*	
callbackcount::callbackcount(vector<GRBVar> const & xinvar,
                             vector<GRBVar> const & xoutvar);
    invar(xinvar),
    outvar(xoutvar),
    keyvar(xkeyvar)
{}

void callbacklast::callback()
{
    try{
        if(where == GRB_CB_MIPSOL)
        {
            uint blockSize = 64;
            //ctrsolution++;

            
            GRBEnv gurobiEnv;
            gurobiEnv.set(GRB_IntParam_OutputFlag, 0);

            GRBModel m(gurobiEnv, "last.mps");
	    m.set(GRB_IntParam_Threads,16);

            for(uint i = 0; i < 64; i++)
            {
                    m.addConstr(m.getVarByName("x_0_"+to_string(i)) == uint(round(getSolution(invar[i]))));
                    m.addConstr(m.getVarByName("x_1_"+to_string(i)) == uint(round(getSolution(outvar[i]))));
                    m.addConstr(m.getVarByName("k_0_"+to_string(i)) == uint(round(getSolution(keyvar[i]))));
            }
                

            
            cout << endl << "counting the actual number of trails..." << endl;
	    
	    vector<uint8_t> keysolution (64);
	    vector<uint8_t> insolution (64);
            vector<uint8_t> outsolution (64);
	    for(uint i = 0; i < blockSize; i++)
	    {
    		keysolution[i] = uint(round(getSolution(keyvar[i])));
		insolution[i] = uint(round(getSolution(invar[i])));
		outsolution[i] = uint(round(getSolution(outvar[i])));
	    }
            
            GRBLinExpr objExpr(0);
            for(uint i = 0; i < blockSize; i++)
	    {
			objExpr += m.getVarByName("x_1_"+to_string(i));
	    }
	    m.setObjective(objExpr);

	    for(uint i = 0; i < 64; i++)
	    {
		    cout << int(insolution[i]);
	    }
	    cout << endl;
				
	    for(uint i = 0; i < 64; i++)
	    {
		    cout << int(keysolution[i]);
	    }
	    cout << endl;

	    for(uint i = 0; i < 64; i++)
	    {
		    cout << int(outsolution[i]);
	    }
	    cout << endl;

	    
	    m.set(GRB_IntParam_PoolSearchMode, 2);
      	    m.set(GRB_IntParam_PoolSolutions, 200000000);
	    m.set(GRB_IntParam_SolutionLimit,10);
	    m.set(GRB_DoubleParam_TimeLimit, 60);


	    m.update();			
	    m.optimize();
            
            bool oddSSB = true;
            if(m.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
            {
                uint nbTrail = m.get(GRB_IntAttr_SolCount);
                cout << endl << nbTrail << "trails" << endl;

                if(nbTrail % 2 == 0) oddSSB = false;
                else
                {
                    cout << endl << "odd solution" << oddSSB << endl;
		}
	    }
			
	    else
            {
		        cout << "reason: " << m.get(GRB_IntAttr_Status) << endl;
			cout << endl << "Too many solutions, aborted" << endl;
			oddSSB = false;
     	    }


            if(!oddSSB)
            {
                GRBLinExpr cutExpr(0);
                for(uint i = 0; i < 64; i++)
                {
                    if(keysolution[i] == 0) cutExpr += keyvar[i];
                    else cutExpr += (1-keyvar[i]);
                }
                for(uint i = 0; i < 64; i++)
                {
                    if(outsolution[i] == 0) cutExpr += outvar[i];
                    else cutExpr += (1-outvar[i]);
                }
                
                addLazy( cutExpr >= 1);
            }
        }
    } catch (GRBException e) {
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Error during callback" << endl;
    }
}

*/
