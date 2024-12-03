#ifndef H_CALLBACK                                                                                                                                                   
#define H_CALLBACK

#include<vector>
#include<iostream>
#include<cmath>
#include<cstdint>

#include<boost/dynamic_bitset.hpp>

#include "gurobi_c++.h"
#include "SearchKey.hpp"

class callbackSearchKeyPattern: public GRBCallback
{
    public:
        unsigned int nbRounds;
        std::vector<GRBVar> inputvar;
        std::vector<GRBVar> outputvar;
        std::vector<std::vector<GRBVar>> allkeyvar;
        //uint64_t ctrsolution;
        //std::vector<std::vector<boost::dynamic_bitset<>>> allSolution;

        callbackSearchKeyPattern(unsigned int const xnbRounds, std::vector<GRBVar> const & xinputvar, std::vector<GRBVar> const & xoutputvar, std::vector<std::vector<GRBVar>> const & xallkeyvar);

    protected:
        void callback();
};







class callbackimSearch: public GRBCallback
{
    public:
        std::vector<std::vector<GRBVar>> invar;
        std::vector<std::vector<GRBVar>> outvar;
        std::vector<std::vector<GRBVar>> keyvar;
	std::set<std::pair<std::vector<uint8_t>, std::vector<uint8_t>>> allsolution;
	int count;

        callbackimSearch( std::vector<std::vector<GRBVar>> const & xinvar, std::vector<std::vector<GRBVar>> const & xoutvar, std::vector<std::vector<GRBVar>> const & xkeyvar);

    protected:
        void callback();
};

class callbackfirst: public GRBCallback
{
    public:
        std::vector<GRBVar> invar;
        std::vector<GRBVar> outvar;
        std::vector<GRBVar> keyvar;
        //uint64_t ctrsolution;
        //std::vector<std::vector<boost::dynamic_bitset<>>> allSolution;

        callbackfirst( std::vector<GRBVar> const & xinvar, std::vector<GRBVar> const & xoutvar, std::vector<GRBVar> const & xkeyvar);

    protected:
        void callback();
};

class callbacklast: public GRBCallback
{
    public:
        std::vector<GRBVar> invar;
        std::vector<GRBVar> outvar;
        std::vector<GRBVar> keyvar;
        //uint64_t ctrsolution;
        //std::vector<std::vector<boost::dynamic_bitset<>>> allSolution;

        callbacklast( std::vector<GRBVar> const & xinvar, std::vector<GRBVar> const & xoutvar, std::vector<GRBVar> const & xkeyvar);

    protected:
        void callback();
};

class callbacksearchKeyall: public GRBCallback

{
    public:
	uint rMax;
        std::vector<GRBVar> invar;
        std::vector<GRBVar> outvar;
	std::vector<std::vector<GRBVar>> keyvar;
        //uint64_t ctrsolution;
        //std::vector<std::vector<boost::dynamic_bitset<>>> allSolution;

        callbacksearchKeyall( uint const xrMax, std::vector<GRBVar> const & xinvar, std::vector<GRBVar> const & xoutvar, std::vector<std::vector<GRBVar>> const & xkeyvar);

    protected:
        void callback();
};

//callback class for the improved dynamic search
class callbackDynamic: public GRBCallback
{
  public:
    uint rMiddle;
    std::vector<uint8_t> const & output;
    std::vector<std::vector<uint8_t>> const & keyVal;
    std::vector<GRBVar> const & inputMiddleVar;
    std::vector<GRBVar> const & keyMiddleVar;
    std::vector<std::vector<std::vector<GRBVar>>> const & inSSBVar;         //input variables for the SSB
    std::vector<std::vector<std::vector<GRBVar>>> const & outSSBVar;        //output variables for the SSB

    callbackDynamic(uint const xrMiddle,
                    std::vector<uint8_t> const & xoutput,
                    std::vector<std::vector<uint8_t>> const & xkeyVal,
                    std::vector<GRBVar> const & xinputMiddleVar,
                    std::vector<GRBVar> const & xkeyMiddleVar,
                    std::vector<std::vector<std::vector<GRBVar>>> const & xinSSBVar,
                    std::vector<std::vector<std::vector<GRBVar>>> const & xoutSSBVar);

  protected:
    void callback();
};

class callbacksearchkey: public GRBCallback
{
  public:

    unsigned int nbRounds;                                          //number of rounds convered by the model
    std::vector<GRBVar> inputVar;                                   //variables for the input
    std::vector<GRBVar> outputVar;                                  //variables for the output
    std::vector<std::vector<GRBVar>> allKeyVar;                     //all key variables
    std::vector<std::vector<std::vector<GRBVar>>> inSSBVar;         //input variables for the SSB
    std::vector<std::vector<std::vector<GRBVar>>> outSSBVar;        //output variables for the SSB
    uint64_t ctrsolutions;                                          //counter for the solutions found
    std::vector<std::vector<boost::dynamic_bitset<>>> allSolutions; //vector to keep all examined solutions

    callbacksearchkey(unsigned int const xnbRounds,
                             std::vector<GRBVar> const & xinputVar,
                             std::vector<GRBVar> const & xoutputVar,
                             std::vector<std::vector<GRBVar>> const & xallKeyVar,
                             std::vector<std::vector<std::vector<GRBVar>>> const & xinSSBVar,
                             std::vector<std::vector<std::vector<GRBVar>>> const & xoutSSBVar);

  protected:
    void callback();
};















#endif

