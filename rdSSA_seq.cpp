
/********************		Gene expression with delay		**************************/

// g++ -std=c++17 -Wall -g -I ~/fastflow -pthread -DNDEBUG -DTRACE_FASTFLOW -finline-functions -O3 -fopenmp frm_ff_delay.cpp -o ff

#include <iostream>
#include <math.h>
#include <limits.h>
#include <vector>
#include <string>
#include <omp.h>
#include <chrono>
#include <utility> 
#include <time.h> 
#include <mutex>
#include <random>
#include <time.h>
#include <algorithm>
#include <stdlib.h>
#include <thread>
#include <atomic>
#include <pthread.h>
#include <omp.h>


#include <ff/ff.hpp>
#include <ff/pipeline.hpp>
#include <ff/parallel_for.hpp>

#include "utimer.cpp"

using namespace std;
using namespace ff;

#define ENDTIME   2			// end of time
#define M 		  3				// number of reaction (should be divisible by 3)
#define N		  4				// number of species (should be divisible by 4)
#define tDelay	  0.1			// delay time for reactions NCD, CD

int num_simulation;
int num_thread_rand;
int num_worker;

double 	rate[N];		    // reaction rates
int		withDelay[N];		//if reaaction is CD = 1 , NCD = 2, ND = 3 

vector<int>	dependencyGraph[M];
vector<pair<int, int>>	reactant[N];
vector<pair<int, int>>	product[N];
	
int 	x[M];	 			// population of chemical species		
int 	x_lower_bound[M];
int		x_upper_bound[M];
double 	p[N];		    	// propensities of reactants
double 	p_lower_bound[N];
double	p_upper_bound[N];
double	total_propensity_upper_bound;
double  rand_seed;


queue<pair<double, int>> q; 
queue<double> randomNumberList; 

double 	t =0.0;
double 	tau;
double	rand1, rand2, rand3;


int idOfElement;
int changedAmount;
double td ;		
int reac_id;
std::pair<double, int> reactionFromQueue;
int minIndex;
bool accepted;


double GenerateSeed()
{
	auto thread_id = std::this_thread::get_id();
	uint64_t* ptr=(uint64_t*) &thread_id;
	double seed = *ptr + time(NULL);		
	srand(seed);
	
	double my_rand = (double)rand();
	return my_rand;
}
void init()
{	
	// reaction rates	
	rate[0] = 1;		 	withDelay[0] = 2; //"NCD";
	rate[1] = 0.01;		 	withDelay[1] = 3; //"ND";
	rate[2] = 5;		 	withDelay[2] = 3; //"ND";
	rate[3] = 0.01;		 	withDelay[3] = 3; //"ND";
	
	
	// S1 --r0--> S1 + S2	
	reactant[0].push_back(make_pair(0, 1)); 	//r0(s1, 1)
	product[0].push_back(make_pair(0, 1));	 	//r0(s1, 1)
	product[0].push_back(make_pair(1, 1));		//r0(s2, 1)		
	
	// S1 + S2 --r1--> S3		
	reactant[1].push_back(make_pair(0, 1)); 	//r1(s1, 1)
	reactant[1].push_back(make_pair(1, 1));	 	//r1(s2, 1)
	product[1].push_back(make_pair(2, 1));		//r1(s3, 1)
	
	// S3 --r2--> S1 + S2		
	reactant[2].push_back(make_pair(2, 1)); 	//r2(s3, 1)
	product[2].push_back(make_pair(0, 1));	 	//r2(s1, 1)
	product[2].push_back(make_pair(1, 1));		//r2(s2, 1)
			
	// S2 --r3--> _		
	reactant[3].push_back(make_pair(1, 1)); 	//r3(s2, 1)		
}

void SR_graph()		// Species-Reaction graph
{
	// S1 --r0--> S1 + S2	
	// S1 + S2 --r1--> S3
	// S3 --r2--> S1 + S2
	// S2 --r3--> _
	dependencyGraph[0].push_back(0);	//G(S1,r0), G(r0, S1)
	dependencyGraph[0].push_back(1);	//G(S1,r1)
	dependencyGraph[0].push_back(2);	//G(r2,S1)
	
	dependencyGraph[1].push_back(1);	//G(S2,r1)
	dependencyGraph[1].push_back(3);	//G(S2,r3)
	dependencyGraph[1].push_back(3);	//G(r2,S2)
	
	dependencyGraph[2].push_back(2);	//G(S3,r2)
	dependencyGraph[2].push_back(1);	//G(r1,S3)
}
	
bool checkSpeciesBound()
{
	for(int i=0; i<M; i++)
	{
		if( x[i]< x_lower_bound[i] ||  x[i]> x_upper_bound[i])
		{
			return false;
		}
	}	
	return true;
};	

double propensity_sum()
{
	 p[0] = rate[0]* x[0];		    		// S1 --r0--> S1 + S2
	 p[1] = rate[1]* x[0]* x[1];				// S1 + S2 --r1--> S3
	 p[2] = rate[2]* x[2];		    		// S3 --r2--> S1 + S2
	 p[3] = rate[3]* x[1];		    		// S2 --r3--> _		    

	double p_sum=0.0;
	
	for(int i=0; i<N; i++)
		p_sum +=  p[i];
	
	return p_sum;
};

void computeSpeciesBounds(int i)
{
	 x_lower_bound[i] =  x[i] - 0.1* x[i];
	if( x_lower_bound[i]<0)
		 x_lower_bound[i]=0;
	 x_upper_bound[i] =  x[i] + 0.1* x[i];
};

void computePropensity(int i)
{
	switch(i)
	{
		case 0:
				{	 p[0] = rate[0]* x[0]; break; }		// S1 --r0--> S1 + S2
		case 1:
				{	 p[1] = rate[1]* x[0]* x[1]; break; }	// S1 + S2 --r1--> S3
		case 2:
				{	 p[2] = rate[2]* x[2]; break; }		// S3 --r2--> S1 + S2
		case 3:
				{	 p[3] = rate[3]* x[1]; break; }		// S2 --r3--> _	
	}				    		
};

void computePropensityBounds(int i)
{
	switch(i)
	{
		case 0:	 
				{	 p_lower_bound[0] = rate[0]* x_lower_bound[0]; 
					 p_upper_bound[0] = rate[0]* x_upper_bound[0];
					break;
				}		
		case 1:
				{	
					 p_lower_bound[1] = rate[1]* x_lower_bound[0]* x_lower_bound[1]; 
					 p_upper_bound[1] = rate[1]* x_upper_bound[0]* x_upper_bound[1];
					break; 
				}
		case 2:
				{	
					 p_lower_bound[2] = rate[2]* x_lower_bound[2];
					 p_upper_bound[2] = rate[2]* x_upper_bound[2];
					break; 
				}		
		case 3:
				{
					 p_lower_bound[3] = rate[3]* x_lower_bound[1];
					 p_upper_bound[3] = rate[3]* x_upper_bound[1];
					break; 
				}
	}
};

double propensityUpperBoundSum()
{
	double p_sum = 0;
	for(int i=0; i<N; i++)
		p_sum +=  p_upper_bound[i];
	
	return p_sum;
}

double randDouble(double seed) 
{
	//thread_local static std::random_device rd;
	thread_local static std::mt19937 rng(seed);
	thread_local std::uniform_real_distribution<double> urd;
	return urd(rng, decltype(urd)::param_type{0,1});
};

int selectMinIndex(double rand2)
{
	double sum=0;
	int i=-1;
	while(sum<=rand2* total_propensity_upper_bound)
	{
		i++;
		sum +=  p_upper_bound[i];		
	}
	
	return i;
};

void delayedRejectionBasedSelection()
{
	//cout<<"inside method"<<endl;
	/************** Generate random number tau ************/		
	 rand1 =  randDouble(rand_seed);
	//cout<<"rand: "<<rand1<<endl;
	 tau = (-1/ total_propensity_upper_bound)*log( rand1);
	//cout<<"tau: "<<tau<<endl;
			
	/************** Checking reactions on the waiting queue ************/
	while(! q.empty())
	{
					
		 reactionFromQueue =  q.front();
		
		 td 		=  reactionFromQueue.first;
		 reac_id  	=  reactionFromQueue.second; 
		//cout<<"Queue is not empty! t=" << t<<endl;
		if( td <  t+ tau)
		{				
			//update propensities of previousReaction
			switch(withDelay[ reac_id])
			{
				case 1: //CD	
							{
								// update the abundance of products
								//cout<<"from queue CD R"<<reac_id<<"fired"<<endl;
								for(int i=0; i<(int)product[ reac_id].size(); i++)
								{
									 idOfElement   = product[ reac_id][i].first;
									 changedAmount = product[ reac_id][i].second;
									 x[ idOfElement] +=  changedAmount;
								}
								break;
							}
				case  2: //NCD	
							{
								// update the abundance of reactants and products
								//cout<<"from queue NCD R"<<reac_id<<"fired"<<endl;
								for(int i=0;i<(int)reactant[ reac_id].size(); i++)
								{
									 idOfElement   = reactant[ reac_id][i].first;
									 changedAmount = reactant[ reac_id][i].second;
									 x[ idOfElement] -=  changedAmount;
								}
								for(int i=0; i<(int)product[ reac_id].size(); i++)
								{
									 idOfElement   = product[ reac_id][i].first;
									 changedAmount = product[ reac_id][i].second;
									 x[ idOfElement] +=  changedAmount;
								}
								break;
							}
			}
						
			 q.pop();
			
			if(checkSpeciesBound()==false)
			{
				//cout<<"stopping firing from queue before t=" << t<<"  and td="<< td;
				 t =  td;
				//cout<<"but now t=" << t<<endl;
				break;					
			}
		}
		else
		{
			break;
		}
	}
	
	 t =  t +  tau;
	//cout<<"after queue we fire next choosen reaction t= "<< t<<endl;
	//sleep(1);
	 rand2 =  randDouble(rand_seed);
	 rand3 =  randDouble(rand_seed);

	
	 minIndex = selectMinIndex(rand2);
	
	 accepted=false;
	if( rand3<=( p_lower_bound[ minIndex]/ p_upper_bound[ minIndex]))
		 accepted=true;
	else
	{
		computePropensity(minIndex);
		if( rand3<=( p[ minIndex]/ p_upper_bound[ minIndex]))
			 accepted=true;
	}
	
	/************** Firing the chosen reaction ************/
	if( accepted)
	{
		//cout<<"trying to fire new reaction ";
		switch(withDelay[ minIndex])
		{
			case 1: //CD 
						{
						//cout<<"CD R"<<minIndex<<" goes to queue"<<endl;
						for(int i=0;i<(int)reactant[ minIndex].size(); i++)
						{
							 idOfElement   = reactant[ minIndex][i].first;
							 changedAmount = reactant[ minIndex][i].second;
							 x[ idOfElement] -=  changedAmount;
						}							
						//cout<<"CD reaction goes to queue becaus of delay."<<endl;
						 q.push(make_pair( t+tDelay,  minIndex));
						break;
						}
			case 2: //NCD	
						{
						//cout<<"NCD R"<<minIndex<<" goes to queue"<<endl;								
						//cout<<"NCD reaction goes to queue becaus of delay."<<endl;
						 q.push(make_pair( t+tDelay,  minIndex));
						break;
						}
			case 3: //ND	
						{
						//cout<<"ND R"<<minIndex<<endl;
						for(int i=0;i<(int)reactant[ minIndex].size(); i++)
						{
							 idOfElement   = reactant[ minIndex][i].first;
							 changedAmount = reactant[ minIndex][i].second;
							 x[ idOfElement] -=  changedAmount;
						}				
						for(int i=0; i<(int)product[ minIndex].size(); i++)
						{
							 idOfElement   = product[ minIndex][i].first;
							 changedAmount = product[ minIndex][i].second;
							 x[ idOfElement] +=  changedAmount;
						}							
						//cout<<"CD reaction goes to queue becaus of delay."<<endl;
						break;
						}
		}
	}

};

void Simulation_Worker()
{
	vector<int> simulation;
	for(int k=0; k<M; k++)
		simulation.push_back(0);
	
	/************** Starting simulations************/	
	for(int s=0; s<num_simulation; s++)
	{
		//cout<<endl<<"Starting sim: "<<s<<endl;
		 rand_seed = GenerateSeed();
		
		 t=0.0;
		 x[0] = 1000;
		 x[1] = 0;	
		 x[2] = 0;
		 
		while(! q.empty())
			 q.pop();
		
		for(int i=0; i<M; i++)
			computeSpeciesBounds(i);
		
		for(int i=0; i<N; i++)
			computePropensityBounds(i);
			
		 total_propensity_upper_bound = propensityUpperBoundSum();
		
		while( t<ENDTIME)
		{
			while(checkSpeciesBound()==true &&  t<ENDTIME)
			{
				//cout<<"Calling delayedRejectionBasedSelection by sim "<<s<<"  t= "<< t<<endl;
				delayedRejectionBasedSelection();					
			}
		
			/**************Updating the propensities of the Species************/
			//cout<<"Updating the propensities"<<endl;
			for(int i=0; i<M; i++)
			{
				if( x[i]< x_lower_bound[i] ||  x[i]> x_upper_bound[i])
				{
					computeSpeciesBounds(i);
					for(int j=0; j<N; j++)
					{
						for(int k=0; k<(int)dependencyGraph[j].size(); k++)
						{
/*------------*/			double old_value =  p_upper_bound[k];
							computePropensityBounds(k);
							 total_propensity_upper_bound +=   p_upper_bound[k] - old_value;
						}
					}
				}
			}
		}
		//cout<<endl<<"Ending sim: "<<s<<endl;
		for(int k=0; k<M; k++)
			simulation[k] += x[k];
	}
		
	//the average result of the simulations 
	for(int k=0; k<M; k++)
	{
		simulation[k] = simulation[k]/num_simulation;
		cout<<simulation[k]<<" ";
	}
	cout<<endl;
}



int main(int argc, char * argv[])
{
	 if (argc<2) 
	 {
        std::cerr << "use: " << argv[0]  << "num_simulation\n";
        return -1;
    }
	
	
	num_simulation	= std::stol(argv[1]);	// The number of simulations	
	
	init();
	SR_graph();
	
	
//************************** Fast Flow ************************	
	auto startReac = chrono::high_resolution_clock::now();
	Simulation_Worker();
	auto elapsedReac = chrono::high_resolution_clock::now() - startReac;	
	cout << "Time	chrono: " << chrono::duration_cast<chrono::milliseconds>(elapsedReac).count() << " milliseconds"<< endl; 
	//std::cout << "Time	pipe.ffTime(): " << farm.ffTime() << "\n";	
//*****************************************************************
	
	return 0;
}

