
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
	

class Variables
{
	public:
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
};

double GenerateSeedPerThread()
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
	
bool checkSpeciesBound(Variables &v)
{
	for(int i=0; i<M; i++)
	{
		if(v.x[i]<v.x_lower_bound[i] || v.x[i]>v.x_upper_bound[i])
		{
			return false;
		}
	}	
	return true;
};	

double propensity_sum(Variables &v)
{
	v.p[0] = rate[0]*v.x[0];		    		// S1 --r0--> S1 + S2
	v.p[1] = rate[1]*v.x[0]*v.x[1];				// S1 + S2 --r1--> S3
	v.p[2] = rate[2]*v.x[2];		    		// S3 --r2--> S1 + S2
	v.p[3] = rate[3]*v.x[1];		    		// S2 --r3--> _		    

	double p_sum=0.0;
	
	for(int i=0; i<N; i++)
		p_sum += v.p[i];
	
	return p_sum;
};

void computeSpeciesBounds(Variables &v, int i)
{
	v.x_lower_bound[i] = v.x[i] - 0.1*v.x[i];
	if(v.x_lower_bound[i]<0)
		v.x_lower_bound[i]=0;
	v.x_upper_bound[i] = v.x[i] + 0.1*v.x[i];
};

void computePropensity(Variables &v, int i)
{
	switch(i)
	{
		case 0:
				{	v.p[0] = rate[0]*v.x[0]; break; }		// S1 --r0--> S1 + S2
		case 1:
				{	v.p[1] = rate[1]*v.x[0]*v.x[1]; break; }	// S1 + S2 --r1--> S3
		case 2:
				{	v.p[2] = rate[2]*v.x[2]; break; }		// S3 --r2--> S1 + S2
		case 3:
				{	v.p[3] = rate[3]*v.x[1]; break; }		// S2 --r3--> _	
	}				    		
};

void computePropensityBounds(Variables &v, int i)
{
	switch(i)
	{
		case 0:	 
				{	v.p_lower_bound[0] = rate[0]*v.x_lower_bound[0]; 
					v.p_upper_bound[0] = rate[0]*v.x_upper_bound[0];
					break;
				}		
		case 1:
				{	
					v.p_lower_bound[1] = rate[1]*v.x_lower_bound[0]*v.x_lower_bound[1]; 
					v.p_upper_bound[1] = rate[1]*v.x_upper_bound[0]*v.x_upper_bound[1];
					break; 
				}
		case 2:
				{	
					v.p_lower_bound[2] = rate[2]*v.x_lower_bound[2];
					v.p_upper_bound[2] = rate[2]*v.x_upper_bound[2];
					break; 
				}		
		case 3:
				{
					v.p_lower_bound[3] = rate[3]*v.x_lower_bound[1];
					v.p_upper_bound[3] = rate[3]*v.x_upper_bound[1];
					break; 
				}
	}
};

double propensityUpperBoundSum(Variables &v)
{
	double p_sum = 0;
	for(int i=0; i<N; i++)
		p_sum += v.p_upper_bound[i];
	
	return p_sum;
}

double randDouble(double seed) 
{
	//thread_local static std::random_device rd;
	thread_local static std::mt19937 rng(seed);
	thread_local std::uniform_real_distribution<double> urd;
	return urd(rng, decltype(urd)::param_type{0,1});
};

void GenerateRandomNumbers(Variables &v, double seed)
{
	double number[100000];
	//cout<<"GenerateRandomNumbers for sim"<<endl;
	ff::ParallelFor pfor(num_thread_rand);
	pfor.parallel_for(0L, 100000, 1, 0, [&v, &seed, &number](const long i)
	{
		//cout<<"here ";
		//cout<<"pf: "<<i<<" ";
		number[i]=randDouble(seed);
	});
	//cout<<endl;	
	
	for(int i=0; i<500; i++)
		v.randomNumberList.push(number[i]);
};

int selectMinIndex(Variables &v, double rand2)
{
	double sum=0;
	int i=-1;
	while(sum<=rand2*v.total_propensity_upper_bound)
	{
		i++;
		sum += v.p_upper_bound[i];		
	}
	
	return i;
};

void delayedRejectionBasedSelection(Variables &v, double seed, int simId)
{
	//cout<<"inside method"<<endl;
	/************** Generate random number tau ************/
	if(v.randomNumberList.empty())
	{
		//cout<<"GenerateRandomNumbers for sim "<<simId<<"   seed = "<<seed<<endl;
		GenerateRandomNumbers(ref(v), seed);
	}		
	v.rand1 = v.randomNumberList.front();
	v.randomNumberList.pop();
	//cout<<"rand: "<<rand1<<endl;
	v.tau = (-1/v.total_propensity_upper_bound)*log(v.rand1);
	//cout<<"tau: "<<tau<<endl;
			
	/************** Checking reactions on the waiting queue ************/
	while(!v.q.empty())
	{
					
		v.reactionFromQueue = v.q.front();
		
		v.td 		= v.reactionFromQueue.first;
		v.reac_id  	= v.reactionFromQueue.second; 
		//cout<<"Queue is not empty! t=" <<v.t<<endl;
		if(v.td < v.t+v.tau)
		{				
			//update propensities of previousReaction
			switch(withDelay[v.reac_id])
			{
				case 1: //CD	
							{
								// update the abundance of products
								//cout<<"from queue CD R"<<reac_id<<"fired"<<endl;
								for(int i=0; i<(int)product[v.reac_id].size(); i++)
								{
									v.idOfElement   = product[v.reac_id][i].first;
									v.changedAmount = product[v.reac_id][i].second;
									v.x[v.idOfElement] += v.changedAmount;
								}
								break;
							}
				case  2: //NCD	
							{
								// update the abundance of reactants and products
								//cout<<"from queue NCD R"<<reac_id<<"fired"<<endl;
								for(int i=0;i<(int)reactant[v.reac_id].size(); i++)
								{
									v.idOfElement   = reactant[v.reac_id][i].first;
									v.changedAmount = reactant[v.reac_id][i].second;
									v.x[v.idOfElement] -= v.changedAmount;
								}
								for(int i=0; i<(int)product[v.reac_id].size(); i++)
								{
									v.idOfElement   = product[v.reac_id][i].first;
									v.changedAmount = product[v.reac_id][i].second;
									v.x[v.idOfElement] += v.changedAmount;
								}
								break;
							}
			}
						
			v.q.pop();
			
			if(checkSpeciesBound(ref(v))==false)
			{
				//cout<<"stopping firing from queue before t=" <<v.t<<"  and td="<<v.td;
				v.t = v.td;
				//cout<<"but now t=" <<v.t<<endl;
				break;					
			}
		}
		else
		{
			break;
		}
	}
	
	v.t = v.t + v.tau;
	//cout<<"after queue we fire next choosen reaction t= "<<v.t<<endl;
	//sleep(1);
	if(v.randomNumberList.empty())
		GenerateRandomNumbers(ref(v), seed);
	v.rand2 = v.randomNumberList.front();
	v.randomNumberList.pop();
	if(v.randomNumberList.empty())
		GenerateRandomNumbers(ref(v), seed);
	v.rand3 = v.randomNumberList.front();
	v.randomNumberList.pop();

	
	v.minIndex = selectMinIndex(ref(v), v.rand2);
	
	v.accepted=false;
	if(v.rand3<=(v.p_lower_bound[v.minIndex]/v.p_upper_bound[v.minIndex]))
		v.accepted=true;
	else
	{
		computePropensity(ref(v), v.minIndex);
		if(v.rand3<=(v.p[v.minIndex]/v.p_upper_bound[v.minIndex]))
			v.accepted=true;
	}
	
	/************** Firing the chosen reaction ************/
	if(v.accepted)
	{
		//cout<<"trying to fire new reaction ";
		switch(withDelay[v.minIndex])
		{
			case 1: //CD 
						{
						//cout<<"CD R"<<minIndex<<" goes to queue"<<endl;
						for(int i=0;i<(int)reactant[v.minIndex].size(); i++)
						{
							v.idOfElement   = reactant[v.minIndex][i].first;
							v.changedAmount = reactant[v.minIndex][i].second;
							v.x[v.idOfElement] -= v.changedAmount;
						}							
						//cout<<"CD reaction goes to queue becaus of delay."<<endl;
						v.q.push(make_pair(v.t+tDelay, v.minIndex));
						break;
						}
			case 2: //NCD	
						{
						//cout<<"NCD R"<<minIndex<<" goes to queue"<<endl;								
						//cout<<"NCD reaction goes to queue becaus of delay."<<endl;
						v.q.push(make_pair(v.t+tDelay, v.minIndex));
						break;
						}
			case 3: //ND	
						{
						//cout<<"ND R"<<minIndex<<endl;
						for(int i=0;i<(int)reactant[v.minIndex].size(); i++)
						{
							v.idOfElement   = reactant[v.minIndex][i].first;
							v.changedAmount = reactant[v.minIndex][i].second;
							v.x[v.idOfElement] -= v.changedAmount;
						}				
						for(int i=0; i<(int)product[v.minIndex].size(); i++)
						{
							v.idOfElement   = product[v.minIndex][i].first;
							v.changedAmount = product[v.minIndex][i].second;
							v.x[v.idOfElement] += v.changedAmount;
						}							
						//cout<<"CD reaction goes to queue becaus of delay."<<endl;
						break;
						}
		}
	}

};

void Simulation_Worker()
{
	Variables sim[num_simulation];
	
	int n=num_simulation/num_worker;
	
	/************** Starting simulations************/	
	ff::ParallelFor pf(num_worker);
	pf.parallel_for(0L, num_simulation, 1, 0, [&](const long s)
	{
		for(int mysim=0; mysim<n; mysim++)
		{
			//cout<<endl<<"Starting sim: "<<s<<endl;
			sim[s].rand_seed = (double)rand();
			
			sim[s].t=0.0;
			sim[s].x[0] = 1000;
			sim[s].x[1] = 0;	
			sim[s].x[2] = 0;
			while(!sim[s].q.empty())
				sim[s].q.pop();
			
			for(int i=0; i<M; i++)
				computeSpeciesBounds(ref(sim[s]), i);
			
			for(int i=0; i<N; i++)
				computePropensityBounds(ref(sim[s]), i);
				
			sim[s].total_propensity_upper_bound = propensityUpperBoundSum(ref(sim[s]));
			
			while(sim[s].t<ENDTIME)
			{
				while(checkSpeciesBound(ref(sim[s]))==true && sim[s].t<ENDTIME)
				{
					//cout<<"Calling delayedRejectionBasedSelection by sim "<<s<<"  t= "<<sim[s].t<<endl;
					delayedRejectionBasedSelection(ref(sim[s]), sim[s].rand_seed, s);					
				}
			
				/**************Updating the propensities of the Species************/
				//cout<<"Updating the propensities"<<endl;
				for(int i=0; i<M; i++)
				{
					if(sim[s].x[i]<sim[s].x_lower_bound[i] || sim[s].x[i]>sim[s].x_upper_bound[i])
					{
						computeSpeciesBounds(ref(sim[s]), i);
						for(int j=0; j<N; j++)
						{
							for(int k=0; k<(int)dependencyGraph[j].size(); k++)
							{
	/*------------*/			double old_value = sim[s].p_upper_bound[k];
								computePropensityBounds(ref(sim[s]), k);
								sim[s].total_propensity_upper_bound +=  sim[s].p_upper_bound[k] - old_value;
							}
						}
					}
				}
			}
			//cout<<endl<<"Ending sim: "<<s<<endl;
		}
	});
	
	
	for(int s=1; s<num_simulation; s++)
		for(int i=0; i<M; i++)
			sim[0].x[i] += sim[s].x[i];
	
	//the average result of the simulations 
	for(int i=0; i<M; i++)
	{
		sim[0].x[i] = sim[0].x[i]/num_simulation;
		cout<<sim[0].x[i]<<" ";
	}
	cout<<endl;
}



int main(int argc, char * argv[])
{
	 if (argc<4) 
	 {
        std::cerr << "use: " << argv[0]  << "num_simulation num_threads_rand num_worker\n";
        return -1;
    }
	
	
	num_simulation	= std::stol(argv[1]);	// The number of simulations	
	num_thread_rand	= std::stol(argv[2]);	// The number of workers
	num_worker 		= std::stol(argv[3]);	// The number of workers
	
	init();
	SR_graph();
	
	
//************************** Fast Flow ************************	
	auto startReac = chrono::high_resolution_clock::now();
	Simulation_Worker();
	auto elapsedReac = chrono::high_resolution_clock::now() - startReac;	
	cout << "Time	chrono: " << chrono::duration_cast<chrono::milliseconds>(elapsedReac).count() << " milliseconds"<< endl; 
	//std::cout << "Time	pipe.ffTime(): " << farm.ffTime() << "\n";	
//*****************************************************************
	
	
//pipe.ffStats(std::cout);
	return 0;
}

