// qaProject.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>

#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <string>
#include <vector>
#include <tuple>

using namespace std;

//declare consts
const int parameterStart = 50;
const int parameterEnd = 300;
const int incrementor = 50;
const double timeDiffThreshold = .005;

// functions
double sequentialRun(int N, int M);
double parrallelOptimizedRun(int N, int M);

int main()
{
    int threads = 2, maxThreads = 5;
    //try for threads 2-5???
    for (threads = 2; threads <= maxThreads; threads++)
    {
        // the number of threads may stay deterministic
        omp_set_num_threads(threads);
        

        // you may add some code here to measure the execution time
        double seqTime, paraTime;
        double seqWins = 0, totalRuns = 0;
        double loseRate = 0;
        vector<tuple<int, int, double, double>> sequentialFasterCases;

        //sequential
        for (int N = parameterStart; N <= parameterEnd; N += incrementor)
        {
            for (int M = parameterStart; M <= parameterEnd; M += incrementor)
            {
                //sequential testing
                double seqStart, seqEnd;

                //get start time
                seqStart = omp_get_wtime();

                //run sequential function
                double result = sequentialRun(N, M);
                //get stop time
                seqEnd = omp_get_wtime();
                //get total time
                seqTime = seqEnd - seqStart;

                //para testing
                double paraStart, paraEnd;
                paraStart = omp_get_wtime();
                double resultParallel = parrallelOptimizedRun(N, M);
                paraEnd = omp_get_wtime();
                paraTime = paraEnd - paraStart;
                
                if (seqTime < paraTime && (paraTime - seqTime) >= timeDiffThreshold)
                {
                    sequentialFasterCases.emplace_back(N, M, seqTime, paraTime);
                    seqWins++;
                }

                totalRuns++;
            }
        }
        loseRate = seqWins / totalRuns;
        printf("\n\nPercentage of Sequential Wins for %d threads is %1.2f\n\n", threads, loseRate);
        
        // handle sequential wins
        for (const auto& caseData : sequentialFasterCases)
        {
            int N_case, M_case;
            double Seq_case, Para_case;
            tie(N_case, M_case, Seq_case, Para_case) = caseData;
            printf("N+M = %d\n", N_case + M_case);
            printf("Sequential is faster for N=%d , M=%d, by %f\n\n", N_case, M_case, Para_case - Seq_case);
        }
    }

    return 0;
}

double sequentialRun(int N, int M)
{
    long i, j;
    long A = 0;
    double B = 0, C = 0, D = 0;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            A += i * j;
        }
    }

    for (i = 1; i < (long)sqrt(A); i++)
        B += 1 / i;

    for (i = 0; i < M * N; i++)
        for (j = 0; j < M; j++)
        {
            D += pow(0.1, i * j);
        }

    for (i = 0; i < (long)B * (N + 1); i++)
        for (j = 1; j < (long)sqrt(D); j++)
        {
            C += i / j;
        }

    return A + B - C / D;
}

double parrallelOptimizedRun(int N, int M)
{
    long i, j;
    long A = 0;
    double B = 0, C = 0, D = 0;

    // handle nested loop with data dependency A
        //fast paralellization
    //if(N * M > 100)
#pragma omp parallel for if(N + M < 350) private(j) reduction(+:A) schedule(dynamic)
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            A += i * j;
        }
    }

    // data dependency
    //if(N > 100)
    #pragma omp parallel for if(N + M < 350) reduction(+:B) schedule(dynamic)
    for (i = 1; i < (long)sqrt(A); i++)
        B += 1 / i;

    //nested loops data dependency D and B
    //if(N * M > 100)
    #pragma omp parallel for if(N + M < 350) private(j) reduction(+:D) schedule(dynamic)
    for (i = 0; i < M * N; i++)
        for (j = 0; j < M; j++)
        {
            D += pow(0.1, i * j);
        }

    //nested loops
    //if((long)B * (N + 1) > 100)
    #pragma omp parallel for if(N + M < 350) reduction(+:C) schedule(dynamic)
    for (i = 0; i < (long)B * (N + 1); i++)
        for (j = 1; j < (long)sqrt(D); j++)
        {
            C += i / j;
        }

    //return result
    return A + B - C / D;
}