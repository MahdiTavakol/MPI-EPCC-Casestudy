#include <string>
#include <cstring>
#include <iostream>
#include <ctime>
#include "pgmio.h"
#define MAXFILENAME 128

using namespace std;
int main(int argc, char** argv)
{
    int iterations[] = {0,1,10,100,1000,10000,50000,100000,200000};
    int M;
    int N;
	double *buff, *new_ptr, *old, *edge;
    string filename;
    char* filename_c;
    char* outputfilename_c;
    char Mchar[3];
    char Nchar[3];

    M = stoi(argv[1]);
    N = stoi(argv[2]);

    filename = "edge" + to_string(M) + "x" + to_string(N) + ".pgm";
    filename_c = new char[filename.length()+1];
    strcpy(filename_c,filename.c_str());
    pgmsize(filename_c, &M, &N);

    buff     = new double[M*N];
    new_ptr  = new double[(M+2)*(N+2)];
    old      = new double[(M+2)*(N+2)];
    edge     = new double[(M+2)*(N+2)];

    clock_t start = clock();
    for (int k = 0; k < sizeof(iterations)/sizeof(iterations[0]); k++)
    {
        int iteration = iterations[k];
        for (int i = 0; i < (M+2)*(N+2) ; i++)
        {
            edge[i] = 255;
            old[i]  = 255;
        }
        pgmread(filename_c,buff,M,N);
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
                {
                    int indx0 = (N+2)-(j+1)-1+(N+2)*(i+1);
                    int indx1 = N - j  -1 + N*i;
                    edge[indx0] = buff[indx1];
                }
        }
        for (int l = 0; l < iteration; l++)
        {
            //cout << k << " step " << l << " out of " << iteration << endl;
            for (int i = 1; i < M+1; i++)
                for (int j = 1; j < N+1; j++)
                {
                    //cout << k << " step " << l << " out of " << iteration << endl;
                    //cout << "Working on " << i << "," << j << endl;
                    int indx0 = (N+2)-  j  -1+(N+2)*i;
                    int indx1 = (N+2)-  j  -1+(N+2)*(i-1);
                    int indx2 = (N+2)-  j  -1+(N+2)*(i+1);
                    int indx3 = (N+2)-(j-1)-1+(N+2)*i;
                    int indx4 = (N+2)-(j+1)-1+(N+2)*i;
                    int indx5 = (N+2)-  j  -1+(N+2)*i;
                    //cout << "Indexes are " << indx0 << "," << indx1 << "," << indx2 << "," << indx3 << "," << indx4 << "," << indx5 << endl;
                    new_ptr[indx0] = 0.25*(old[indx1]+old[indx2]+old[indx3]+old[indx4]-edge[indx5]);
                }
            for (int i = 1; i < M+1; i++)
                for (int j = 1; j < N+1; j++)
                {
                    int indx = (N+2)-j-1+(N+2)*i;
                    old[indx] = new_ptr[indx];
                }
        }

        for (int i = 1; i < M+1; i++)
            for (int j = 1; j < N+1; j++)
            {
                int indx0 = N-(j-1)-1+N*(i-1);
                int indx1 = (N+2) - j -1 + (N+2)*i;
                buff[indx0] = old[indx1];  
            }

        string outputfilename = "edge" + to_string(M) + "x" + to_string(N) + "-iteration-" + to_string(iteration) + ".pgm";
        outputfilename_c = new char[outputfilename.length()+1];
        strcpy(outputfilename_c, outputfilename.c_str());
        pgmwrite(outputfilename_c,buff,M,N);
    }
    clock_t end = clock();
    float numIters = 0;

    for (int i = 0; i < sizeof(iterations)/sizeof(iterations[0]); i++)
        numIters = numIters + (float)iterations[i];
    float timing = (float(end - start)/ CLOCKS_PER_SEC)*1000.0 / numIters;
	printf("It takes %f ms on 1 CPU\n",timing);

    delete [] buff;
    delete [] new_ptr;
    delete [] old;
    delete [] edge;

    return 0;
}
