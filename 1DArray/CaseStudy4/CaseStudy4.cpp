#include <string>
#include <cstring>
#include <iostream>
#include "mpi.h"
#include "pgmio.h"
#define MAXFILENAME 128

using namespace std;
int main(int argc, char** argv)
{
    int iterations[] = {0,1,10,100,1000,10000,50000,100000,200000};
    int M;
    int N;
    int Mp;
    int Np;
	double *masterbuff, *buff, *new_ptr, *old, *edge;
    string filename;
    char* filename_c;
    char* outputfilename_c;

    M = stoi(argv[1]);
    N = stoi(argv[2]);

    int rank;
    int size;

    MPI_Comm world = MPI_COMM_WORLD;
    MPI_Request request;
    MPI_Status status;

    MPI_Init(NULL,NULL);
    MPI_Comm_rank(world,&rank);
    MPI_Comm_size(world,&size);


    filename = "edge" + to_string(M) + "x" + to_string(N) + ".pgm";
    filename_c = new char[filename.length()+1];
    strcpy(filename_c,filename.c_str());
    pgmsize(filename_c, &M, &N);

    Mp = M/size;
    Np = N;

    masterbuff = new double[M*N];
    buff       = new double[Mp*Np];
    new_ptr    = new double[(Mp+2)*(Np+2)];
    old        = new double[(Mp+2)*(Np+2)];
    edge       = new double[(Mp+2)*(Np+2)];


    MPI_Barrier(world);
    float start = MPI_Wtime();
    for (int k = 0; k < sizeof(iterations)/sizeof(iterations[0]); k++)
    {
        int iteration = iterations[k];
        if (rank == 0) 
            cout << "Starting the iteration " << iteration  << endl;
        for (int i = 0; i < (Mp+2)*(Np+2) ; i++)
        {
            edge[i] = 255;
            old[i]  = 255;
        }
        if (rank == 0)
        {
            pgmread(filename_c,masterbuff,M,N);
            cout << "Finished reading " << filename  << endl;
        }
        MPI_Scatter(masterbuff,Mp*Np,MPI_DOUBLE,buff,Mp*Np,MPI_DOUBLE,0,world);
        for (int i = 0; i < Mp; i++)
        {
            for (int j = 0; j < Np; j++)
                {
                    int indx0 = (Np+2)-(j+1)-1+(Np+2)*(i+1);
                    int indx1 = Np - j  -1 + Np*i;
                    edge[indx0] = buff[indx1];
                }
        }
        for (int l = 0; l < iteration; l++)
        {
            //cout << k << " step " << l << " out of " << iteration << endl;
            for (int i = 1; i < Mp+1; i++)
                for (int j = 1; j < Np+1; j++)
                {
                    //cout << k << " step " << l << " out of " << iteration << endl;
                    //cout << "Working on " << i << "," << j << endl;
                    int indx0 = (Np+2)-  j  -1+(Np+2)*i;
                    int indx1 = (Np+2)-  j  -1+(Np+2)*(i-1);
                    int indx2 = (Np+2)-  j  -1+(Np+2)*(i+1);
                    int indx3 = (Np+2)-(j-1)-1+(Np+2)*i;
                    int indx4 = (Np+2)-(j+1)-1+(Np+2)*i;
                    int indx5 = (Np+2)-  j  -1+(Np+2)*i;
                    //cout << "Indexes are " << indx0 << "," << indx1 << "," << indx2 << "," << indx3 << "," << indx4 << "," << indx5 << endl;
                    new_ptr[indx0] = 0.25*(old[indx1]+old[indx2]+old[indx3]+old[indx4]-edge[indx5]);
                }
            for (int i = 1; i < Mp+1; i++)
                for (int j = 1; j < Np+1; j++)
                {
                    int indx = (Np+2)-j-1+(Np+2)*i;
                    old[indx] = new_ptr[indx];
                }
        }

        for (int i = 1; i < Mp+1; i++)
            for (int j = 1; j < Np+1; j++)
            {
                int indx0 = Np-(j-1)-1+Np*(i-1);
                int indx1 = (Np+2) - j -1 + (Np+2)*i;
                buff[indx0] = old[indx1];  
            }

        MPI_Gather(buff, Mp*Np, MPI_DOUBLE,masterbuff,Mp*Np,MPI_DOUBLE,0,world);
        if (rank == 0)
        {
            string outputfilename = "edge" + to_string(M) + "x" + to_string(N) + "-iteration-" + to_string(iteration) + ".pgm";
            outputfilename_c = new char[outputfilename.length()+1];
            strcpy(outputfilename_c, outputfilename.c_str());
            pgmwrite(outputfilename_c,masterbuff,M,N);
        }
    }

    MPI_Barrier(world);
    float end = MPI_Wtime();
    float numIters = 0;

    for (int i = 0; i < sizeof(iterations)/sizeof(iterations[0]); i++)
        numIters = numIters + (float)iterations[i];
    float timing = (end - start)*1000.0 / numIters;
    if (rank == 0)
		    printf("It takes %f ms on %d CPUs\n",timing,size);
    
    delete [] buff;
    delete [] new_ptr;
    delete [] old;
    delete [] edge;
    MPI_Finalize();

    return 0;
}
