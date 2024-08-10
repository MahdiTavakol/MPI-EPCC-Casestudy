#include <string>
#include <iostream>
#include "mpi.h"
#include "pgmio.h"

using namespace std;
int main(int argc, char** argv)
{
    int iterations[] = {0,1,10,100,1000,10000,50000,100000,200000};
    int M;
    int N;
    int Mp;
    int Np;
    double **masterbuff, **buff, **new_ptr, **old, **edge;
    std::string filename;
    std::string outputfilename;

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
    
    // MPI topology
    MPI_Comm topology;
    int dim = size;
    int period = 0;
    int reorder = 0;
    int left_proc;
    int right_proc;

    MPI_Cart_create(world,1,&dim,&period,reorder,&topology);
    MPI_Cart_shift(topology,0,1,&left_proc,&right_proc);



    filename = "edge" + to_string(M) + "x" + to_string(N) + ".pgm";
    pgmsize(filename, &M, &N);

    Mp = M/size;
    Np = N;

    masterbuff = allocate(masterbuff,M,N);
    buff       = allocate(buff,Mp,Np);
    new_ptr    = allocate(new_ptr,Mp+2,Np+2);
    old        = allocate(old,Mp+2,Np+2);
    edge       = allocate(edge,Mp+2,Np+2);
    


    MPI_Barrier(world);
    float start = MPI_Wtime();
    for (int k = 0; k < sizeof(iterations)/sizeof(iterations[0]); k++)
    {
        int iteration = iterations[k];
        if (rank == 0) 
            cout << "Starting the iteration " << iteration  << endl;
        for (int i = 0; i < Mp+2 ; i++)
	    for (int j = 0; j < Np+2 ; j++)
	    {
		edge[i][j] = 255;
		old[i][j] = 255;
	    }
	   

        if (rank == 0)
        {
            pgmread(filename,masterbuff,M,N);
            cout << "Finished reading " << filename  << endl;
        }
        MPI_Scatter(masterbuff[0],Mp*Np,MPI_DOUBLE,buff[0],Mp*Np,MPI_DOUBLE,0,world);
        for (int i = 0; i < Mp; i++)
            for (int j = 0; j < Np; j++)
		edge[i+1][j+1] = buff[i][j];
	    
        for (int l = 0; l < iteration; l++)
        {
            //cout << k << " step " << l << " out of " << iteration << endl;
            for (int i = 1; i < Mp+1; i++)
                for (int j = 1; j < Np+1; j++)
                {
                    //cout << k << " step " << l << " out of " << iteration << endl;
                    //cout << "Working on " << i << "," << j << endl;
		    new_ptr[i][j] = 0.25*(old[i+1][j]+old[i-1][j]+old[i][j+1]+old[i][j-1]-edge[i][j]);
                }
            for (int i = 1; i < Mp+1; i++)
                for (int j = 1; j < Np+1; j++)
                    old[i][j] = new_ptr[i][j];
                    
            //Halo swapping
            int tag1  = 0;
            int indx0 = 1;
            int indx1 = (Mp+1);
            MPI_Issend(old[indx0],Np,MPI_DOUBLE,left_proc,tag1,world,&request);
            MPI_Recv(old[indx1],Np,MPI_DOUBLE,right_proc,tag1,world,MPI_STATUS_IGNORE);
            int tag2  = 1;
            int indx2 = Mp;
            int indx3 = 0;
            MPI_Issend(old[indx2],Np,MPI_DOUBLE,right_proc,tag2,world,&request);
            MPI_Recv(old[indx3],Np,MPI_DOUBLE,left_proc,tag2,world,MPI_STATUS_IGNORE);
        }

        for (int i = 1; i < Mp+1; i++)
            for (int j = 1; j < Np+1; j++)
                buff[i-1][j-1] = old[i][j];  

        MPI_Gather(buff[0], Mp*Np, MPI_DOUBLE,masterbuff[0],Mp*Np,MPI_DOUBLE,0,world);
        if (rank == 0)
        {
            string outputfilename = "edge" + to_string(M) + "x" + to_string(N) + "-iteration-" + to_string(iteration) + ".pgm";
            pgmwrite(outputfilename,masterbuff,M,N);
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

    deallocate(masterbuff);
    deallocate(buff);
    deallocate(new_ptr);
    deallocate(old);
    deallocate(edge);
    MPI_Finalize();

    return 0;
}
