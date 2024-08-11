#include <string>
#include <cstring>
#include <iostream>
#include <cmath>
#include "mpi.h"
#include "pgmio.h"
#define MAXFILENAME 128

using namespace std;
int main(int argc, char** argv)
{
    int maxiteration = 2000000;
    int lastloop;
    double maxdelta = 0.0001; //0.1 0.01 0.001 0.0001 0.00001
    int M;
    int N;
    int Mp;
    int Np;
	double *masterbuff, *buff, *new_ptr, *old, *edge;
    double delta, masterdelta;
    string filename;
    char* filename_c;
    char* outputfilename_c;

    M = stoi(argv[1]);
    N = stoi(argv[2]);

    int rank;
    int size;

    // MPI parameters
    MPI_Comm world = MPI_COMM_WORLD;
    MPI_Request request;
    MPI_Status status;

    // MPI initialization
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


    int iteration = maxiteration;
    if (rank == 0) 
        cout << "Starting the iteration " << iteration  << endl;
    for (int i = 0; i < (Mp+2)*(Np+2) ; i++)
    {
        edge[i] = 255;
        old[i]  = 255;
    }


    MPI_Barrier(world);
    float start = MPI_Wtime();
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
        delta = 0.0;
        for (int i = 1; i < Mp+1; i++)
            for (int j = 1; j < Np+1; j++)
            {
                int indx = (Np+2)-j-1+(Np+2)*i;
                double currentdelta = abs(old[indx]-new_ptr[indx]);
                if (currentdelta > delta)
                    delta = currentdelta;
                old[indx] = new_ptr[indx];
            }
        // Maximum change in the value
        if (l > 0)
        {
            MPI_Allreduce(&delta,&masterdelta,1,MPI_DOUBLE,MPI_MAX,world);
            if (masterdelta < maxdelta)
            {
                lastloop = l;
                break;
            }
        }

        //Halo swapping
        int tag1  = 0;
        int indx0 = (Np+2)*1;
        int indx1 = (Np+2)*(Mp+1);
        MPI_Issend(&old[indx0],Np,MPI_DOUBLE,left_proc,tag1,world,&request);
        MPI_Recv(&old[indx1],Np,MPI_DOUBLE,right_proc,tag1,world,MPI_STATUS_IGNORE);
        int tag2  = 1;
        int indx2 = (Np+2)*Mp;
        int indx3 = 0;
        MPI_Issend(&old[indx2],Np,MPI_DOUBLE,right_proc,tag2,world,&request);
        MPI_Recv(&old[indx3],Np,MPI_DOUBLE,left_proc,tag2,world,MPI_STATUS_IGNORE);
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
        string outputfilename = "edge" + to_string(M) + "x" + to_string(N) + "-iteration-" + to_string(lastloop) + ".pgm";
        outputfilename_c = new char[outputfilename.length()+1];
        strcpy(outputfilename_c, outputfilename.c_str());
        pgmwrite(outputfilename_c,masterbuff,M,N);
    }

    MPI_Barrier(world);
    float end = MPI_Wtime();

    float timing = (end - start)*1000.0 / lastloop;
    if (rank == 0)
		    printf("It takes %f ms on %d CPUs\n",timing,size);

    delete [] buff;
    delete [] new_ptr;
    delete [] old;
    delete [] edge;
    MPI_Finalize();

    return 0;
}
