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
    int * count = new int[3];
    int * blklength = new int[3];
    int * stride = new int[3];
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
    MPI_Comm         world = MPI_COMM_WORLD;
    MPI_Request * requests = new MPI_Request[4];
    MPI_Status    * status = new MPI_Status[4];
    MPI_Datatype * mpitype = new MPI_Datatype[4];

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

    Mp = ceil((double)M/(double)size);
    Np = N;


    masterbuff = new double[Mp*Np*size];
    buff       = new double[Mp*Np];
    new_ptr    = new double[(Mp+2)*(Np+2)];
    old        = new double[(Mp+2)*(Np+2)];
    edge       = new double[(Mp+2)*(Np+2)];


    int iteration = maxiteration;
    if (rank == 0) 
        cout << "Starting the iteration " << iteration  << endl;


    // MPI derived types
    count[0] = 1;
    blklength[0] = N;
    stride[0] = N;

    count[1] = Mp;
    blklength[1] = Np;
    stride[1] = Np+2;

    MPI_Type_vector(count[0],blklength[0],stride[0],MPI_DOUBLE,&mpitype[0]);
    MPI_Type_commit(&mpitype[0]);

    MPI_Type_vector(count[1],blklength[1],stride[1],MPI_DOUBLE,&mpitype[1]);
    MPI_Type_commit(&mpitype[1]);

    // MPI derived type for the halo
    count[2] = 1;
    blklength[2] = Np + 2;
    stride[2] = Np + 2;



    MPI_Type_vector(count[2],blklength[2],stride[2],MPI_DOUBLE,&mpitype[2]);
    MPI_Type_commit(&mpitype[2]);



    for (int i = 0; i < (Mp+2)*(Np+2) ; i++)
    {
        edge[i] = 255;
        old[i] = 255;
    }


    MPI_Barrier(world);
    float start = MPI_Wtime();
    if (rank == 0)
    {
        pgmread(filename_c,masterbuff,M,N);
        cout << "Finished reading " << filename  << endl;
    }
    
    
    



    MPI_Scatter(masterbuff,Mp,mpitype[0],&edge[Np+3],1,mpitype[1],0,world);


    for (int l = 0; l < iteration; l++)
    {
        //Halo swapping
        int tag1  = 0;
        int indx0 = (Np+2)*1;
        int indx1 = (Np+2)*(Mp+1);
        MPI_Issend(&old[indx0],1,mpitype[2],left_proc,tag1,world,&requests[0]); 
        MPI_Irecv(&old[indx1],1,mpitype[2],right_proc,tag1,world,&requests[1]);
        int tag2 = 1;
        int indx2 = (Np+2)*Mp;
        int indx3 = (Np+2)*(0);
        MPI_Issend(&old[indx2],1,mpitype[2],right_proc,tag2,world,&requests[2]);
        MPI_Irecv(&old[indx3],1,mpitype[2],left_proc,tag2,world,&requests[3]); 


        for (int i = 2; i < Mp; i++)
            for (int j = 2; j < Np; j++)
            {
                int indx0 = (Np+2)-  j  -1+(Np+2)*i;
                int indx1 = (Np+2)-  j  -1+(Np+2)*(i-1);
                int indx2 = (Np+2)-  j  -1+(Np+2)*(i+1);
                int indx3 = (Np+2)-(j-1)-1+(Np+2)*i;
                int indx4 = (Np+2)-(j+1)-1+(Np+2)*i;
                int indx5 = (Np+2)-  j  -1+(Np+2)*i;
                new_ptr[indx0] = 0.25*(old[indx1]+old[indx2]+old[indx3]+old[indx4]-edge[indx5]);
            }


        MPI_Wait(&requests[0],&status[0]);
        MPI_Wait(&requests[1],&status[1]);
        MPI_Wait(&requests[2],&status[2]);
        MPI_Wait(&requests[3],&status[3]);



        int is[] = {1,Mp};
        int js[] = {1,Np};

        for (int i: is)   
            for (int j = 2; j < Np; j++)
            {
                int indx0 = (Np+2)-  j  -1+(Np+2)*i;
                int indx1 = (Np+2)-  j  -1+(Np+2)*(i-1);
                int indx2 = (Np+2)-  j  -1+(Np+2)*(i+1);
                int indx3 = (Np+2)-(j-1)-1+(Np+2)*i;
                int indx4 = (Np+2)-(j+1)-1+(Np+2)*i;
                int indx5 = (Np+2)-  j  -1+(Np+2)*i;
                new_ptr[indx0] = 0.25*(old[indx1]+old[indx2]+old[indx3]+old[indx4]-edge[indx5]);
            }
        for (int i = 1; i < Mp+1; i++)
            for (int j: js)
            {
                int indx0 = (Np+2)-  j  -1+(Np+2)*i;
                int indx1 = (Np+2)-  j  -1+(Np+2)*(i-1);
                int indx2 = (Np+2)-  j  -1+(Np+2)*(i+1);
                int indx3 = (Np+2)-(j-1)-1+(Np+2)*i;
                int indx4 = (Np+2)-(j+1)-1+(Np+2)*i;
                int indx5 = (Np+2)-  j  -1+(Np+2)*i;
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

    }

    MPI_Gather(&old[Np+3], 1, mpitype[1],masterbuff,Mp,mpitype[0],0,world); // It should work since it is opposite the MPI_Scatter....
    if (rank == 0)
    {
        string outputfilename = "edge" + to_string(M) + "x" + to_string(N) + "-iteration-" + to_string(lastloop) + ".pgm";
        outputfilename_c = new char[outputfilename.length()+1];
        strcpy(outputfilename_c, outputfilename.c_str());
        pgmwrite(outputfilename_c,masterbuff,Mp*size,N);
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
