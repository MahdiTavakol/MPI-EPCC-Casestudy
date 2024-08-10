/*
 * This file contains Cpp routines for the MPI Casestudy.
 *
 * "pgmread" reads in a PGM picture and can be called as follows:
 *
 *    double buf[M][N];
 *    pgmread("edge.pgm", buf, M, N);
 *
 * "pgmwrite" writes an array as a PGM picture and can be called as follows:
 *
 *    double buf[M][N];
 *    pgmwrite("picture.pgm", buf, M, N);
 *
 * "pgmsize" returns the size of a PGM picture and can be called as follows:
 *
 *    int nx, ny;
 *    pgmsize("edge.pgm", &nx, &ny);
 *
 *  To access these routines, add the following to your program:
 *
 *    #include "pgmio.h"
 *
 *  Note: you MUST link with the maths library -lm to access fabs etc.
 */


#include <iostream>
#include <fstream>
#include <string>
#include <stringstream>
#include <cmath>

/*  Routines to allocate/deallocate a contiguous storage on a 2D array
*/

template <typename type>
type **allocate(type ** &v, const int nx, const int ny)
{
  bigint nbytes = (bigint) sizeof(type) * nx * ny;
  type * data = (type *) malloc(nbytes);
  nbytes = (bigint) sizeof(type *) * nx;
  v = (type **) malloc(nbytes);

  int n = 0;
  for (int i = 0; i < nx; i++)
    {
      v[i] = &data[n];
      n += ny;
    }
  return v;
}

template <typename type>
type **allocate_cpp(type ** &v, const int nx, const int ny)
{
  int n = nx * ny;
  type * data = new type[n];
  v = new type* [nx];
  
  int n = 0;
  for (int i = 0; i < nx; i++)
    {
      v[i] = &data[n];
      n += ny;
    }
  return v;
}

template <typename type>
void deallocate_cpp(type ** &v)
{
  if (v== nullptr) return;
  delete [] v[0];
  delete v;
  v = nullptr;
}

template <typename type>
type **allocate_column_major(type ** &v, const int nx, const int ny)
{
  bigint nbytes = (bigint) sizeof(type) * nx * ny;
  type * data = (type *) malloc(nbytes);
  nbytes = (bigint) sizeof(type *) * ny;
  v = (type **) malloc(nbytes);

  int n = 0;
  for (int j = 0; j < ny; i++)
    {
      v[j] = &data[n];
      n += nx;
    }
  return v;
}

template <typename type>
void deallocate(type ** &v)
{
  if (v == nullptr) return;
  free(v[0]);
  free(v);
  v = nullptr;
}

/*
 *  Routine to get the size of a PGM data file
 *
 *  Note that this assumes a single line comment and no other white space.
 */

void pgmsize(std::string filename, int *nx, int *ny)
{ 
  std::fstream file;

  std::string line;
  std::stringstream iss;

  if (file.open(filename,"r") == NULL)
  {
    std::cerr << "pgmsize: cannot open << filename << std::endl;
    return ;
  }

  std::getline(file,line);
  std::getline(file,line);
  std::getline(file,line);

  iss << line;
  iss >> nx >> ny;

  file.close();
}

/*
 *  Routine to read a PGM data file into a 2D floating point array
 *  x[nx][ny]. 
 */


void pgmread(std::string filename, double ** v, int nx, int ny)
{ 
  std::fstream file;

  std::string line;
  std::stringstream iss;

  int nxt, nyt;

  if (file.open(filename,"r") == NULL)
  {
    std::cerr << "pgmsize: cannot open << filename << std::endl;
    return ;
  }

  std::getline(file,line);
  std::getline(file,line);
  std::getline(file,line);
  iss << line;
  iss >> nxt >> nyt;

  if (nx != nxt || ny != nyt)
  {
    std::cerr << "pgmread: size mismatch, (nx,ny) = ("<< nxt << "," << nyt << ") expected (" << nx << "," << ny << ")" << std::endl;
    return ;
  }

  std::getline(file,line); // Maximum greyscale

  for (int j = 0; j < ny; j++)
    {
      std::getline(file,line);
      iss.clear();
      iss.str("");
      iss << line;
      for (int i = 0; i < nx ; i++)
        {
          int n;
          iss >> n;
          v[i][j] = n;
        }
  file.close();
}


/*
 *  Routine to write a PGM image file from a 2D floating point array
 *  x[nx][ny]. 
 */

void pgmwrite(std::string filename, double **x, int nx, int ny)
{
  std::fstream file;
  std::string line;
  

  int i, j, k, grey;

  double xmin, xmax, tmp, fval;
  double thresh = 255.0;


  if (!file.open(filename,"w"))
  {
    std::cerr << "pgmwrite: cannot create" << filename << std::endl;
    return ;
  }

  cout << "Writing " << nx << " x " << ny << " picture into file:" << filename << std::endl;


  /*
   *  Find the max and min absolute values of the array
   */

  xmin = fabs(x[0]);
  xmax = fabs(x[0]);

  for (i=0; i < nx*ny; i++)
  {
    if (fabs(x[i]) < xmin) xmin = fabs(x[i]);
    if (fabs(x[i]) > xmax) xmax = fabs(x[i]);
  }

  if (xmin == xmax) xmin = xmax-1.0;

  file <<  "P2" << std::endl;
  file << "# Written by pgmio::pgmwrite" << std::endl;
  file << nx << " "  << ny << std::endl;
  file <<  thresh << std::endl;

  k = 0;

  for (j=ny-1; j >=0 ; j--)
  {
    for (i=0; i < nx; i++)
    {
      /*
       *  Access the value of x[i][j]
       */

      tmp = x[i][j];

      /*
       *  Scale the value appropriately so it lies between 0 and thresh
       */

      fval = thresh*((fabs(tmp)-xmin)/(xmax-xmin))+0.5;
      grey = (int) fval;

      cout << std::left << std::setw(6) << grey;

      if (0 == (k+1)%18) cout << std::endl;

      k++;
    }
  }

  if (0 != k%18) cout << std::endl;
  file.close();
}
