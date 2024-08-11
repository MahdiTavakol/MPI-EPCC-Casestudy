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
#include <sstream>
#include <cmath>
#include <iomanip>


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
  
  file.open(filename,std::fstream::in);

  if (!file.is_open())
  {
    std::cerr << "pgmsize: cannot open" << filename << std::endl;
    return ;
  }

  std::getline(file,line);
  std::getline(file,line);
  std::getline(file,line);

  iss << line;
  iss >> *nx >> *ny;

  file.close();
}

/*
 *  Routine to read a PGM data file into a 2D floating point array
 *  x[nx][ny]. 
 */
 
void pgmread_offset(std::string filename, double ** v, int nx, int ny)
{ 
  std::fstream file;

  std::string line;
  std::stringstream iss;

  int nxt, nyt;

  file.open(filename,std::fstream::in);

  if (!file.is_open())
  {
    std::cerr << "pgmsize: cannot open" << filename << std::endl;
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

  iss.clear();
  iss.str("");
  std::getline(file,line);
  iss << line;

  for (int j = ny-1; j >= 0; j--)
  {

    for (int i = 0; i < nx ; i++)
    {
      if (iss.eof())
      {
         std::getline(file,line);
         iss.clear();
         iss.str("");
         iss << line;
      }
      int n;
      iss >> n;
      double m = (double) n;
      v[i+1][j+1] = n;
    }
  }
  file.close();
}


void pgmread(std::string filename, double ** v, int nx, int ny)
{ 
  std::fstream file;

  std::string line;
  std::stringstream iss;

  int nxt, nyt;

  file.open(filename,std::fstream::in);

  if (!file.is_open())
  {
    std::cerr << "pgmsize: cannot open" << filename << std::endl;
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

  iss.clear();
  iss.str("");
  std::getline(file,line);
  iss << line;

  for (int j = ny-1; j >= 0; j--)
  {

    for (int i = 0; i < nx ; i++)
    {
      if (iss.eof())
      {
         std::getline(file,line);
         iss.clear();
         iss.str("");
         iss << line;
      }
      int n;
      iss >> n;
      double m = (double) n;
      v[i][j] = n;
    }
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


  file.open(filename,std::fstream::out);

  if (!file.is_open())
  {
    std::cerr << "pgmwrite: cannot create" << filename << std::endl;
    return ;
  }

  std::cout << "Writing " << nx << " x " << ny << " picture into file:" << filename << std::endl;


  /*
   *  Find the max and min absolute values of the array
   */

  xmin = fabs(x[0][0]);
  xmax = fabs(x[0][0]);

  for (i=0; i < nx; i++)
    for (j=0; j < ny; j++)
    {
       if (fabs(x[i][j]) < xmin) xmin = fabs(x[i][j]);
       if (fabs(x[i][j]) > xmax) xmax = fabs(x[i][j]);
    }

  if (xmin == xmax) xmin = xmax-1.0;

  file <<  "P2" << std::endl;
  file << "# Written by pgmio::pgmwrite" << std::endl;
  file << " " << nx << " "  << ny << std::endl;
  file << " " << thresh << std::endl;

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

      file << std::right << std::setw(6) << grey;

      if (0 == (k+1)%12) file << std::endl;

      k++;
    }
  }

  if (0 != k%12) file << std::endl;
  file.close();
}
