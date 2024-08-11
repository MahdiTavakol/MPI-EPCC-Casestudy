#include <string>

/*  Routines to allocate/deallocate a contiguous storage on a 2D array
*/

template <typename type>
type **allocate(type ** &v, const int nx, const int ny)
{
  long int nbytes = (long int) sizeof(type) * nx * ny;
  type * data = (type *) malloc(nbytes);
  nbytes = (long int) sizeof(type *) * nx;
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
type **allocate_cpp(type ** &v, int nx, int ny)
{
  int n = nx * ny;
  type * data = new type[n];
  v = new type* [nx];
  
  n = 0;
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
type **allocate_row_major(type ** &v, int nx, int ny)
{
  long int nbytes = (long int) sizeof(type) * nx * ny;
  type * data = (type *) malloc(nbytes);
  nbytes = (long int) sizeof(type *) * ny;
  v = (type **) malloc(nbytes);

  int n = 0;
  for (int j = 0; j < ny; j++)
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


void pgmsize (std::string filename, int *nx, int *ny);
void pgmread_offset(std::string filename, double ** v, int nx, int ny);
void pgmread (std::string filename, double **vx, int nx, int ny);
void pgmwrite(std::string filename, double **vx, int nx, int ny);
