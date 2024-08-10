#include <string>

template <typename type>
void deallocate(type ** &v);
template <typename type>
type **allocate(type ** &v, const int nx, const int ny);
template <typename type>
type **allocate_cpp(type ** &v, const int nx, const int ny)
template <typename type>
void deallocate_cpp(type ** &v);

void pgmsize (std::string filename, int *nx, int *ny);
void pgmread (std::string filename, void *vx, int nx, int ny);
void pgmwrite(std::string filename, void *vx, int nx, int ny);
