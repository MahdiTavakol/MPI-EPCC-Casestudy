template <typename type>
void deallocate(type ** &v);
template <typename type>
type **allocate(type ** &v, const int nx, const int ny);
template <typename type>
type **allocate_cpp(type ** &v, const int nx, const int ny)
template <typename type>
void deallocate_cpp(type ** &v);

void pgmsize (char *filename, int *nx, int *ny);
void pgmread (char *filename, void *vx, int nx, int ny);
void pgmwrite(char *filename, void *vx, int nx, int ny);
