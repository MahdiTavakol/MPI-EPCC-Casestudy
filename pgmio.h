template (typename type)
void deallocate(type ** &v);
template (typename type)
type **allocate(type ** &v, const int nx, const int ny);
void pgmsize (char *filename, int *nx, int *ny);
void pgmread (char *filename, void *vx, int nx, int ny);
void pgmwrite(char *filename, void *vx, int nx, int ny);
