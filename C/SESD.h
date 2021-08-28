void igrid(int *grid, int size);
void fill(int *grid, int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe, int nthreads);
void ses(int *grid, int nx, int ny, int nz, double step, double probe, int nthreads);
int define_surface_points(int *grid, int nx, int ny, int nz, int i, int j, int k);
void filter_surface(int *grid, int nx, int ny, int nz, int nthreads);
void detect(int *grid, int size, int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe, int is_ses, int nthreads, int verbose);