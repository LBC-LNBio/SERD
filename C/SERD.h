/* Grid initialization */
void igrid(int *grid, int size);

/* Grid filling */
void fill(int *grid, int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe, int nthreads);

/* Biomolecular surface representation */
int check_protein_neighbours(int *grid, int nx, int ny, int nz, int i, int j, int k);
void ses(int *grid, int nx, int ny, int nz, double step, double probe, int nthreads);

/* Surface points detection */
int define_surface_points(int *grid, int nx, int ny, int nz, int i, int j, int k);
void filter_surface(int *grid, int nx, int ny, int nz, int nthreads);

/* Filter noise points */
int remove_noise_points(int *grid, int nx, int ny, int nz, int i, int j, int k);
void filter_noise_points(int *grid, int nx, int ny, int nz, int nthreads);

/* Enclosed points removal - flood and fill algorithm */
int check_unclustered_neighbours(int *grid, int nx, int ny, int nz, int i, int j, int k);
void flood_and_fill(int *grid, int nx, int ny, int nz, int i, int j, int k, int tag);
void filter_enclosed_regions(int *grid, int nx, int ny, int nz, double step, int nthreads);

/* Solvent-exposed surface detection */
void _surface(int *grid, int size, int nx, int ny, int nz, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe, int is_ses, int nthreads, int verbose);

/* Solvent-exposed residues detection */
typedef struct node
{
    int pos;
    struct node *next;
} res;
res *create(int pos);
void insert(res **head, res *res_new);
char **_interface(int *grid, int nx, int ny, int nz, char **pdb, double *atoms, int natoms, int xyzr, double *reference, int ndims, double *sincos, int nvalues, double step, double probe, int nthreads, int verbose);
