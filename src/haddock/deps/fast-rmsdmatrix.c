#define  _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <limits.h>
#include <stdbool.h>

typedef struct traj{
    int frames;             /*!< number of frames in the trajectory */
    double **traj_coords;   /*!< 2D array of xyz coordinates */
    int n_at;               /*!< number of atoms in the atomistic structure */
    int pairs;              /*!< number of possible pairs of structures */
}traj;

typedef struct alignments{
    double *rmsd_mat;           /*!< condensed pairwise RMSD matrix */
    double **coms;              /*!< array of centers of mass */
    int *ref_structs;           /*!< array of reference structures */
    int *mod_structs;           /*!< array of model structures */
}alignments;

void failed(char message[]) {
    printf("\n FAILED *** %s ***\n \n", message);
    exit(1);
}

struct arguments;

void zero_vec_d(double *a, int dim) {

    int i;
    for (i = 0; i < dim; i++) {
        a[i] = 0;
    }
}

int *i1t(int n1) {
    int *p, *a;
    int i;
    if ((p = (int *) malloc((size_t) n1 * sizeof(int))) == NULL)
        failed("i1t: failed");
    for (i = 0, a = p; i < n1; i++)
        *a++ = 0;
    return p;
}

double *d1t(int n1) {
    double *p, *a;
    int i;
    if ((p = (double *) malloc((size_t) n1 * sizeof(double))) == NULL)
        failed("d1t: failed n1 ");
    for (i = 0, a = p; i < n1; i++)
        *a++ = 0;
    return p;
}

double **d2t(int n1, int n2) {
    double **p, *a;
    int i;
    if ((p = (double **) malloc((size_t) n1 * sizeof(double *))) == NULL)
        failed("d2t: failed n1");
    if ((p[0] = (double *) malloc((size_t) n1 * n2 * sizeof(double))) == NULL)
        failed("d2t: failed n2");
    for (i = 0; i < n1 - 1; i++)
        p[i + 1] = p[i] + n2;
    for (i = 0, a = p[0]; i < n1 * n2; i++)
        *a++ = 0;
    return p;
}

void free_d2t(double **p) {

    free(p[0]);
    free(p);
}

void free_d1t(double *p) {

    free(p);
}

double scal_d(double *a, double *b, int dim) {

    int i;
    double temp;

    temp = 0.0;
    for (i = 0; i < dim; i++) {
        temp += a[i] * b[i];
    }
    return (temp);
}

double norm_d(double *a, int dim) {

    return (sqrt(scal_d(a, a, dim)));
}

void normalize_d(double *a, int dim) {
    int i;
    double temp;

    temp = norm_d(a, dim);
    for (i = 0; i < dim; i++) {
        a[i] = a[i] / temp;
    }
}

void vecprod_d(double *a, double *b, double *c) {

    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

void myjacobi(double a[][3], int n, double *d, double v[][3], int *nrot) {

    int j, iq, ip, i;
    double tresh, theta, tau, t, sm, s, h, g, c;
    double b[3], z[3];

#define ROTATE(a, i, j, k, l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

    for (ip = 0; ip <= (n - 1); ip++) {
        for (iq = 0; iq <= (n - 1); iq++)
            v[ip][iq] = 0.0;
        v[ip][ip] = 1.0;
    }
    for (ip = 0; ip <= (n - 1); ip++) {
        b[ip] = d[ip] = a[ip][ip];
        z[ip] = 0.0;
    }
    *nrot = 0;
    for (i = 1; i <= 500; i++) {
        sm = 0.0;
        for (ip = 0; ip <= n - 2; ip++) {
            for (iq = ip + 1; iq <= (n - 1); iq++)
                sm += fabs(a[ip][iq]);
        }
        if (sm == 0.0) {
            return;
        }
        if (i < 4)
            tresh = 0.2 * sm / (n * n);
        else
            tresh = 0.0;
        for (ip = 0; ip <= n - 2; ip++) {
            for (iq = ip + 1; iq <= (n - 1); iq++) {
                g = 100.0 * fabs(a[ip][iq]);
                if (i > 4 && (fabs((fabs(d[ip]) + g) - fabs(d[ip])) < 1.0e-6)
                    && (fabs((fabs(d[iq]) + g) - fabs(d[iq])) < 1.0e-6))
                    a[ip][iq] = 0.0;
                else if (fabs(a[ip][iq]) > tresh) {
                    h = d[iq] - d[ip];
                    if (fabs((fabs(h) + g) - fabs(h)) < 1.0e-6)
                        t = (a[ip][iq]) / h;
                    else {
                        theta = 0.5 * h / (a[ip][iq]);
                        t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                        if (theta < 0.0)
                            t = -t;
                    }
                    c = 1.0 / sqrt(1 + t * t);
                    s = t * c;
                    tau = s / (1.0 + c);
                    h = t * a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq] = 0.0;
                    for (j = 0; j <= ip - 1; j++) {
                        ROTATE (a, j, ip, j, iq)
                    }
                    for (j = ip + 1; j <= iq - 1; j++) {
                        ROTATE (a, ip, j, j, iq)
                    }
                    for (j = iq + 1; j <= (n - 1); j++) {
                        ROTATE (a, ip, j, iq, j)
                    }
                    for (j = 0; j <= (n - 1); j++) {
                        ROTATE (v, j, ip, j, iq)
                    }
                    ++(*nrot);
                }
            }
        }
        for (ip = 0; ip <= (n - 1); ip++) {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }
    printf("Too many iterations in routine JACOBI %lf\n", sm);
    /*  exit (1); */
#undef ROTATE
}

int columns(char* string){
    /**
    * routine that returns the number of columns for each row inside the file chosen.
    *
    * Parameter
    * ----------
    *
    * `string`   : string token in account
    */
 
    int counter = 0;

    while(string){
        counter++;
        string = strtok(NULL, " \t\v\r\n");
    }
    return counter;
}

int n_rows (FILE *f){
    /**
    * routine that returns the number of rows in a file. It counts correctly this number even if the last row does not present \n
    *
    * Parameter
    * ----------
    *
    * `f` : FILE structure that represents the file opened.
    */
 
    size_t rows_i = 0, n = 0;
    char *buf = NULL;

    fseek(f, 0, SEEK_SET);

    while(getline(&buf, &n, f) != -1) rows_i++;

    free(buf);

    return rows_i;
}

int get_pair(int nmodels, int idx, int ij_arr[2]){
    if (nmodels < 0 || idx < 0){
        printf("get_pair cannot accept negative numbers");
        printf("Input is %d , %d", nmodels, idx);
        exit(EXIT_FAILURE);
    }
    // solve the second degree equation
    int b = 1 - (2 * nmodels);
    int i = (-b - sqrt(b * b - 8 * idx)) / 2;
    int j = idx + i * (b + i + 2) / 2 + 1;
    ij_arr[0] = i;
    ij_arr[1] = j;
    return 0;
} 

double optimal_alignment(double **x, double **y, int cgnum, double u[][3]) {
    /**
    * routine that computes the Kabsch alignment and the rmsd between two configurations
    * 
    *
    * Parameters
    * ----------
    * 
    * `x`, `y` : CG structures
    * 
    * `cgnum` : length of CG mapping
    * 
    * `u` : rotation matrix
    */
    void myjacobi(double a[][3], int n, double *d, double v[][3], int *nrot);
    int i, j, k, sign[3], order[3], nrot;
    double e, e0;
    double r[3][3], rt[3][3], temp; //, **x, **y;
    double a[3][3], eval[3], evec[3][3];
    double eigenvalues[3], eigenvectors[3][3], b[3][3];
    int speak = 1;

    e0 = 0.0;
    for (i = 0; i < cgnum; i++) {
        e0 += 0.5 * norm_d(x[i], 3) * norm_d(x[i], 3);
        e0 += 0.5 * norm_d(y[i], 3) * norm_d(y[i], 3);
    }

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            r[i][j] = 0.0;
            for (k = 0; k < cgnum; k++) {
                r[i][j] += y[k][i] * x[k][j];
            }
            rt[j][i] = r[i][j];
        }
    }

    if (isnan(e0) == 1) printf("Found a NaN in Kabsch alignment at checkpoint 1\n");

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            a[i][j] = 0;
            for (k = 0; k < 3; k++) {
                a[i][j] += rt[i][k] * r[k][j];
            }
        }
    }

    myjacobi(a, 3, eval, evec, &nrot);
    /* we add small quantities in order to remove potentially dangerous degeneracies */
    if (isnan(eval[0]) == 1) printf("Found a NaN in Kabsch alignment at checkpoint 2\n");
    if (isnan(eval[1]) == 1) printf("Found a NaN in Kabsch alignment at checkpoint 3\n");
    if (isnan(eval[2]) == 1) printf("Found a NaN in Kabsch alignment at checkpoint 4\n");

    eval[0] += 0.0000000000000001;
    eval[1] += 0.00000000000000001;
    eval[2] += 0.00000000000000002;

    if ((eval[0] < eval[1]) && (eval[0] < eval[2])) {
        order[0] = 1;
        order[1] = 2;
        order[2] = 0;
    }

    if ((eval[1] < eval[0]) && (eval[1] < eval[2])) {
        order[0] = 0;
        order[1] = 2;
        order[2] = 1;
    }

    if ((eval[2] < eval[0]) && (eval[2] < eval[1])) {
        order[0] = 0;
        order[1] = 1;
        order[2] = 2;
    }

    for (i = 0; i < 3; i++) {
        eigenvalues[i] = eval[order[i]];
        for (j = 0; j < 3; j++) {
            eigenvectors[i][j] = evec[j][order[i]];
        }
    }

    normalize_d(eigenvectors[0], 3);
    normalize_d(eigenvectors[1], 3);
    vecprod_d(eigenvectors[0], eigenvectors[1], eigenvectors[2]);
    normalize_d(eigenvectors[2], 3);

    for (i = 0; i < 2; i++) {
        for (j = 0; j < 3; j++) {
            b[i][j] = 0;
            for (k = 0; k < 3; k++) {
                b[i][j] += r[j][k] * eigenvectors[i][k];
            }
        }
        normalize_d(b[i], 3);
    }

    vecprod_d(b[0], b[1], b[2]);
    normalize_d(b[2], 3);

    temp = 0.0;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            temp += b[2][i] * r[i][j] * eigenvectors[2][j];
        }
    }
    sign[2] = +1;
    if (temp < 0)
        sign[2] = -1;
    sign[0] = sign[1] = 1;

    if (fabs(eigenvalues[2]) < 0.0000001) {
        e = e0 - sqrt(eigenvalues[0]) - sqrt(eigenvalues[1]);
    } else {
        e = e0 - sqrt(eigenvalues[0]) - sqrt(eigenvalues[1]) - sign[2] * sqrt(eigenvalues[2]);
    }

    if (isnan(e) == 1) {
        printf("Found a NaN in Kabsch alignment at checkpoint 5 | \n");
        printf("e %lf e0 %lf e1 %lf e2 %lf e3 %lf\n", e, e0, eigenvalues[0], eigenvalues[1], eigenvalues[2]);
    }
/********************/
    e = 2.0 * e / cgnum;
    if (e < 0.0) {
        if (fabs(e) < 1.0e-3) {
            if (speak == 1) {
                //printf("Warning. In Kabsch alignment found slightly negative value of e (%e). Identical structures? Setting it to 0.000001.\n",
                //       e);
                e = 0.000001;
            }
        }
            /* occasionally, when dealing with two practically identical configurations
                the value of e may be slightly negative due to the small offsets and roundoff errors.
                In this case we set it equal to zero. */
        else {
            FILE *fe; 
            fe = fopen("error.dat", "w");
            fprintf(fe, "Error. In Kabsch alignment found negative value of e: (%lf)\n", e);
            printf("Error. In Kabsch alignment found negative value of e: (%lf)\n", e);
            fclose(fe);
            exit(EXIT_FAILURE);
        }
    }
    e = sqrt(e);
    if (isnan(e) == 1) printf("Found a NaN in Kabsch alignment at checkpoint 6\n");
/********************/
    // filling rotation_matrix
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            u[i][j] = 0.0;
            for (k = 0; k < 3; k++) {
                u[i][j] += b[k][i] * eigenvectors[k][j];
            }
        }
    }
    return (e);
}

void cycle_ilrmsd(traj *Trajectory, traj *Ligand_Trajectory, alignments *align, int start_pair[2]) {
    /**
    * routine that cycles over pairs of frames in a trajectory, aligning on the receptor
    * and computing the interface ligand rmsd
    * 
    * Parameters
    * ----------
    * 
    * `Trajectory` : traj object
    * 
    * `align` : alignments object
    * 
    * `start_pair` : pair of structures to begin with
    */
    int n;
    int i, j;
    double **x, **y, **xlig, **ylig, **ylig_rot;
    double u[3][3], cm1[3], cm2[3];
    int rot_idx;
    // allocating vectors
    x = d2t(Trajectory->n_at, 3);
    xlig = d2t(Ligand_Trajectory->n_at, 3);
    y = d2t(Trajectory->n_at, 3);
    ylig = d2t(Ligand_Trajectory->n_at, 3);
    ylig_rot = d2t(Ligand_Trajectory->n_at, 3);
    int ref = start_pair[0];
    int mod = start_pair[1];
    int old_ref = ref - 1;

    for (n = 0; n < Trajectory->pairs; n++) {
        if (ref != old_ref) {
            // recalculate com and translate structure
            zero_vec_d(cm1, 3);
            for (i = 0; i < Trajectory->n_at; i++) {
                for (j = 0; j < 3; j++) {
                    cm1[j] += Trajectory->traj_coords[ref][i * 3 + j] / Trajectory->n_at;
                }
            }
            // filling structure
            for (i = 0; i < Trajectory->n_at; i++) {
                for (j = 0; j < 3; j++) {
                    x[i][j] = Trajectory->traj_coords[ref][i * 3 + j] - cm1[j];
                }
            }
            // now the ligand
            for (i = 0; i < Ligand_Trajectory->n_at; i++) {
                for (j = 0; j < 3; j++) {
                    xlig[i][j] = Ligand_Trajectory->traj_coords[ref][i * 3 + j] - cm1[j];
                }
            }
            old_ref += 1;
        }
        
        // second structure
        zero_vec_d(cm2, 3);
        for (i = 0; i < Trajectory->n_at; i++) {
                for (j = 0; j < 3; j++) {
                    cm2[j] += Trajectory->traj_coords[mod][i * 3 + j] / Trajectory->n_at;
            }
        }
        for (i = 0; i < Trajectory->n_at; i++) {
            for (j = 0; j < 3; j++) {
                y[i][j] = Trajectory->traj_coords[mod][i * 3 + j] - cm2[j];
            }
        }
        // now ylig
        for (i = 0; i < Ligand_Trajectory->n_at; i++) {
            for (j = 0; j < 3; j++) {
                ylig[i][j] = Ligand_Trajectory->traj_coords[mod][i * 3 + j] - cm2[j];
            }
        }

        // Alignment AND RMSD computation
        align->rmsd_mat[n] = optimal_alignment(x, y, Trajectory->n_at, u);
        align->ref_structs[n] = ref + 1;
        align->mod_structs[n] = mod + 1;
        
        // now correcting the coordinates of the ligand
        for (i = 0; i < Ligand_Trajectory->n_at; i++) {
            for (j = 0; j < 3; j++) {
                ylig_rot[i][j] = 0.0;
                for (rot_idx = 0; rot_idx < 3; rot_idx++) {
                  ylig_rot[i][j] += u[rot_idx][j] * ylig[i][rot_idx];
                }
            }
        }
        // calculate rmsd between ligands
        double ilrmsd = 0.0;
        double delta;
        for (i = 0; i < Ligand_Trajectory->n_at; i++) {
            for (j = 0; j < 3; j++) {
                delta = ylig_rot[i][j] - xlig[i][j];
                ilrmsd += delta * delta;
            }
        }

        ilrmsd = sqrt(ilrmsd / Ligand_Trajectory->n_at);
        // the value before was the receptor interface rmsd, now let's change it to the ligand interface rmsd
        align->rmsd_mat[n] = ilrmsd;

        // update idx
        if (mod == Trajectory->frames - 1) {
            ref += 1;
            mod = ref + 1;
        }
        else {
            mod += 1;
        }
        
    }
    free_d2t(x);
    free_d2t(y);
    free_d2t(xlig);
    free_d2t(ylig);
    free_d2t(ylig_rot);
}

void cycle_rmsd(traj *Trajectory, alignments *align, int start_pair[2]) {
    /**
    * routine that cycles over all pairs of frames in a trajectory, calling `optimal_alignment`
    * 
    * Parameters
    * ----------
    * 
    * `Trajectory` : traj object
    * 
    * `align` : alignments object
    * 
    * `start_pair` : pair of structures to begin with
    */ 
    int n;
    int i, j;
    double **x, **y;
    double u[3][3], cm1[3], cm2[3];
    printf("calculating rotation matrices..\n");
    // allocating vectors
    x = d2t(Trajectory->n_at, 3);
    y = d2t(Trajectory->n_at, 3);
    int ref = start_pair[0];
    int mod = start_pair[1];
    int old_ref = ref - 1;
    //for (n = 0; n < Trajectory->frames; n++) {
    for (n = 0; n < Trajectory->pairs; n++) {
        if (ref != old_ref) {
            // recalculate com and translate structure

            zero_vec_d(cm1, 3);
            for (i = 0; i < Trajectory->n_at; i++) {
                for (j = 0; j < 3; j++) {
                    cm1[j] += Trajectory->traj_coords[ref][i * 3 + j] / Trajectory->n_at;
                }
            }
            // filling structure
            for (i = 0; i < Trajectory->n_at; i++) {
                for (j = 0; j < 3; j++) {
                    x[i][j] = Trajectory->traj_coords[ref][i * 3 + j] - cm1[j];
                }
            }
            old_ref += 1;
        }
        // second structure
        zero_vec_d(cm2, 3);
        for (i = 0; i < Trajectory->n_at; i++) {
                for (j = 0; j < 3; j++) {
                    cm2[j] += Trajectory->traj_coords[mod][i * 3 + j] / Trajectory->n_at;
            }
        }
        for (i = 0; i < Trajectory->n_at; i++) {
            for (j = 0; j < 3; j++) {
                y[i][j] = Trajectory->traj_coords[mod][i * 3 + j] - cm2[j];
            }
        }
        // Alignment AND RMSD computation
        align->rmsd_mat[n] = optimal_alignment(x, y, Trajectory->n_at, u);
        align->ref_structs[n] = ref + 1;
        align->mod_structs[n] = mod + 1;
        
        // update idx
        if (mod == Trajectory->frames - 1) {
            ref += 1;
            mod = ref + 1;
        }
        else {
            mod += 1;
        }
        
    }
    free_d2t(x);
    free_d2t(y);
}

void read_TrajectoryFile(char *TrajFileName, traj *Trajectory){
    /**
    * routine that reads the input xyz coordinate file
    *
    * Parameters
    * ----------
    *
    * `TrajFileName` : trajectory filename
    *
    * `Trajectory` : traj object
    */
 
    FILE *ft;                   // trajectory file; 
    FILE *fe;                   // declaring error file;  

    int rows_i, cols_i, i, p, frame_idx, j, count; 
    size_t line_buf_size = 0;

    char *token;
    char *token_arr[100];  //*string, *line;
    char *string = NULL; 
    char *line;  

    /* Opening and check if it exist and if it is empty */
    ft = fopen(TrajFileName, "r");
    printf("Reading Trajectory FILE %s\n", TrajFileName);

    //check_empty_file(ft, TrajFileName);

    rows_i = n_rows(ft);                                // Number of rows in trajectory file "ft".
    fseek(ft, 0, SEEK_SET);
   
    line   = (char *) malloc (200); 					     // Allocate char line with size (lenght) = 200; 
									     // We are sure that the lenght of each line is less than 200

    /* Initialize the 2D-array Trajectory->traj_coords[][] */
    for(i = 0; i < Trajectory->frames; i++){
        for(j = 0; j < Trajectory->n_at; j++){
            Trajectory->traj_coords[i][3 * j + 0] = 0.0;
            Trajectory->traj_coords[i][3 * j + 1] = 0.0;
            Trajectory->traj_coords[i][3 * j + 2] = 0.0;
        }
    }

    /* Checking for empty rows; 
     * Checking that there are ONLY rows with one column and three columns 
     * Checking that the rows with one column correspond to an integer number i.e. the number of atoms; 
     * Checking that the rows with three columns are float corresponding to x,y,z trajectory coordinate. */

    frame_idx = 0;                                      // Initialize frame index to  0 
    p  = 0;                                             // Initialize counter "p" to 0 

    for(i=0; i<rows_i; i++){
        if(p>=3*Trajectory->n_at)                       // p increases from 0 to 3*atomnum i.e. p=[0;3*230) i.e. p = [0;689]
            p = 0; 

        getline(&string, &line_buf_size, ft);           // Reading entire line;   
   
        strcpy(line, string);                           // Copying "string" in "line"; we need it in case of three columns. 

        //if( i != (Trajectory->n_at + 2)*(frame_idx-1) + 1)  // Checking for empty rows except for the 2nd row of each frame representing the title
        //    check_empty_rows(string);                   // that could also be an empty string   

        token = strtok(string, " \t\v\r\n");            // Splitting the string-line in columns separated by spaces, tabs, or \n  or \v or \r

        cols_i = columns(token);                        // Counting the number of columns in each row of file.  

         
        if(i!=(Trajectory->n_at + 2)*(frame_idx-1) + 1){   // exclude the 2nd row of each frame 

            if(cols_i == 1){

                if(i != (Trajectory->n_at + 2)*frame_idx) {
                    printf("Error. The %dth frame contains a different number of rows. Each frame must have %d rows\n", frame_idx+1, Trajectory->n_at + 2);
                    exit(EXIT_FAILURE);
                } 
                else{
                    frame_idx++;
                    if(frame_idx>Trajectory->frames){
                        fe = fopen("error.dat", "w");
                        fprintf(fe, "Error. The number of trajectory frames is higher than the declared one in parameter file (%d).\nAborting\n", Trajectory->frames);
                    	fclose(fe);
                    	printf("Error. The number of trajectory frames is higher than the declared one in parameter file (%d).\nAborting\n", Trajectory->frames);
                    	exit(EXIT_FAILURE);
                    }
                    //check_int_string(string, i, TrajFileName);             // Checking that each row is an INTEGER number. 
                    if(atoi(string) != Trajectory->n_at){
                        fe = fopen("error.dat","w");
                    	fprintf(fe,"Error. The number of atoms at %dth row has length (%d) different from atomnum(%d). Aborting\n",\
                                i+1,atoi(string),Trajectory->n_at);
                   	 printf("Error. The number of atoms at %dth row has length (%d) different from atomnum(%d). Aborting\n",\
                                     i+1,atoi(string),Trajectory->n_at);
                    	fclose(fe); 
                    	exit(EXIT_FAILURE);
                    }
                }
            }
    
            if(cols_i == 2){
 
                if(i != (Trajectory->n_at + 2)*frame_idx + 1){
                    fe = fopen("error.dat", "w");
                    fprintf(fe,"Error. Each row must not contain 2 columns (except the title in the 2nd row of each frame that can contain N columns).\n\
                      	        ONLY 1 or 4 columns are allowed in Trajectory (%dth row has 2 columns)\n", i+1);

                    printf("Error. Each row must not contain 2 columns (except the title in the 2nd row of each frame that can contain N columns)\n\
                            ONLY 1 or 4 columns are allowed in Trajectory (%dth row has 2 columns.\n", i+1);
                    fclose(fe); 
                    exit(EXIT_FAILURE);
                }
            }

            if(cols_i == 3){
 
                if(i != (Trajectory->n_at + 2)*frame_idx + 1){
                    fe = fopen("error.dat", "w");
                    fprintf(fe,"Error. Each row must not contain 3 columns (except the title in the 2nd row of each frame that can contain N columns). \
                                ONLY 1 or 4 columns are allowed in Trajectory (%dth row has 3 columns)\n", i+1);

                    printf("Error. Each row must not contain 3 columns (except the title in the 2nd row of each frame that can contain N columns) \
                            ONLY 1 or 4 columns are allowed in Trajectory (%dth row has 3 columns.\n", i+1);
                    fclose(fe); 
                    exit(EXIT_FAILURE);
                }
            }

            if(cols_i == 4){ 
                if(i != (Trajectory->n_at + 2)*frame_idx){
                    token_arr[0] = strtok(line, " \t\v\r\n");

                    count = 0;  
                    while(token_arr[count]){
                        count++; 
                        token_arr[count] = strtok(NULL, " \t\v\r\n");
                    }

                    for (j=1; j <= count-1 ; j++){
                        //check_float_string(token_arr[j], i, TrajFileName);                     // Checking that each row is an FLOAT number.
                        Trajectory->traj_coords[frame_idx-1][p] = atof(token_arr[j]);          // Assigning each float-string to traj_coords[i][j] 2D-array 
                        p++;
                    }

                }

                else{ 
                    printf("Error. Each row of the frame (except the 1st row containing the number of atoms and the 2nd row containing the title),\n\
                    must contain 4 columns corresponding to at_name, x, y, and z coodinates. Check also if there is one extra row in %dth frame\n",
                            frame_idx); 
                    exit(EXIT_FAILURE);
                }
            }

            if(cols_i > 4 ){ 
                if(i != (Trajectory->n_at + 2)*frame_idx){ 
                    fe = fopen("error.dat","w"); 
                    fprintf(fe,     "Error. The maximum number of columns allowed is 4. The %dth row of your file contains %d columns\n", i+1, cols_i);
                    printf("Error. The maximum number of columns allowed is 4. The %dth row of your file contains %d columns\n", i+1, cols_i);
                    fclose(fe); 
                    exit(EXIT_FAILURE);
                }
            }
        }

    }  // END FOR LOOP 
    
    free(string);
    //free(token);
    free(line);
    fclose(ft);// Close trajectory file. 
    // 1st check: frame_idx must be = frames
    if (frame_idx !=  Trajectory->frames) {
        fe = fopen("error.dat", "w");
        fprintf(fe, "Error. Frames completed: %d, while declared frames in parameter file are %d. Trajectory incomplete.\nAborting\n",frame_idx,\
                     Trajectory->frames);
        fclose(fe);
        printf("Error. Frames completed: %d, while declared frames in parameter file are %d. Trajectory incomplete.\nAborting\n",frame_idx,\
                Trajectory->frames);
        exit(EXIT_FAILURE);
    }

    // final check: frame_idx corresponds with what we expect but the number of row is not the same => It means that 
    //              the frame completed are frame_idx - 1 
    else{
        if(rows_i !=  (Trajectory->n_at + 2)*Trajectory->frames){
            fe = fopen("error.dat", "w");
            fprintf(fe,"Error. Total number of rows: %d, while expected number of rows are %d. Trajectory incomplete.\nAborting\n",rows_i,\
                        (Trajectory->n_at + 2)*Trajectory->frames);
            fclose(fe);
            printf("Error. Total number of rows: %d, while expected number of rows are %d. Trajectory incomplete.\nAborting\n", rows_i,\
                    (Trajectory->n_at + 2)*Trajectory->frames);
            exit(EXIT_FAILURE);
        }
    }
}


int main(int argc, char *argv[]) {
    /**
    * main file of the program
    */
    time_t seconds;
    time_t seconds_ref = time(NULL);
    printf("Starting program at %s\n", ctime(&seconds_ref));

    printf("input trajectory %s\n", argv[1]);
    printf("core index is %s\n", argv[2]);
    printf("number of pairs is %s\n", argv[3]);
    printf("starting pair is %s-%s\n", argv[4], argv[5]);
    traj *Trajectory = malloc (sizeof(traj));

    int npairs = atoi(argv[3]);
    int start_pair[2];
    start_pair[0] = atoi(argv[4]);
    start_pair[1] = atoi(argv[5]);

    int rec_frames = atoi(argv[6]);
    int rec_atomnum = atoi(argv[7]);
     
    Trajectory->frames = rec_frames;
    Trajectory->n_at = rec_atomnum;
    Trajectory->traj_coords = d2t(rec_frames, 3 * rec_atomnum);
    Trajectory->pairs = npairs;
    
    printf("frames = %d\n", Trajectory->frames);
    printf("overall pairs = %d\n", Trajectory->pairs);
    
    // read trajectory
    read_TrajectoryFile(argv[1], Trajectory);
    alignments *align = malloc (sizeof(alignments));
    align->rmsd_mat = d1t(Trajectory->pairs);
    align->ref_structs = i1t(Trajectory->pairs);
    align->mod_structs = i1t(Trajectory->pairs);
    align->coms = d2t(rec_frames, 3);

    char out_filename[100];
    if (argc == 10){
        // if there are arguments 8 to 10, it is the ligand trajectory, frames and atomnum
        traj *Ligand_Trajectory = malloc (sizeof(traj));
        printf("loading ligand data\n");

        int ligand_atomnum = atoi(argv[9]);
        Ligand_Trajectory->frames = rec_frames;
        Ligand_Trajectory->n_at = ligand_atomnum;
        Ligand_Trajectory->traj_coords = d2t(rec_frames, 3 * ligand_atomnum);
        read_TrajectoryFile(argv[8], Ligand_Trajectory);
        cycle_ilrmsd(Trajectory, Ligand_Trajectory, align, start_pair);
        sprintf(out_filename, "ilrmsd_%s.matrix", argv[2]);
    }
    else{
        cycle_rmsd(Trajectory, align, start_pair);
        sprintf(out_filename, "rmsd_%s.matrix", argv[2]);
    }
    // write the rmsd matrix
    int i;
    // check if file exists
    if (access(out_filename, F_OK) != -1){
        printf("Warning: file %s already exists.\n", out_filename);
    }
    // write to file
    FILE *frmsd;
    frmsd = fopen(out_filename, "w");
    for (i = 0; i < Trajectory->pairs; i++){
        fprintf(frmsd, "%d %d %.3lf\n",align->ref_structs[i],align->mod_structs[i],align->rmsd_mat[i]);
    }
    
    seconds = time(NULL);
    printf("\nOverall execution time: %ld seconds\n", seconds-seconds_ref);
    free_d2t(align->coms);
    free_d1t(align->rmsd_mat);
    free(align);
    free_d2t(Trajectory->traj_coords);
    free(Trajectory);
    return 0;
}
