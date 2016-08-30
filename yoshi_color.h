#ifndef _YOSHI_COLOR_H_
#define _YOSHI_COLOR_H_
#define Word 100

void read_data (char* data_file, int nech, int Nlam, int Nlines, char column_name[nech+1][Word], double Lambda_data[Nlines], double ME_data[Nlines][nech]);

void read_data_quanty (char* data, int Nlam, int Nlines, char column_name[3][Word], double Lambda_data[Nlines], double ME_data[Nlines][1]);

extern double eps_to_trans (double epsilon);

extern double eV_to_nm(double energy);

extern double yoshi_spectrum(double wave_length, int echantillon, int nlam, int nech, double ME_expected [nlam][nech]);

extern double interpolation (double lambda1, double val1, double lambda2, double val2, double lambda_expected);

extern double match_value (double lambda1, double val1, double lambda2, double val2, double lambda_expected);

struct index_values {    
 int index;
 int size;
 double* values;
};
typedef struct index_values index_values;


extern void print_1Dtable (int l1, double t[l1]);

extern void print_2Dtable (int l1, int l2, double t[l1][l2]);

extern void initialize_Lambda (int nlam, double Lambda[nlam]);

extern bool belongs_interval (double lambda1, double lambda2, double lambda_expected); 

extern struct index_values search_values (int lam_idx, int idx, int nlam, int nech, int nlam_data, double Lambda[nlam], double Lambda_data[nlam_data], double ME_expected [nlam][nech], double ME_data[nlam_data][nech]);

extern void ME_interpolation (int nlam, int nech, int nlam_data, double Lambda[nlam], double Lambda_data[nlam_data], double ME_expected [nlam][nech], double ME_data[nlam_data][nech]);

extern void spectrum_to_xyz(double (*spec_intens)(double wavelength),
                     double *x, double *y, double *z);


extern void spectrum_to_XYZ(int nlam, int nech, double ME_expected[nlam][nech], double *_X_, double *_Y_, double *_Z_, int echantillon);

extern double f (double t);


extern void XYZ_to_lab (double (*f)(double t), double *_X_, double *_Y_,
double *_Z_, double *L, double *a, double *b);

extern void Draw_Node_Graph (char* FILE_OUT, char* x_label, char*
y_label,char* x_range,double Lab[30][3] );

extern void write_Lab (int nech, double Phys[30][3],char* FILE_OUT, char column_name[nech+1][Word]);

extern void write_xyz (int nech, double Phys[30][3],char* FILE_OUT, char column_name[nech+1][Word]);

extern void header ();


extern int mygeti(int *result);

extern int main_yoshi (int nlines, int nech, char* data_file, bool is_quanty);

#endif
