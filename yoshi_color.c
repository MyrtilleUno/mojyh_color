/*
       Yoshi Color
 
 Programme écrit en langage C permettant de calculer les coordonnées Lab d'une couleur
 à partir du spectre de transmission de l'échantillon.
 
 Les données d'entrées doivent se présenter sous forme de colones :
 - 1er colone : longueurs d'onde 
 - colones suivantes : données en transmission (entre 0 et 1)
 
 les spectres doivent couvrir la gamme 380 - 780 nm / 5nm.
 
 le programme calcul les coordinnées Lab avec pour illuminant D65 et observer CIE 1931 2°;
 
 
 2012 - Myrtille Hunault
 
 ********
 Upgrade - 2014 - Myrtille Hunault
 
 - added the sample name in the output files
 - removed the plotting in the main
 - add a 4th argument in write_Lab as the number of samples form the value "Ech" of the main
 - duplicate the write_Lab into write_wyz function to definite columns names in each output file
 
 
 ********

*/

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#define Nech 1
#define Word 100

/*
void get_parameters (char* data_file){

  FILE* file_input = NULL;
  int l, ech;
  file_input = fopen(data_file, "r");
  if (file_input != NULL){}
  else {
    printf("impossible to read file\n");
    return;
  }
  char string1[Word]; char string2[Word];
  double*dbl;
  b = true;
  while (true){
    if b {fscanf(file_input,"%s", string1); b = false;} else {fscanf(file_input,"%s", string2); b = true;}
    if (sscanf(string1, "%f", &d) == 1 && sscanf(string2, "%f", &dbl) == 1) {break;}

   }
   fscanf(file_input,"%s", string2);
     if sscanf(string1, "%f", &nombrefloat) == 1 && 
  }


}

\*

/*---------------------------*/
void read_data (char* data_file, int nech, int Nlam, int Nlines, char column_name[nech+1][Word], double Lambda_data[Nlines], double ME_data[Nlines][nech]){

  FILE* file_input = NULL;
  int l, ech;
  file_input = fopen(data_file, "r");
  if (file_input != NULL){}
  else {
    printf("impossible to read file\n");
    return;
  }


  for (ech = 0; ech < nech +1; ech++){
    fscanf(file_input,"%s", column_name[ech]);
    printf( "%s\n", column_name[ech]) ;
  }
  for (l=0; l< Nlines;l++){
    fscanf(file_input,"%lf", &Lambda_data[l]);
    for (ech = 0; ech < nech ;ech++){
      fscanf(file_input,"%lf", &ME_data[l][ech]);
    }
  }
  fclose(file_input);
  printf("data have been read correctly\n");
  printf("-------------------------------\n");
  /*printf( "%lf\n", ME_data[1][3]) ;*/
}
    /* Old Stuff
       for (c = 0; c< Ech+1 ;c++){
       fscanf(file_input,"%lf",&dummy_ME[c]);
       }
       int current_lambda = dummy_ME[0];
       if (current_lambda>=380 && current_lambda<=780 && (current_lambda % 5)==0){
       int index_ME_data = (current_lambda -380)/5;

       for (c = 0; c< Ech ;c++){
       ME_data[index_ME_data][c] = dummy_ME[c+1];
       ME_expected[index_ME_data][c] = dummy_ME[c+1];
       }
       End of old stuff*/


double yoshi_spectrum(double wave_length, int echantillon, int nlam, int nech, double ME_expected [nlam][nech])
{
    int lambda_index = (wave_length -380.0) / 5;
/*    printf( "%lf\n", ME_data [lambda_index][echantillon]) ;*/
    return ME_expected [lambda_index][echantillon];

}
/*                          SPECTRUM_TO_XYZ

    Calculate the CIE X, Y, and Z coordinates corresponding to
    a light source with spectral distribution given by  the
    function SPEC_INTENS, which is called with a series of
    wavelengths between 380 and 780 nm (the argument is
    expressed in meters), which returns emittance at  that
    wavelength in arbitrary units.  The chromaticity
    coordinates of the spectrum are returned in the x, y, and z
    arguments which respect the identity:

            x + y + z = 1.
*/
/*---------------------------*/

double interpolation (double lambda1, double val1, double lambda2, double val2, double lambda_expected)
{
  return (val2 - val1) / (lambda2 - lambda1) * (lambda_expected - lambda1) + val1;
} 

/* prerequisites : lambda_expected belongs to [lambda1, lambda2] */
double match_value (double lambda1, double val1, double lambda2, double val2, double lambda_expected)
{
  double result = 0.0;
  if (lambda1 == lambda_expected) {result = val1;}
  else if (lambda2 == lambda_expected) {result = val2;}
  else {result = interpolation(lambda1, val1, lambda2, val2, lambda_expected);}
  //printf("%lf\t %lf\t %lf\t %lf\t %lf\t %lf\n", lambda1, lambda2, lambda_expected, val1, val2, result);
  return result;
}

struct index_values {    
 int index;
 double values[Nech];
};

void print_1Dtable (int l1, double t[l1])
{
  int i = 0;
  for (i = 0; i < l1; i++){printf("%lf\n", t[i]);}
  return;
}

void print_2Dtable (int l1, int l2, double t[l1][l2])
{
  int i1 = 0; int i2 = 0;
  for (i1 = 0; i1 < l1; i1++){
    for (i2 = 0; i2 < l2-1; i2++){
      printf("%lf\t", t[i1][i2]);
    }
    printf("%lf\n", t[i1][l2 - 1]);
  }
  return;
}

void initialize_Lambda (int nlam, double Lambda[nlam])
{
  double lambda = 380.0; int lam_idx = 0;
  for (lam_idx = 0; lam_idx < nlam; lam_idx++) {
    Lambda[lam_idx] = lambda;
    lambda += 5.0;
  }
  return;
}

bool belongs_interval (double lambda1, double lambda2, double lambda_expected) {
 return (lambda1 <= lambda_expected &&  lambda_expected <= lambda2) || (lambda2 <= lambda_expected &&  lambda_expected <= lambda1);
}

struct index_values search_values (int lam_idx, int idx, int nlam, int nech, int nlam_data, double Lambda[nlam], double Lambda_data[nlam_data], double ME_expected [nlam][nech], double ME_data[nlam_data][nech]){
 struct index_values r;
 double lambda_expected = Lambda[lam_idx]; int i = 0;
 double lambda1, lambda2 = 0.0;
 for (i = idx; i < nlam_data - 1; i++){
   lambda1 = Lambda_data[i]; lambda2 = Lambda_data[i + 1];
   if (belongs_interval (lambda1, lambda2, lambda_expected)){ 
     int ech = 0;
     for (ech = 0; ech < nech; ech++) {
       double val1 = ME_data[i][ech]/100.0; double val2 = ME_data[i+1][ech]/100.0;
       r.values[ech] = match_value (lambda1, val1, lambda2, val2, lambda_expected);
     }
     //r.index = i + 1; // Hack to accelerate search : only works when data are sorted with increasing lambda
     r.index = 0;
     break;
   }
 }
 return r;
}

void ME_interpolation (int nlam, int nech, int nlam_data, double Lambda[nlam], double Lambda_data[nlam_data], double ME_expected [nlam][nech], double ME_data[nlam_data][nech])
{
  int lam_idx = 0; int idx = 0; double value = 0.0;
  for (lam_idx = 0; lam_idx < nlam; lam_idx++){
    struct index_values couple = search_values (lam_idx, idx, nlam, nech, nlam_data, Lambda, Lambda_data, ME_expected, ME_data);
    idx = couple.index;
    int ech = 0;
    for (ech = 0; ech < nech; ech++) {
      ME_expected [lam_idx][ech] = couple.values[ech];
    }
  }
  return;
}

void spectrum_to_xyz(double (*spec_intens)(double wavelength),
                     double *x, double *y, double *z)
{
    int i;
    double lambda, X = 0, Y = 0, Z = 0, XYZ;

    /* CIE colour matching functions xBar, yBar, and zBar for
       wavelengths from 380 through 780 nanometers, every 5
       nanometers.  For a wavelength lambda in this range:

            cie_colour_match[(lambda - 380) / 5][0] = xBar
            cie_colour_match[(lambda - 380) / 5][1] = yBar
            cie_colour_match[(lambda - 380) / 5][2] = zBar

  To save memory, this table can be declared as floats
  rather than doubles; (IEEE) float has enough
  significant bits to represent the values. It's declared
  as a double here to avoid warnings about "conversion
  between floating-point types" from certain persnickety
  compilers. */

    static double cie_colour_match[81][3] = {
        {0.0014,0.0000,0.0065}, {0.0022,0.0001,0.0105}, {0.0042,0.0001,0.0201},
        {0.0076,0.0002,0.0362}, {0.0143,0.0004,0.0679}, {0.0232,0.0006,0.1102},
        {0.0435,0.0012,0.2074}, {0.0776,0.0022,0.3713}, {0.1344,0.0040,0.6456},
        {0.2148,0.0073,1.0391}, {0.2839,0.0116,1.3856}, {0.3285,0.0168,1.6230},
        {0.3483,0.0230,1.7471}, {0.3481,0.0298,1.7826}, {0.3362,0.0380,1.7721},
        {0.3187,0.0480,1.7441}, {0.2908,0.0600,1.6692}, {0.2511,0.0739,1.5281},
        {0.1954,0.0910,1.2876}, {0.1421,0.1126,1.0419}, {0.0956,0.1390,0.8130},
        {0.0580,0.1693,0.6162}, {0.0320,0.2080,0.4652}, {0.0147,0.2586,0.3533},
        {0.0049,0.3230,0.2720}, {0.0024,0.4073,0.2123}, {0.0093,0.5030,0.1582},
        {0.0291,0.6082,0.1117}, {0.0633,0.7100,0.0782}, {0.1096,0.7932,0.0573},
        {0.1655,0.8620,0.0422}, {0.2257,0.9149,0.0298}, {0.2904,0.9540,0.0203},
        {0.3597,0.9803,0.0134}, {0.4334,0.9950,0.0087}, {0.5121,1.0000,0.0057},
        {0.5945,0.9950,0.0039}, {0.6784,0.9786,0.0027}, {0.7621,0.9520,0.0021},
        {0.8425,0.9154,0.0018}, {0.9163,0.8700,0.0017}, {0.9786,0.8163,0.0014},
        {1.0263,0.7570,0.0011}, {1.0567,0.6949,0.0010}, {1.0622,0.6310,0.0008},
        {1.0456,0.5668,0.0006}, {1.0026,0.5030,0.0003}, {0.9384,0.4412,0.0002},
        {0.8544,0.3810,0.0002}, {0.7514,0.3210,0.0001}, {0.6424,0.2650,0.0000},
        {0.5419,0.2170,0.0000}, {0.4479,0.1750,0.0000}, {0.3608,0.1382,0.0000},
        {0.2835,0.1070,0.0000}, {0.2187,0.0816,0.0000}, {0.1649,0.0610,0.0000},
        {0.1212,0.0446,0.0000}, {0.0874,0.0320,0.0000}, {0.0636,0.0232,0.0000},
        {0.0468,0.0170,0.0000}, {0.0329,0.0119,0.0000}, {0.0227,0.0082,0.0000},
        {0.0158,0.0057,0.0000}, {0.0114,0.0041,0.0000}, {0.0081,0.0029,0.0000},
        {0.0058,0.0021,0.0000}, {0.0041,0.0015,0.0000}, {0.0029,0.0010,0.0000},
        {0.0020,0.0007,0.0000}, {0.0014,0.0005,0.0000}, {0.0010,0.0004,0.0000},
        {0.0007,0.0002,0.0000}, {0.0005,0.0002,0.0000}, {0.0003,0.0001,0.0000},
        {0.0002,0.0001,0.0000}, {0.0002,0.0001,0.0000}, {0.0001,0.0000,0.0000},
        {0.0001,0.0000,0.0000}, {0.0001,0.0000,0.0000}, {0.0000,0.0000,0.0000}
    };

    for (i = 0, lambda = 380; lambda < 780.1; i++, lambda += 5) {
        double Me;

        Me = (*spec_intens)(lambda);
        X += Me * cie_colour_match[i][0];
        Y += Me * cie_colour_match[i][1];
        Z += Me * cie_colour_match[i][2];
    }
    XYZ = (X + Y + Z);
    *x = X / XYZ;
    *y = Y / XYZ;
    *z = Z / XYZ;
}

/* Cette nouvelle fonction différente de la précédente ne rend pas x y
et z tels qu'ils sont définis dans les conventions
de notations mais rend les valeur _X_ etc.  qui sont le rapport : X/Xn */

/*---------------------------*/
void spectrum_to_XYZ(int nlam, int nech, double ME_expected[nlam][nech], double *_X_, double *_Y_, double *_Z_, int echantillon)
{
    int i;
    double lambda, Xn = 0, Yn = 0, Zn = 0,  X = 0, Y = 0, Z = 0, K=0, K2=0 ;

    /* CIE colour matching functions xBar, yBar, and zBar for
       wavelengths from 380 through 780 nanometers, every 5
       nanometers.  For a wavelength lambda in this range:

            cie_colour_match[(lambda - 380) / 5][0] = xBar
            cie_colour_match[(lambda - 380) / 5][1] = yBar
            cie_colour_match[(lambda - 380) / 5][2] = zBar

  To save memory, this table can be declared as floats
  rather than doubles; (IEEE) float has enough
  significant bits to represent the values. It's declared
  as a double here to avoid warnings about "conversion
  between floating-point types" from certain persnickety
  compilers. */


	static double D65[81][1] = {
{49.9755}, {52.3118}, {54.6482}, 
{68.7015}, {82.7549}, {87.1204}, 
{91.486}, {92.4589}, {93.4318}, 
{90.057}, {86.6823}, {95.7736}, 
{104.865}, {110.936}, {117.008}, 
{117.41}, {117.812}, {116.336}, 
{114.861}, {115.392}, {115.923}, 
{112.367}, {108.811}, {109.082}, 
{109.354}, {108.578}, {107.802}, 
{106.296}, {104.79}, {106.239}, 
{107.689}, {106.047}, {104.405}, 
{104.225}, {104.046}, {102.023}, 
{100.000}, {98.1671}, {96.3342}, 
{96.0611}, {95.788}, {92.2368}, 
{88.6856}, {89.3459}, {90.0062},
{89.8026}, {89.5991}, {88.6489},
{87.6987}, {85.4936}, {83.2886},
{83.4939}, {83.6992}, {81.863},
{80.0268}, {80.1207}, {80.2146},
{81.2462}, {82.2778}, {80.281}, 
{78.2842}, {74.0027}, {69.7213},
{70.6652}, {71.6091}, {72.979},
{74.349}, {67.9765}, {61.604}, 
{65.7448}, {69.8856}, {72.4863},
{75.087}, {69.3398}, {63.5927},
{55.0054}, {46.4182}, {56.6118}, 
{66.8054}, {65.0941}, {63.3828}
};




static double cie_colour_match[81][3] = {
{0.001368,0.000039,0.00645}, {0.002236,0.000064,0.01055}, 
{0.004243,0.00012,0.02005}, {0.00765,0.000217,0.03621}, {0.01431,0.000396,0.06785}, 
{0.02319,0.00064,0.1102}, {0.04351,0.00121,0.2074}, {0.07763,0.00218,0.3713}, 
{0.13438,0.004,0.6456}, {0.21477,0.0073,1.03905}, {0.2839,0.0116,1.3856}, 
{0.3285,0.01684,1.62296}, {0.34828,0.023,1.74706}, {0.34806,0.0298,1.7826}, 
{0.3362,0.038,1.77211}, {0.3187,0.048,1.7441}, {0.2908,0.06,1.6692}, 
{0.2511,0.0739,1.5281}, {0.19536,0.09098,1.28764}, {0.1421,0.1126,1.0419}, 
{0.09564,0.13902,0.81295}, {0.05795,0.1693,0.6162}, {0.03201,0.20802,0.46518}, 
{0.0147,0.2586,0.3533}, {0.0049,0.323,0.272}, {0.0024,0.4073,0.2123}, 
{0.0093,0.503,0.1582}, {0.0291,0.6082,0.1117}, {0.06327,0.71,0.07825}, 
{0.1096,0.7932,0.05725}, {0.1655,0.862,0.04216}, {0.22575,0.91485,0.02984}, 
{0.2904,0.954,0.0203}, {0.3597,0.9803,0.0134}, {0.43345,0.99495,0.00875}, 
{0.51205,1.0,0.00575}, {0.5945,0.995,0.0039}, {0.6784,0.9786,0.00275}, 
{0.7621,0.952,0.0021}, {0.8425,0.9154,0.0018}, {0.9163,0.87,0.00165}, 
{0.9786,0.8163,0.0014}, {1.0263,0.757,0.0011}, {1.0567,0.6949,0.001}, 
{1.0622,0.631,0.0008}, {1.0456,0.5668,0.0006}, {1.0026,0.503,0.00034}, {0.9384,0.4412,0.00024}, 
{0.85445,0.381,0.00019}, {0.7514,0.321,0.0001}, {0.6424,0.265,0.00005}, 
{0.5419,0.217,0.00003}, {0.4479,0.175,0.00002}, {0.3608,0.1382,0.00001}, 
{0.2835,0.107,0.0}, {0.2187,0.0816,0.0}, {0.1649,0.061,0.0}, 
{0.1212,0.04458,0.0}, {0.0874,0.032,0.0}, {0.0636,0.0232,0.0}, 
{0.04677,0.017,0.0}, {0.0329,0.01192,0.0}, {0.0227,0.00821,0.0}, 
{0.01584,0.005723,0.0}, {0.011359,0.004102,0.0}, {0.008111,0.002929,0.0}, 
{0.00579,0.002091,0.0}, {0.004109,0.001484,0.0}, {0.002899,0.001047,0.0}, 
{0.002049,0.00074,0.0}, {0.00144,0.00052,0.0}, {0.001,0.000361,0.0}, 
{0.00069,0.000249,0.0}, {0.000476,0.000172,0.0}, {0.000332,0.00012,0.0}, 
{0.000235,0.000085,0.0}, {0.000166,0.00006,0.0}, {0.000117,0.000042,0.0}, 
{0.000083,0.00003,0.0}, {0.000059,0.000021,0.0}, {0.000042,0.000015,0.0}
};

    for (i = 0; i < nlam; i++) {
        double Me;

        Me = ME_expected[i][echantillon];
        X += 5* Me * cie_colour_match[i][0] * D65[i][0] ;
        Y += 5* Me * cie_colour_match[i][1] * D65[i][0];
        Z += 5* Me * cie_colour_match[i][2] * D65[i][0];
		K += 5* cie_colour_match[i][1] * D65[i][0];
    }
  
    /*for (i = 0, lambda = 380; lambda < 780.1; i++, lambda += 5) {

        Xn += cie_colour_match[i][0];
        Yn += cie_colour_match[i][1];
        Zn += cie_colour_match[i][2];
    }*/
 //XYZ = (X + Y + Z);
K2   = 100 / K;
*_X_ = K2 * X / 95.0429669;
*_Y_ = K2 * Y / 100;
*_Z_ = K2 * Z / 108.88005470;
printf("%f\n", K2*X);
printf("%f\n", K2*Y);
printf("%f\n", K2*Z);
   // *_X_ = X / XYZ;
   // *_Y_ = Y / XYZ;
   // *_Z_ = Z / XYZ;
}
/* Fonction f qui permet de convertir XYZ vers LAB */
/*---------------------------*/
double f (double t )
{
double borne=0, mario =0;
borne = pow(6.0/29.0,3.0);

 if (t > borne) {
  return  mario = pow(t,(1.0/3.0));
    }
  else {
  return  mario = (1.0/3.0) * pow((29.0/6.0),2.0) * t + (4.0/29.0) ;
  }
}




/*---------------------------*/

/*                  XYZ_to_LAB
 Fonction qui calcule les coordonnées La*b* à partir de XYZ */

 void XYZ_to_lab (double (*f)(double t), double *_X_, double *_Y_,
double *_Z_, double *L, double *a, double *b)
{
  double Yn, Xn, Zn;
//  Yn = 100.000 ;
  //Xn =  95.047;
 // Zn =108.83;
  Xn=Yn=Zn=1.0;
  *L = 116.0 * (*f)(1.00 * ( *_Y_ /Yn) ) - 16.0;
  *a = 500.0 * (( *f)( *_X_ /Xn ) - (*f)( *_Y_ /Yn ));
  *b = 200.0 * (( *f)( *_Y_/Yn  )- (*f)( *_Z_/Zn ));
}

/*---------------------------*/
void Draw_Node_Graph (char* FILE_OUT, char* x_label, char*
y_label,char* x_range,double Lab[30][3] ){


//printf("Drawing \"random load.png\"\n");
FILE *cmd=popen("gnuplot","w");
char tmp[256]={0x0};
/*fprintf(cmd,"set term png font
\"/usr/share/fonts/truetype/ttf-liberation/LiberationSerif-Regular.ttf,20\"\n");*/
fprintf(cmd,"set terminal svg \n");
fprintf(cmd,"set output \"Lab_B13.svg\"\n");
fprintf(cmd,"set xlabel \'%s\'\n",x_label);
fprintf(cmd,"set ylabel \'%s\'\n",y_label);
int label;
for (label=0;label<30;label++){
fprintf(cmd,"set label %d \"%d\" at %le,%le, 0 left norotate back nopoint offset character 0, 0, 0\n",label+1,label+1,Lab[label][1],Lab[label][2]);

}
fprintf(cmd,"set xrange %s\n",x_range);
fprintf(cmd,"set size 1.0,1.0\n");
fprintf(cmd,"set pointsize 1\n");
fprintf(cmd,"plot \'%s\' using 1:2 notitle\n",FILE_OUT);
fprintf(cmd,"exit \n");
/*
while (fgets(tmp,sizeof(tmp),cmd)!=NULL)
{
    printf("%s\n", tmp);
}*/
pclose(cmd);
}

/*---------------------------*/
void write_Lab (int nech, double Phys[30][3],char* FILE_OUT, char column_name[nech+1][Word]) {

  int j=0;
  FILE* file_out = NULL;
  file_out = fopen(FILE_OUT, "w");
	 fprintf(file_out,"sample_name L a b\n");
  for (j=0;j<nech;j++) {

    //fprintf(file_out,"%.4le %20.15le\n",freq(j),Fourier_Phys[j]);
    //fprintf(file_out,"%d %20.15le\n",j,Fourier_Phys[j]);
    fprintf(file_out,"%s %.4lf %.4lf %.4lf\n",column_name[j+1],Phys[j][0],Phys[j][1],Phys[j][2]);
    //printf("freq = %le\n",freq(j));
    //printf("freq = %le\n",2.0*PI*j/N);
  }

  fclose(file_out);

}

/*---------------------------*/
void write_xyz (int nech, double Phys[30][3],char* FILE_OUT, char column_name[nech+1][Word]) {
	
	int j=0;
	FILE* file_out = NULL;
	file_out = fopen(FILE_OUT, "w");
	fprintf(file_out,"sample_name x y z\n");
	for (j=0;j<nech;j++) {
		
		//fprintf(file_out,"%.4le %20.15le\n",freq(j),Fourier_Phys[j]);
		//fprintf(file_out,"%d %20.15le\n",j,Fourier_Phys[j]);
		fprintf(file_out,"%s %.4lf %.4lf %.4lf\n",column_name[j+1],Phys[j][0],Phys[j][1],Phys[j][2]);
		//printf("freq = %le\n",freq(j));
		//printf("freq = %le\n",2.0*PI*j/N);
	}
	
	fclose(file_out);
	
}
/*---------------------------*/
void header () {
	
	printf("                          @@@  @@@                 \n");
	printf("                          @@@  @@@                 \n");
	printf("                        @@...@@...@@               \n");
	printf("                        @@...  .....@              \n");
	printf("                       @....      ..@              \n");
	printf("                 @@@@@@#..` @@  @@  @              \n");
	printf("                 @#####@..  @@  @@  @              \n");
	printf("                @#'''''@    @@  @@  @              \n");
	printf("                @@'''#@@            @@@@@          \n");
	printf("                @@'''@@@            @@@@@          \n");
	printf("               `@@@@@.....     .  @@.....@         \n");
	printf("              @@@@@@@.....    `.` @@.....@         \n");
	printf("              @@''@................,......@@       \n");
	printf("             @''''@................,.@@..@@@@@     \n");
	printf("             @''''@..................@@..@@@@@     \n");
	printf("             @''''@.........................@@     \n");
	printf("              @#@@@.......    ..............@@     \n");
	printf("              @@@@@.....`      .............@@     \n");
	printf("                @@@..,..       @@...........@@     \n");
	printf("                @@'@@...       @@...........@@     \n");
	printf("                @#'@@...     @@@@...........@@     \n");
	printf("                  @@#@..`    @@@@.......,...@@     \n");
	printf("                  @#@@....     @@.........@@       \n");
	printf("                  @'#@.....     @@.....@@@         \n");
	printf("             @@@@@@@@.,.....    @@....,@@@    \n");     
	printf("             @@@@@@@@.......   @@@@@@@@@@  \n");        
	printf("           @#'''@@@@@......`   @@@@@@@@@@          \n");
	printf("           ##'''@@@@@..,...    @                   \n");
	printf("         @@'''#@@@...@@...     @                   \n");
	printf("    @@@@@  ###......@...`      @                   \n");
	printf("    @@@@@  @@@.,` ` @...       @                   \n");
	printf("    @@...@@@@...    .@@        @                   \n");
	printf("    @@...@@@@...    @@@        @                   \n");
	printf("    @@ ..@@.....    @@@        @                   \n");
	printf("    @@ ..@@.........,@@        @                   \n");
	printf("    @@   @@........@@        @@                    \n");
	printf("      @@ @@........@@      @@                      \n");
	printf("      @@   @@@@@@@@        @@                      \n");
	printf("      @@   @@@@@@@@       @                        \n");
	printf("        @@@@@@@@@   `     @                        \n");
	printf("        @@@@@@@@@         @                        \n");
	printf("        @+++''''@   @@@@@@                         \n");
	printf("        @++'''''@   @#@@@@                         \n");
	printf("        @++'''''@@@@#'@@@#@                        \n");
	printf("      @@++++'`'''#@++''''` #@                      \n");
	printf("      @@++++' '''#@++''''` #@                      \n");
	printf("      @@+++++''''#@+++'''''#@                      \n");
	printf("      @@+++++''''#@+++'''''@@                      \n");
	printf("      @@@@@@@@@@@@@@@@@@@@@@@    \n");    	
	
}

/*---------------------------*/

/*---------------------------*/
int mygeti(int *result)
{
	char buff [ 13 ]; /* signed 32-bit value, extra room for '\n' and '\0' */
	return fgets(buff, sizeof buff, stdin) && sscanf(buff, "%d", result) == 1;
}


/*---------------------------*/
int main(int argc, char* argv[])
{ double Lab[30][3];
 double XYZ[30][3];

   printf("-------------------------------\n");
   printf("  STARTING PROGRAM\n");
	  printf("-------------------------------\n");
	header();
   printf("-------------------------------\n");
   printf("This program calculates xyz and L*a*b* CIE values from a transmission spectrum \n");
   printf("The input spectra should be transmission data (between 0 and 1)  \n");
   printf("covering 380nm - 780nm range by 5nm step\n");
   printf("Data should be presented in column format with the first column containing wavelength values in (nm) \n");
   printf("D65 illuminant and CIE 1931 2° observer will be used\n");
   printf("-------------------------------\n");
/*printf("File is: %s\n",argv[1]);*/
  int Ech;
  do {
	fputs("Enter the number of spectra in your input file: ", stdout);
	fflush(stdout);
	} while ( !mygeti(&Ech) );
  printf("number of spectra = %d\n", Ech);

  int Nlam = 81; int Nlines = 2048;
  char column_name[Ech+1][Word];
  double Lambda [Nlam];
  double Lambda_data [Nlines];
  double ME_data [Nlines][Ech];
  double ME_expected [Nlam][Ech];

  initialize_Lambda (Nlam, Lambda);
  read_data(argv[1],Ech, Nlam, Nlines, column_name, Lambda_data, ME_data);
  print_1Dtable(Nlines, Lambda_data);
  print_2Dtable(Nlines, Ech, ME_data);
	ME_interpolation (Nlam, Ech, Nlines, Lambda, Lambda_data, ME_expected, ME_data);

  /*  printf("  %5.5lf x  %5.5lf y   %5.5lf z\n",_X_,_Y_,_Z_);*/
  int i =0;
  for (i=0; i<Ech; i++){
      double _X_, _Y_, _Z_, L, a,b;
	  spectrum_to_XYZ(Nlam, Ech, ME_expected, &_X_, &_Y_, &_Z_, i);
	  XYZ_to_lab ((*f), &_X_, &_Y_, &_Z_, &L, &a, &b);
	  Lab[i][0]= L;
	  Lab[i][1]= a;
	  Lab[i][2]= b;
	  XYZ[i][0]= _X_;
	  XYZ[i][1]= _Y_;
	  XYZ[i][2]= _Z_;
	  printf("  %5.5lf _X_  %5.5lf _Y_   %5.5lf _Z_ \n",_X_,_Y_,_Z_);
	  printf("  %5.5lf L  %5.5lf a   %5.5lf b\n",L,a,b);
    }
  for (i=Ech; i<29; i++){
	  Lab[i+1][0]= 0;
	  Lab[i+1][1]= 0;
	  Lab[i+1][2]= 0;
  }
	printf("-------------------------------\n");
	printf("  Generate 'lab_data' file containing L*a*b* values ::: DONE \n");
	write_Lab(Ech, Lab, "Lab_data", column_name);
	printf("  Generate 'xyz_data' file containing xyz values  ::: DONE \n");
	write_xyz(Ech, XYZ, "xyz_data", column_name);
	/*Draw_Node_Graph("Lab_data", "a", "b", "[:]",Lab);
/*Draw_Node_Graph("Lab_data_B13", "L", "a", "[:]",Lab);
Draw_Node_Graph("Lab_data_B13", "L", "b", "[:]",Lab);*/
printf("-------------------------------\n");
printf("   NORMAL END\n");
printf("-------------------------------\n");
return 0;
}
