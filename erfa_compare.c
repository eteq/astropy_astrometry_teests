// Compilation instructions:
// First set up an environment variable $ERFA_SRC pointing to the ERFA source  
// code (the src directory, not the erfa dir itself).  Then do:
// gcc -I $ERFA_SRC -o erfa_compare erfa_compare.c `find $ERFA_SRC -name "*.c" ! -name t_erfa_c.c` -lm
// The resulting executable takes a ra/dec list in stdin and outputs the altaz:
// ./erfa_compare < radec_icrs.dat > altaz_erfa.dat

#include <stdio.h>
#include <math.h>
#include "erfa.h"

double jd_2010 = 2455197.5;
double dut1_2010 = 0.1140749;

//PMs in radians
double pm_x = 4.786371548122005e-07;
double pm_y = 9.347644104104853e-07;

double lon = -95. * M_PI / 180.;
double lat = 48. * M_PI / 180.;


int main() {
    FILE * fp;
    double radeg, decdeg, azob, zenob, hob, dob, rob, eo, altob;

    //fp = fopen("radec_icrs.dat", "r");
    fp = stdin;
    if (fp == NULL ) perror("Error opening file");
    else {
        while (fscanf(fp, " %lf %lf ", &radeg, &decdeg) != EOF) {
            eraAtco13(radeg*M_PI/180., decdeg*M_PI /180., 
              0, 0, 0, 0, 
              jd_2010, 0, dut1_2010,
              lon, lat, 0, pm_x, pm_y,
              0, 0, 0, 1,
              &azob, &zenob, &hob, &dob, &rob, &eo);
            altob = M_PI/2. - zenob;
            printf("%.10f %.10f\n", altob*180/M_PI, azob*180/M_PI);
        }

        fclose(fp);
    }
    return 0;
}

