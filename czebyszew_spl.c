#include <stdlib.h>
#include <math.h>
#include "makespl.h"
#include "piv_ge_solver.h"

double gen_baz(double x, int deg) {
    double km1, k, kp1;
    if (deg == 0)
        return 1;
    if (deg == 1)
        return x;
    k = x;
    km1 = 1;
    for (int i = 2; i <= deg; i++) {
        kp1 = 2 * x * k - km1;
        km1 = k;
        k = kp1;
    }
    return kp1;
}

double getf1(double *wspz, double x, int n) {
    double sum = 0;
    for (int i = 1; i < n; i++)
        sum += gen_baz(x, i - 1) * wspz[i] * i;
    return sum;
}

double getf2(double *wspz, double x, int n) {
    double sum = 0;
    for (int i = 2; i < n; i++)
        sum += gen_baz(x, i - 2) * wspz[i] * i;
    return sum;
}

double getf3(double *wspz, double x, int n) {
    double sum = 0;
    for (int i = 3; i < n; i++)
        sum += gen_baz(x, i - 3) * wspz[i] * i;
    return sum;
}

double getf(double *wsp, double x, int n) {
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += gen_baz(x, i) * wsp[i];
    return sum;
}

double calc_dx(double x, double *wsp,  int n, double s, double t){
    double delda=1.0e-6;
    return (getf(wsp, s + t * (x+delda), n) - getf(wsp, s + t * (x-delda), n)) / delda;
}

double calc_ddx(double x, double *wsp,  int n, double s, double t){
    double delda=1.0e-6;
    return (calc_dx((x+delda),wsp,n,s,t)-calc_dx((x-delda),wsp,n,s,t))/delda;
}

void make_spl(points_t *pts, spline_t *spl) {


    if (pts->n == 0)
        exit(50);

    double max_x = *pts->x, min_x = *pts->x;


    for (int i = 1; i < pts->n; i++) {
        if (pts->x[i] < min_x)
            min_x = pts->x[i];
        if (pts->x[i] > max_x)
            max_x = pts->x[i];

    }
//x*=s+tx
    double t = (max_x - min_x) / 2;
    double s = (max_x + min_x) / 2;

/*
 * Do zapisu:
 * x, y - wartosci punktów
 * f to wynik wielomianu czebyszewa aproksymującego
 * f1-f3 wyniki pochodnych ww. wielomianu
 */


    matrix_t *mat = make_matrix(pts->n, pts->n + 1);
    for (int i = 0; i < pts->n; i++) {//rzad
        int j;
        for (j = 0; j < pts->n; j++) {//kol
            add_to_entry_matrix(mat, i, j, gen_baz(s + t * pts->x[i], j));
        }
        add_to_entry_matrix(mat, i, j, pts->y[i]);
    }
    write_matrix(mat, stdout);
    if (piv_ge_solver(mat)) {
        spl->n = 0;
        return;
    }
    write_matrix(mat, stdout);

    double *wspcz = malloc(pts->n * sizeof *wspcz);
    for (int i = 0; i < pts->n; i++)
        wspcz[i] = get_entry_matrix(mat, i, pts->n);

    alloc_spl(spl, pts->n);
    spl->n = pts->n;
    int i;
    for (i = 0; i < pts->n - 1; i++) {
        double dx = pts->x[i + 1] - pts->x[i];
        spl->x[i] = pts->x[i];
        spl->f[i] = getf(wspcz, s + t * pts->x[i], pts->n);
        double delda=1.0e-6;
        spl->f1[i] = calc_dx(pts->x[i],wspcz,pts->n, s, t);
        spl->f2[i] = calc_ddx(pts->x[i],wspcz,pts->n, s, t);
        //spl->f3[i] = pow((getf(wspcz, s + t * pts->x[i + 1], pts->n) - getf(wspcz, s + t * pts->x[i], pts->n)), 3) /
                     pow(dx, 3);

    }
    spl->x[i] = pts->x[i];
    spl->f[i] = getf(wspcz, s + t * pts->x[i], pts->n);
    spl->f1[i] = spl->f1[i-1];
    spl->f2[i] = spl->f2[i-1];
    spl->f3[i] = spl->f3[i-1];


    printf("----\n");
    FILE *oo=fopen("out.tt","w");
    for(double dt=5.03;dt<5.96;dt+=0.01) {

        fprintf(oo,"%g %g\n", dt, getf(wspcz,s+t*dt,pts->n));
    }


}