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


double getf(double *wsp, double x, int n) {
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += gen_baz(x, i) * wsp[i];
    return sum;
}

typedef struct x {
    double wsp;
    int st;
} x_t;

typedef struct w {
    x_t *list;
    int size;
} wielomian_t;


void wstaw_w(wielomian_t *w, wielomian_t *nowy) {

    for (int i = 0; i < nowy->size; i++) {
        int dodane=0;
        /*
         * Szukamy w *w elementu którego stopień jest taki sam (w->list[i].st) jak w nowym
         * Jeżeli jest to sumujemy w *w współcz.
         * Jeżeli nie to dodajmey nowy element x_t
         */
        for(int j=0;j<w->size;j++){
            if (w->list[j].st == nowy->list[i].st)
            {
                w->list[j].wsp += nowy->list[i].wsp;
                dodane=1;
            }
        }


        if(dodane==0){
            w->list=realloc(w->list,(w->size +1)*sizeof(*w->list));
            w->size++;
            w->list[w->size-1].st=nowy->list[i].st;
            w->list[w->size-1].wsp=nowy->list[i].wsp;
        }
    }


}

wielomian_t *gen_b(int n) {
    wielomian_t *w = malloc(sizeof *w);
    wielomian_t *k = malloc(sizeof *k);
    wielomian_t *km1 = malloc(sizeof *km1);
    wielomian_t *km2 = malloc(sizeof *km2);
    w->size = km1->size = km2->size= k->size = 0;
    w->list= malloc(sizeof *w->list);
    k->list= malloc(sizeof *k->list);
    km1->list= malloc(sizeof *km1->list);
    km2->list= malloc(sizeof *km2->list);
    for (int i = 0; i <= n; i++) {
        if (i == 0) {
            km2->list = malloc(sizeof *km2->list);
            km2->size = 1;
            km2->list->wsp = 1;
            km2->list->st = 0;
            wstaw_w(w, km2);

        } else if (i == 1) {
            km1->list = malloc(sizeof *km1->list);
            km1->size = 1;
            km1->list->wsp = 1;
            km1->list->st = 1;
            wstaw_w(w, km1);
        } else {
            for (int j = 0; j < km1->size; j++) {
                km1->list[j].st++;
                km1->list[j].wsp *= 2.;
            }
            for (int j = 0; j < km2->size; j++) {
                km2->list[j].wsp *= -1.;
            }
            wstaw_w(k, km1);
            wstaw_w(k,km2);
            wstaw_w(w,k);
            for (int j = 0; j < km1->size; j++) {
                km1->list[j].st--;
                km1->list[j].wsp /= 2.;
            }
            for (int j = 0; j < km2->size; j++) {
                km2->list[j].wsp *= -1.;
            }
            km2=km1;
            km1=k;
            k=malloc(sizeof *k);
            k->list= malloc(sizeof *k->list);
            k->size=0;
        }
    }
    return w;
}

wielomian_t *gen_b(int n) {
    wielomian_t *w = malloc(sizeof *w);
    wielomian_t *k = malloc(sizeof *k);
    wielomian_t *km1 = malloc(sizeof *km1);
    wielomian_t *km2 = malloc(sizeof *km2);
    w->size = km1->size = km2->size= k->size = 0;
    w->list= malloc(sizeof *w->list);
    k->list= malloc(sizeof *k->list);
    km1->list= malloc(sizeof *km1->list);
    km2->list= malloc(sizeof *km2->list);
    for (int i = 0; i <= n; i++) {
        if (i == 0) {
            km2->list = malloc(sizeof *km2->list);
            km2->size = 1;
            km2->list->wsp = 1;
            km2->list->st = 0;
            wstaw_w(w, km2);

        } else if (i == 1) {
            km1->list = malloc(sizeof *km1->list);
            km1->size = 1;
            km1->list->wsp = 1;
            km1->list->st = 1;
            wstaw_w(w, km1);
        } else {
            for (int j = 0; j < km1->size; j++) {
                km1->list[j].st++;
                km1->list[j].wsp *= 2.;
            }
            for (int j = 0; j < km2->size; j++) {
                km2->list[j].wsp *= -1.;
            }
            wstaw_w(k, km1);
            wstaw_w(k,km2);
            wstaw_w(w,k);
            for (int j = 0; j < km1->size; j++) {
                km1->list[j].st--;
                km1->list[j].wsp /= 2.;
            }
            for (int j = 0; j < km2->size; j++) {
                km2->list[j].wsp *= -1.;
            }
            km2=km1;
            km1=k;
            k=malloc(sizeof *k);
            k->list= malloc(sizeof *k->list);
            k->size=0;
        }
    }
    return w;
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
    double t = 2 / (max_x - min_x);
    double s = -(min_x + max_x) / (max_x - min_x);

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
    //write_matrix(mat, stdout);
    if (piv_ge_solver(mat)) {
        spl->n = 0;
        return;
    }
    //write_matrix(mat, stdout);

    double *wsp = malloc(pts->n * sizeof *wsp);
    for (int i = 0; i < pts->n; i++)
        wsp[i] = get_entry_matrix(mat, i, pts->n);


    alloc_spl(spl, pts->n);
    spl->n = pts->n;
    int i;
    puts("Churski\t|\tNasze");
    for (i = 0; i < pts->n - 1; i++) {
        //printf("%g\t|\t%g\n", fi(pts->n, s + t * pts->x[i]), getf(wsp, s + t * pts->x[i], pts->n));


        double w = get_entry_matrix(mat, i, pts->n);
        spl->x[i] = pts->x[i];
        spl->f[i] = gen_baz(s + t * pts->x[i], pts->n);
        //spl->f1[i] = w*dfi(pts->n,s+t*pts->x[i]);
        //spl->f2[i] = w*d2fi(pts->n,s+t*pts->x[i]);
        //spl->f3[i] = w*d3fi(pts->n,s+t*pts->x[i]);


    }
    spl->x[i] = pts->x[i];
    spl->f[i] = getf(wsp, s + t * pts->x[i], pts->n);
    spl->f1[i] = spl->f1[i - 1];
    spl->f2[i] = spl->f2[i - 1];
    spl->f3[i] = spl->f3[i - 1];


    printf("----\n");
    FILE *oo = fopen("out.tt", "w");
    for (double dt = 5.03; dt < 5.96; dt += 0.01) {

        fprintf(oo, "%g %g\n", dt, getf(wsp, s + t * dt, pts->n));
    }


}