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
        int dodane = 0;
        /*
         * Szukamy w *w elementu którego stopień jest taki sam (w->list[i].st) jak w nowym
         * Jeżeli jest to sumujemy w *w współcz.
         * Jeżeli nie to dodajmey nowy element x_t
         */
        for (int j = 0; j < w->size; j++) {
            if (w->list[j].st == nowy->list[i].st) {
                w->list[j].wsp += nowy->list[i].wsp;
                dodane = 1;
            }
        }


        if (dodane == 0) {
            w->list = realloc(w->list, (w->size + 1) * sizeof(*w->list));
            w->size++;
            w->list[w->size - 1].st = nowy->list[i].st;
            w->list[w->size - 1].wsp = nowy->list[i].wsp;
        }
    }


}

wielomian_t *gen_b(int n) {
    wielomian_t *w = malloc(sizeof *w);
    wielomian_t *k = malloc(sizeof *k);
    wielomian_t *km1 = malloc(sizeof *km1);
    wielomian_t *km2 = malloc(sizeof *km2);
    w->size = km1->size = km2->size = k->size = 0;
    w->list = malloc(sizeof *w->list);
    k->list = malloc(sizeof *k->list);
    km1->list = malloc(sizeof *km1->list);
    km2->list = malloc(sizeof *km2->list);

        km2->list = malloc(sizeof *km2->list);
        km2->size = 1;
        km2->list->wsp = 1;
        km2->list->st = 0;



    km1->list = malloc(sizeof *km1->list);
    km1->size = 1;
    km1->list->wsp = 1;
    km1->list->st = 1;

    if(n==0)
        return km2;
    if(n==1)
        return km1;
    for (int i = 2; i <= n; i++) {

            for (int j = 0; j < km1->size; j++) {
                km1->list[j].st++;
                km1->list[j].wsp *= 2.;
            }
            for (int j = 0; j < km2->size; j++) {
                km2->list[j].wsp *= -1.;
            }
            wstaw_w(k, km1);
            wstaw_w(k, km2);
            wstaw_w(w, k);
            for (int j = 0; j < km1->size; j++) {
                km1->list[j].st--;
                km1->list[j].wsp /= 2.;
            }
            for (int j = 0; j < km2->size; j++) {
                km2->list[j].wsp *= -1.;
            }
            km2 = km1;
            km1 = k;
            k = malloc(sizeof *k);
            k->list = malloc(sizeof *k->list);
            k->size = 0;

    }
    return w;
}


wielomian_t *czeb(int n, double *wsp) {
    wielomian_t *w = malloc(sizeof *w);
    wielomian_t *nowy;

    w->size = 0;
    w->list = malloc(sizeof *w->list);

    for (int i = 0; i < n; i++) {
        nowy = gen_b(i);
        for (int j = 0; j < nowy->size; j++) {
            nowy->list[j].wsp *= wsp[i];
        }
        wstaw_w(w, nowy);

    }

    return w;
}


wielomian_t *calc_dx(wielomian_t *w) {
    wielomian_t *nowy = malloc(sizeof *nowy);
    int n = w->size;
    nowy->size = w->size;
    nowy->list = malloc(nowy->size * sizeof *nowy->list);
    for (int i = 0; i < n; i++) {

        if (w->list[i].st == 0) {

            nowy->list[i].wsp = 0;
            nowy->list[i].st = 0;
        } else {
            nowy->list[i].wsp = (w->list[i].st) * (w->list[i].wsp);
            nowy->list[i].st = w->list[i].st - 1;

        }

    }
    return nowy;

}


double wylicz(wielomian_t *w, double x) {
    double wynik = 0;
    for (int i = 0; i < w->size; i++) {
        wynik += pow(x, w->list[i].st) * (w->list[i].wsp);
    }
    return wynik;

}

double mapt(double a, double b, double x){ //t(x)
    return ((a+b)/2)+(((b-a)/2)*x);
}

double mapx(double a, double b, double t){ //x(t)
    return (2/(b-a))*t-((a+b)/(b-a));
}


void make_spl(points_t *pts, spline_t *spl) {


    if (pts->n == 0)
        exit(50);

    double a = *pts->x, b = *pts->x;


    for (int i = 1; i < pts->n; i++) {
        if (pts->x[i] < b)
            b = pts->x[i];
        if (pts->x[i] > a)
            a = pts->x[i];

    }
double t=mapx(a,b,a);
    double d=mapx(a,b,b);

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
            wielomian_t *gb=gen_b(j);
            double val=wylicz(gb,pts->x[i]);
            add_to_entry_matrix(mat, i, j, val );
        }
        add_to_entry_matrix(mat, i, j, pts->y[i]);
    }
    write_matrix(mat, stdout);
    if (piv_ge_solver(mat)) {
        spl->n = 0;
        return;
    }
    write_matrix(mat, stdout);

    double *wsp = malloc(pts->n * sizeof *wsp);
    for (int i = 0; i < pts->n; i++)
        wsp[i] = get_entry_matrix(mat, i, pts->n);


    wielomian_t *interpol = czeb(pts->n, wsp);
wielomian_t *dint= calc_dx(interpol);
wielomian_t *ddint= calc_dx(dint);
wielomian_t *dddint= calc_dx(ddint);
    //będzie trzeba elementy skalować ponownie

    alloc_spl(spl, pts->n);
    spl->n = pts->n;
    int i;
    for (i = 0; i < pts->n; i++) {
        double w = get_entry_matrix(mat, i, pts->n);
        spl->x[i] = pts->x[i];
        spl->f[i] = pts->y[i];
        spl->f1[i] =wylicz(dint,pts->x[i]);
        spl->f2[i] =wylicz(ddint,pts->x[i]);
        spl->f3[i] = wylicz(dddint,pts->x[i]);


    }



}