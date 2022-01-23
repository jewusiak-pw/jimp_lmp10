#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "piv_ge_solver.h"
#include "makespl.h"


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

    if (n == 0)
        return km2;
    if (n == 1)
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


/* UWAGA: liczbę używanych f. bazowych można ustawić przez wartość
          zmiennej środowiskowej APPROX_BASE_SIZE
*/

/*
 * Funkcje bazowe: n - liczba funkcji a,b - granice przedzialu aproksymacji i
 * - numer funkcji x - wspolrzedna dla ktorej obliczana jest wartosc funkcji
 */
double
fi(wielomian_t *w, double x) {
    return wylicz(w, x);
}

/* Pierwsza pochodna fi */
double
dfi(wielomian_t *w, double x) {
    return wylicz(calc_dx(w), x);
}

/* Druga pochodna fi */
double
d2fi(wielomian_t *w, double x) {
    return wylicz(calc_dx(calc_dx(w)), x);
}

/* Trzecia pochodna fi */
double
d3fi(wielomian_t *w, double x) {
    return wylicz(calc_dx(calc_dx(calc_dx(w))), x);
}

/* Pomocnicza f. do rysowania bazy */
double
xfi(double a, double b, int n, int i, FILE *out) {
    double h = (b - a) / (n - 1);
    double h3 = h * h * h;
    int hi[5] = {i - 2, i - 1, i, i + 1, i + 2};
    double hx[5];
    int j;

    for (j = 0; j < 5; j++)
        hx[j] = a + h * hi[j];

    fprintf(out, "# nb=%d, i=%d: hi=[", n, i);
    for (j = 0; j < 5; j++)
        fprintf(out, " %d", hi[j]);
    fprintf(out, "] hx=[");
    for (j = 0; j < 5; j++)
        fprintf(out, " %g", hx[j]);
    fprintf(out, "]\n");
}

double mapx(double a, double b, double t) { //x(t)
    return (2 / (b - a)) * t - ((a + b) / (b - a));
}

void
make_spl(points_t *pts, spline_t *spl) {

    matrix_t *eqs = NULL;
    double *x = pts->x;
    double *y = pts->y;
    double a = x[0];
    double b = x[pts->n - 1];
    int i, j, k;
    int nb = pts->n  > 10 ? 10: pts->n ;
    char *nbEnv = getenv("APPROX_BASE_SIZE");

    if (nbEnv != NULL && atoi(nbEnv) > 0)
        nb = atoi(nbEnv);

    eqs = make_matrix(nb, nb + 1);



    for (j = 0; j < nb; j++) {
        for (i = 0; i < nb; i++)
            for (k = 0; k < pts->n; k++)
                add_to_entry_matrix(eqs, j, i, wylicz(gen_b(i), mapx(a, b, x[k])) * wylicz(gen_b(j), mapx(a, b, x[k])));

        for (k = 0; k < pts->n; k++)
            add_to_entry_matrix(eqs, j, nb, y[k] * wylicz(gen_b(j), mapx(a,b,x[k])));
    }


    if (piv_ge_solver(eqs)) {
        spl->n = 0;
        return;
    }

    if (alloc_spl(spl, nb) == 0) {
        for (i = 0; i < spl->n; i++) {
            double xx = spl->x[i] = a + i * (b - a) / (spl->n - 1);
            xx += 10.0 * DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
            spl->f[i] = 0;
            spl->f1[i] = 0;
            spl->f2[i] = 0;
            spl->f3[i] = 0;
            for (k = 0; k < nb; k++) {
                double ck = get_entry_matrix(eqs, k, nb);
                wielomian_t *bz = gen_b(k);
                xx=mapx(a,b,xx);
                spl->f[i] += ck * fi(bz, xx);
                spl->f1[i] += ck * dfi(bz, xx);
                spl->f2[i] += ck * d2fi(bz, xx);
                spl->f3[i] += ck * d3fi(bz, xx);
            }
        }
    }



}
