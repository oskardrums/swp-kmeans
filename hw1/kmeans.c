#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct km_ctx_s {
    size_t   dimension; /* dimension                                           */
    size_t   nclusters; /* # of clusters                                       */
    size_t   nobserves; /* total # of observations                             */
    size_t   max_iters; /* MAX_ITER                                            */
    size_t   inits_set; /* number initial mean values inserted so far          */
    size_t * ob_clusts; /* array of size N, observation index -> cluster index */
    size_t * cardinals; /* cardinalities of clusters, array of size K          */
    double * mean_vals; /* mean values, array of K * d                         */
    double * data_vals; /* input values, array of N * d                        */
};

void km_dump(struct km_ctx_s * ctx)
{
    size_t i, j;
    for (i = 0; i < ctx->nclusters; ++i) {
        for (j = 0; j < ctx->dimension - 1; ++j) {
            printf("%.2f,", *(ctx->mean_vals + (ctx->dimension * i) + j));
        }
        printf("%.2f\n", *(ctx->mean_vals + (ctx->dimension * i) + j));
    }
}

void km_destroy(struct km_ctx_s * ctx)
{
    if (ctx != NULL) {
        if (ctx->cardinals != NULL) {
            free(ctx->cardinals);
        }
        if (ctx->mean_vals != NULL) {
            free(ctx->mean_vals);
        }
        if (ctx->data_vals != NULL) {
            free(ctx->data_vals);
        }
        free(ctx);
    }
}

struct km_ctx_s * km_create(size_t d, size_t k, size_t n, size_t m)
{
    struct km_ctx_s * ctx = (struct km_ctx_s *) malloc (sizeof(struct km_ctx_s));
    if (ctx == NULL) {
        perror("km_create: allocate context");
        return NULL;
    }
    memset(ctx, 0, sizeof(struct km_ctx_s));

    ctx->dimension = d;
    ctx->nclusters = k;
    ctx->nobserves = n;
    ctx->max_iters = m;
    ctx->inits_set = 0;

    ctx->cardinals = (size_t *) malloc (ctx->nclusters * sizeof(size_t));
    if (ctx->cardinals == NULL) {
        perror("km_create: allocate cardinalities vector");
        km_destroy(ctx);
        return NULL;
    }
    memset(ctx->cardinals, 0, sizeof(size_t) * ctx->nclusters);

    ctx->mean_vals = (double *) malloc (ctx->nclusters * ctx->dimension * sizeof(double));
    if (ctx->mean_vals == NULL) {
        perror("km_create: allocate centroids matrix");
        km_destroy(ctx);
        return NULL;
    }
    memset(ctx->mean_vals, 0, sizeof(double) * ctx->nclusters * ctx->dimension);

    ctx->data_vals = (double *) malloc (ctx->nobserves * ctx->dimension * sizeof(double));
    if (ctx->data_vals == NULL) {
        perror("km_create: allocate observation matrix");
        km_destroy(ctx);
        return NULL;
    }
    memset(ctx->data_vals, 0, sizeof(double) * ctx->nobserves * ctx->dimension);

    return ctx;
}

double km_distance_squared(struct km_ctx_s * ctx, size_t cluster_id, double * w)
{
    double d = 0, a = 0;
    size_t j;

    for (j = 0; j < ctx->dimension; ++j) {
        a = (*(ctx->mean_vals + (ctx->dimension * cluster_id) + j)) - (*(w + j));
        d += a * a;
    }

    return d;
}

size_t km_cluster(struct km_ctx_s * ctx, double * w)
{
    size_t i, s = 0;
    double last, curr = 0;
    last = km_distance_squared(ctx, 0, w);
    for (i = 1; i < ctx->nclusters; ++i) {
        curr = km_distance_squared(ctx, i, w);
        if (curr < last) {
            s = i;
            last = curr;
        }
    }
    return s;
}

void km_update(struct km_ctx_s * ctx, size_t i, double * w)
{
    size_t j;
    if (i < ctx->nclusters) {
        for (j = 0; j < ctx->dimension; ++j) {
            *(ctx->mean_vals + (ctx->dimension * i) + j) = w[j];
        }
    }

    for (j = 0; j < ctx->dimension; ++j) {
        *(ctx->data_vals + (ctx->dimension * i) + j) = w[j];
    }
}

int km_scan_input(struct km_ctx_s * ctx)
{
    char c = 0;
    double f = 0;
    size_t i, j;
    double * w = NULL;
    
    w = (double *) malloc (sizeof(double) * ctx->dimension);
    if (w == NULL) {
        perror("km_scan_input: malloc");
        return -1;
    }

    for (i = 0; i < ctx->nobserves; ++i) {
        for (j = 0; j < ctx->dimension; ++j) {
            if (scanf("%lf%c", &f, &c) != 2) {
                perror("km_scan_input: bad input");
                return -1;
            }
            *(w + j) = f;
        }
        km_update(ctx, i, w);
    }

    free(w);
    w = NULL;

    return 0;
}

int km_iterate(struct km_ctx_s * ctx)
{
    size_t i, j, s = 0;
    double * new_means = NULL, * temp = NULL;
    size_t * new_cards = NULL;

    new_means = (double *) malloc (ctx->nclusters * ctx->dimension * sizeof(double));
    if (new_means == NULL) {
        perror("km_iterate allocate new centroids matrix");
        return -1;
    }
    memset(new_means, 0, sizeof(double) * ctx->nclusters * ctx->dimension);

    new_cards = (size_t *) malloc (ctx->nclusters * sizeof(size_t));
    if (new_cards == NULL) {
        perror("km_iterate allocate new cardinalities vector");
        return -1;
    }
    memset(new_cards, 0, sizeof(size_t) * ctx->nclusters);

    /* iterate over all observavtions */
    for (i = 0; i < ctx->nobserves; ++i) {
        /* calculate cluster index for the current observation */
        s = km_cluster(ctx, ctx->data_vals + (i * ctx->dimension));

        /* iterate over the observation's entries */
        for (j = 0; j < ctx->dimension; ++j) {
            /* new_means[s][j] = (new_means[s][j] * new_cards[s] + observation[j]) / (new_cards[s] + 1) */
            *(new_means + (s * ctx->dimension) + j) *= new_cards[s];
            *(new_means + (s * ctx->dimension) + j) += *(ctx->data_vals + (i * ctx->dimension) + j);
            *(new_means + (s * ctx->dimension) + j) /= new_cards[s] + 1;
        }
        /* new_cards[s] += 1 */
        new_cards[s]++;
    }

    /* compare new_means to the current ones */
    for (i = 0; i < ctx->nclusters; ++i) {
        for (j = 0; j < ctx->dimension; ++j) {
            if ( *(new_means + (i * ctx->dimension) + j) != *(ctx->mean_vals + (i * ctx->dimension) + j) ) {
                temp = ctx->mean_vals;
                ctx->mean_vals = new_means;
                free(temp);
                free(new_cards);
                return 0;
            }
        }
    }
    free(new_means);
    free(new_cards);
    return 1;
}

int km_converge(struct km_ctx_s * ctx)
{
    int r = -1;
    size_t iter;

    for (iter = 0; iter < ctx->max_iters; ++iter) {
        if ((r = km_iterate(ctx)) == 1) {
            return 1;
        } else if (r < 0) {
            perror("km_converge: km_iterate failed");
            return  -1;
        }
    }

    return 0;
}

int main(int argc, char *argv[]) {
    size_t d = 0, k = 0, n = 0, m = 0;
    struct km_ctx_s * ctx = NULL;

    if (argc != 5) {
        printf("usage: %s K N d MAX_ITER < infile > outfile\n", argv[0]);
        return 1;
    }

    k = atoi(argv[1]);
    n = atoi(argv[2]);
    d = atoi(argv[3]);
    m = atoi(argv[4]);

    if (k <= 0 || n <= k || d <= 0 || m <= 0) {
        printf("%s: bad argument\n", argv[0]);
        return 2;
    }

    ctx = km_create(d, k, n, m);
    if (ctx == NULL) {
        perror("main: km_create failed");
        return -1;
    }

    if (km_scan_input(ctx) < 0) {
        perror("main: km_scan_input failed");
        return -2;
    }

    if (km_converge(ctx) < 0) {
        perror("main: km_converge failed");
        return -3;
    }

    km_dump(ctx);

    km_destroy(ctx);

    return 0;
}
