#include <stdlib.h>
#include <string.h>
#include <stdio.h>
struct km_ctx_s {
    size_t   dimension; // dimension
    size_t   nclusters; // # of clusters
    size_t   nobserves; // total # of observations
    size_t   max_iters; // MAX_ITER
    size_t   inits_set; // number initial mean values inserted so far
    size_t * cardinals; // cardinalities of clusters, array of size k
    double * mean_vals; // mean values, array of k * d
};

void km_dump(struct km_ctx_s * ctx)
{
    int i, j;
    for (i = 0; i < ctx->nclusters; ++i) {
        for (j = 0; j < ctx->dimension - 1; ++j) {
            printf("%.2f,", *(ctx->mean_vals + (ctx->dimension * i) + j));
        }
        printf("%.2f\n", *(ctx->mean_vals + (ctx->dimension * i) + j + 1));
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

    ctx->cardinals = (size_t *) calloc (sizeof(size_t), ctx->nclusters);
    if (ctx->cardinals == NULL) {
        perror("km_create: allocate cardinalities vector");
        free(ctx);
        return NULL;
    }
    memset(ctx->cardinals, 0, sizeof(size_t) * ctx->nclusters);

    ctx->mean_vals = (double *) calloc (sizeof(double), ctx->nclusters * ctx->dimension);
    if (ctx->mean_vals == NULL) {
        perror("km_create: allocate observation matrix");
        free(ctx->cardinals);
        free(ctx);
        return NULL;
    }
    memset(ctx->mean_vals, 0, sizeof(double) * ctx->nclusters * ctx->dimension);

    return ctx;
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
        free(ctx);
    }
}

double km_distance(struct km_ctx_s * ctx, size_t cluster_id, double * w)
{
    double d = 0, a = 0;
    int j;

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
    last = km_distance(ctx, 0, w);
    for (i = 1; i < ctx->nclusters; ++i) {
        curr = km_distance(ctx, i, w);
        if (curr < last) {
            s = i;
            last = curr;
        }
    }
    return s;
}

void km_update_at(struct km_ctx_s * ctx, size_t cluster_id, double * w)
{
    double * u = NULL;
    size_t j, c = 0;
    for (j = 0; j < ctx->dimension; ++j) {
        c = (*(ctx->cardinals + cluster_id))++;
        u = ctx->mean_vals + (ctx->dimension * cluster_id) + j;
        *u = ((*u)*c + (*(w + j))) / (c + 1);
    }
}

void km_update(struct km_ctx_s * ctx, double * w)
{
    size_t s = 0;
    if (ctx->inits_set == ctx->nclusters) {
        s = km_cluster(ctx, w);
    } else {
        s = ctx->inits_set++;
    }
    km_update_at(ctx, s, w);
}

int main(int argc, char *argv[]) {
    size_t d = 0, k = 0, n = 0, m = 0;
    struct km_ctx_s * ctx = NULL;
    d = atoi(argv[1]);
    k = atoi(argv[2]);
    n = atoi(argv[3]);
    m = atoi(argv[4]);
    ctx = km_create(d, k, n, m);
    if (ctx == NULL) {
        return 1;
    }
    km_dump(ctx);
    km_destroy(ctx);
    return 0;
}
