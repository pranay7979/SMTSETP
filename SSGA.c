/* Approach 1 - set 1 */

/*

Initialize random population
Evaluate population
Store best solution

while termination not met:
    Select Parent 1 (binary tournament)
    Select Parent 2 (binary tournament)

    Offspring = crossover(Parent1, Parent2)
    if Offspring not unique:
        continue

    Mutant = mutation(Offspring)
    if Mutant not unique:
        continue

    BestNew = min(Offspring, Mutant)

    Replace worst individual with BestNew

    Apply API local search on BestNew

    Update global best solution

*/

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <windows.h>

#define MAXN 500 /*maximum no. of jobs*/
#define MAX_T 9  /*maximum no. of tardiness due dates*/
#define MAX_E 2  /*maximum no. of earliness due dates */
#define STAGNATION_LIMIT 2000

#define POPSIZE 100
#define MAX_ITERS 10000
#define THETA_RATIO 0.3 /*crpssover*/
#define INDEPENDENT_RUNS 10
#define STAGNATION_LIMIT 2000

#define BASE_PATH "C:\\Thesis Project\\SET_I"

/* ================= DATA STRUCTURES ================= */

typedef struct
{
    int n, mT, mE;
    int p[MAXN];             /*processing time*/
    int e_dead[MAX_E][MAXN]; /* one job- multiple deadlines so we use 2d array*/
    int e_cost[MAX_E][MAXN];
    int d[MAX_T][MAXN];
    int w[MAX_T][MAXN];
} Instance;

typedef struct
{
    int perm[MAXN]; /*one job schedule*/
    long long cost;
} Individual;

/* ================= GLOBALS ================= */

Instance I;
Individual pop[POPSIZE];

/* ================= UTILITIES ================= */

int rand_int(int a, int b)
{
    return a + rand() % (b - a + 1);
}

void swap(int *a, int *b)
{
    int t = *a;
    *a = *b;
    *b = t;
}

void shuffle(int *a, int n)
{
    for (int i = n - 1; i > 0; i--)
    {
        int j = rand_int(0, i);
        swap(&a[i], &a[j]);
    }
}

/* ================= OBJECTIVE FUNCTION ================= */
/* calculating the total stepwise cost incurred for a schedule */

long long evaluate(const int *perm)
{
    long long time = 0, cost = 0;

    for (int pos = 0; pos < I.n; pos++)
    {
        int j = perm[pos];
        time += I.p[j];
        int C = (int)time;

        int early = 0;
        for (int e = 0; e < I.mE; e++)
        {
            if (C <= I.e_dead[e][j])
            {
                cost += I.e_cost[e][j];
                early = 1;
                break;
            }
        }
        if (early)
            continue; /*if job is early so it will not be late so, we skip tardy part*/

        if (I.mT > 0 && C > I.d[0][j])
        {
            int applied = 0;
            for (int k = 1; k < I.mT; k++)
            {
                if (C <= I.d[k][j])
                {
                    cost += I.w[k - 1][j];
                    applied = 1;
                    break;
                }
            }
            if (!applied)
                cost += I.w[I.mT - 1][j];
        }
    }
    return cost;
}

/* ================= UNIQUENESS CHECK ================= */

int is_duplicate(const int *perm)
{
    for (int i = 0; i < POPSIZE / 2; i++)
    {
        if (memcmp(pop[i].perm, perm, sizeof(int) * I.n) == 0)
            return 1; /*duplicate*/
    }
    return 0; /* unique*/
}

/* ================= BINARY TOURNAMENT SELECTION ================= */

Individual tournament_selection()
{
    /*selecting two parents for crossover*/
    int a = rand_int(0, POPSIZE - 1);
    int b = rand_int(0, POPSIZE - 1);
    while (b == a)
        b = rand_int(0, POPSIZE - 1);
    return (pop[a].cost < pop[b].cost) ? pop[a] : pop[b];
}

/* ================= θ-POSITION UOB CROSSOVER ================= */
/* selecting random positions from parent 1 and remaining from parent 2 */

Individual crossover(const Individual *p1, const Individual *p2)
{
    Individual child;
    int used[MAXN] = {0};

    for (int i = 0; i < I.n; i++)
        child.perm[i] = -1;

    int theta = (int)(THETA_RATIO * I.n);
    if (theta < 1)
        theta = 1;

    int selected = 0;
    while (selected < theta)
    {
        int pos = rand_int(0, I.n - 1);
        if (child.perm[pos] == -1)
        {
            int job = p1->perm[pos];
            child.perm[pos] = job;
            used[job] = 1;
            selected++;
        }
    }

    int idx = 0;
    for (int i = 0; i < I.n; i++)
    {
        int job = p2->perm[i];
        if (!used[job])
        {
            while (child.perm[idx] != -1)
                idx++;
            child.perm[idx] = job;
        }
    }

    child.cost = evaluate(child.perm);
    return child;
}

/* ================= MUTATION: DOUBLE SWAP ================= */

void mutation(Individual *ind)
{
    int x = rand_int(0, I.n - 1);
    int y = rand_int(0, I.n - 1);
    int z = rand_int(0, I.n - 1);

    while (y == x)
        y = rand_int(0, I.n - 1);
    while (z == x || z == y)
        z = rand_int(0, I.n - 1);

    swap(&ind->perm[x], &ind->perm[y]);
    swap(&ind->perm[x], &ind->perm[z]);

    ind->cost = evaluate(ind->perm);
}

/* ================= API LOCAL SEARCH ================= */

/* Adajacent Pairwise Interchange - swap adjacent position and check cost, if improved keep it, if cost same check for completion time. */
void api_local_search(Individual *ind)
{
    int improved = 1;

    while (improved)
    {
        improved = 0;

        for (int i = 0; i < I.n - 1; i++)
        {
            int temp[MAXN];
            memcpy(temp, ind->perm, sizeof(int) * I.n);

            swap(&temp[i], &temp[i + 1]);

            long long new_cost = evaluate(temp);

            if (new_cost < ind->cost)
            {
                memcpy(ind->perm, temp, sizeof(int) * I.n);
                ind->cost = new_cost;
                improved = 1;
            }
        }
    }
}

/* ================= POPULATION SORTING ================= */
/* comparator for sorting population in ascending order of cost */

int compare_individuals(const void *a, const void *b)
{
    const Individual *x = (const Individual *)a;
    const Individual *y = (const Individual *)b;

    if (x->cost < y->cost)
        return -1;
    if (x->cost > y->cost)
        return 1;
    return 0;
}

/* ================= REPLACEMENT ================= */
/* population maintained in sorted order and replacing with worst one */

void replace_worst(const Individual *child)
{
    if (child->cost >= pop[POPSIZE - 1].cost)
        return;

    pop[POPSIZE - 1] = *child;
    qsort(pop, POPSIZE, sizeof(Individual), compare_individuals);
}

/* ================= INSTANCE READER ================= */

int read_instance(const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f)
        return 0;

    fscanf(f, "%d %d %d", &I.n, &I.mT, &I.mE);

    for (int j = 0; j < I.n; j++)
    {
        fscanf(f, "%d", &I.p[j]);

        for (int e = 0; e < I.mE; e++)
            fscanf(f, "%d", &I.e_dead[e][j]);

        for (int k = 0; k < I.mT; k++)
            fscanf(f, "%d", &I.d[k][j]);

        for (int e = 0; e < I.mE; e++)
            fscanf(f, "%d", &I.e_cost[e][j]);

        for (int k = 0; k < I.mT; k++)
            fscanf(f, "%d", &I.w[k][j]);
    }

    fclose(f);
    return 1;
}

/* ================= RECURSIVE DIRECTORY PROCESSING ================= */

void process_directory(const char *base, FILE *csv)
{
    WIN32_FIND_DATAA fd; /* a Windows-specific data structure used for directory traversal. */
    char search[MAX_PATH];
    sprintf(search, "%s\\*", base);

    HANDLE h = FindFirstFileA(search, &fd);
    if (h == INVALID_HANDLE_VALUE)
        return;

    do
    {
        if (!strcmp(fd.cFileName, ".") || !strcmp(fd.cFileName, ".."))
            continue;

        char path[MAX_PATH];
        sprintf(path, "%s\\%s", base, fd.cFileName);

        if (fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
        {
            process_directory(path, csv);
        }
        else
        {

            /* Print which instance file is being processed */
            printf("Processing Instance: %s\n", path);

            if (!read_instance(path))
            {
                printf("Failed to read instance: %s\n", path);
                continue;
            }

            for (int run = 1; run <= INDEPENDENT_RUNS; run++)
            {

                printf("   Run %d / %d started...\n", run, INDEPENDENT_RUNS);

                srand((unsigned)(time(NULL) ^ (run * 7919))); /* All future calls to rand() depend on this seed , different for each independent runs*/

                clock_t start = clock();

                /* Initialize random population */
                for (int i = 0; i < POPSIZE; i++)
                {
                    for (int j = 0; j < I.n; j++)
                        pop[i].perm[j] = j;    // [1,2,3,4,5]
                    shuffle(pop[i].perm, I.n); // [2,5,3,1,4]
                    pop[i].cost = evaluate(pop[i].perm);
                }

                qsort(pop, POPSIZE, sizeof(Individual), compare_individuals); /* sort population according to cost */

                long long init_best = pop[0].cost;
                long long best = init_best;

                /* Evolution loop */
                /* Evolution loop with stagnation-based termination */

                best = pop[0].cost;
                int no_improve = 0;

                for (int it = 0; it < MAX_ITERS && no_improve < STAGNATION_LIMIT; it++)
                {

                    Individual p1 = tournament_selection();
                    Individual p2 = tournament_selection();

                    Individual c = crossover(&p1, &p2);
                    mutation(&c);

                    if (rand() % 100 < 30)
                        api_local_search(&c);

                    if (is_duplicate(c.perm))
                        continue;

                    replace_worst(&c);

                    /* Check for global best improvement */
                    if (pop[0].cost < best)
                    {
                        best = pop[0].cost;
                        no_improve = 0; /* reset stagnation counter */
                    }
                    else
                    {
                        no_improve++; /* increase stagnation counter */
                    }
                }

                double cpu =
                    (double)(clock() - start) / CLOCKS_PER_SEC;

                printf("   Run %d completed | Final Best = %lld | Time = %.3f sec\n",
                       run, best, cpu);

                fprintf(csv, "%s,%d,%d,%lld,%lld,%.4f\n",
                        path, run, I.n, init_best, best, cpu);
            }

            printf("Finished Instance: %s\n", path);
        }

    } while (FindNextFileA(h, &fd));

    FindClose(h);
}

/* ================= MAIN ================= */

int main()
{
    FILE *csv = fopen("corrected_code_results.csv", "w");
    fprintf(csv,
            "InstancePath,RunID,n,InitialBestCost,FinalBestCost,CPUTimeSeconds\n");
    process_directory(BASE_PATH, csv);
    fclose(csv);
    printf("\nAll experiments completed successfully.\n");
    return 0;
}
