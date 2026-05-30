/*approach 2_set2 - Initial population generation using heuristic*/

/*

moore-
M-Moore heuristic incrementally builds a schedule by inserting jobs in increasing due-date order and removes the longest job whenever a deadline violation occurs, repeating the process for multiple due dates.

NEH-
Generate a random permutation of jobs ω.

Select the first two jobs from ω and evaluate the objective values of the two alternative permutations formed by these two jobs

Select the permutation which is having minimum objective value.

Select the kth job Jωk and insert it on k available places in the incomplete solution.
The process results in k incomplete permutations of jobs.

Choose the permutation which incurs the least objective value and repeat this process until a complete solution is generated.

*/

/*
architecture-

Initialize population:
    1 Moore heuristic solution
    1 earliness moore
    3 NEH heuristic solutions
    95 Random solutions

Evaluate population
Sort population
Store best solution

while termination not met:

    Select Parent 1 (binary tournament)
    Select Parent 2 (binary tournament)

    Offspring = crossover(Parent1, Parent2)
    if not unique: continue

    Mutant = mutation(Offspring)
    if not unique: continue

    BestNew = min(Offspring, Mutant)

    Replace worst individual

    Apply API local search

    Update global best

*/
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <windows.h>

#define MAXN 500
#define MAX_T 9
#define MAX_E 2
#define POPSIZE 100
#define MAX_ITERS 10000
#define THETA_RATIO 0.3
#define INDEPENDENT_RUNS 10
#define BASE_PATH "C:\\Thesis Project\\SET_II"

typedef struct
{
    int n, mT, mE;
    int p[MAXN];             // processing time of jobs
    int e_dead[MAX_E][MAXN]; // earliness due dates
    int e_cost[MAX_E][MAXN]; // earliness penalties
    int d[MAX_T][MAXN];      // tardiness due dates
    int w[MAX_T][MAXN];      // tardiness penalties
} Instance;

typedef struct
{
    int perm[MAXN]; // one schedule of jobs
    long long cost;
} Individual;

Instance I;
Individual pop[POPSIZE];

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

/*fitness evalaution function*/

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
            continue; /*if the job is early then it will be not late*/

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

/* finding the unique solution or not */

int is_duplicate(const int *perm, long long current_cost)
{
    for (int i = 0; i < POPSIZE; i++)
    {
        // Only compare sequences if the costs are identical (performance boost)
        if (pop[i].cost == current_cost)
        {
            if (memcmp(pop[i].perm, perm, sizeof(int) * I.n) == 0)
                return 1; // True duplicate
        }
    }
    return 0; 
}

/* binary tournament selection function*/

Individual tournament_selection()
{
    int a = rand_int(0, POPSIZE - 1);
    int b = rand_int(0, POPSIZE - 1);
    while (b == a)
        b = rand_int(0, POPSIZE - 1);

    return (pop[a].cost < pop[b].cost) ? pop[a] : pop[b];
}

/* crossover theta positions */

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

/* 2 swap mutation */

Individual mutate(const Individual *parent)
{
    Individual child = *parent;

    int x = rand_int(0, I.n - 1);
    int y = rand_int(0, I.n - 1);
    int z = rand_int(0, I.n - 1);

    while (y == x)
        y = rand_int(0, I.n - 1);
    while (z == x || z == y)
        z = rand_int(0, I.n - 1);

    swap(&child.perm[x], &child.perm[y]);
    swap(&child.perm[x], &child.perm[z]);

    child.cost = evaluate(child.perm);
    return child;
}

/* local search - api */

void api_local_search(Individual *ind)
{
    int max_passes = 2;   // limit depth
    int pass = 0;

    while (pass < max_passes)
    {
        int improved = 0;

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

        if (!improved)
            break;

        pass++;
    }
}


/* comparator for our sorting function */
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

/* replacing the worst solution in population with current solution*/

void replace_worst(const Individual *child)
{
    if (child->cost >= pop[POPSIZE - 1].cost)
        return;

    pop[POPSIZE - 1] = *child;
    qsort(pop, POPSIZE, sizeof(Individual), compare_individuals);
}

/* neh heuristic*/
long long evaluate_partial(const int *perm, int length)
{
    long long time = 0, cost = 0;

    for (int pos = 0; pos < length; pos++)
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
            continue;

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
Individual neh_heuristic_optimized()
{
    Individual result;
    int hold[MAXN];
    int omega[MAXN];
    
    // Step 1: Create the initial priority order
    for (int i = 0; i < I.n; i++) omega[i] = i;
    shuffle(omega, I.n);

    // Step 2: Start with the first job
    hold[0] = omega[0];
    int current_size = 1;

    // Step 3: Iteratively insert the remaining jobs
    for (int k = 1; k < I.n; k++)
    {
        int job_to_insert = omega[k];
        long long best_cost = LLONG_MAX;
        int best_pos = 0;

        // Try inserting 'job_to_insert' at every possible position (0 to current_size)
        for (int pos = 0; pos <= current_size; pos++)
        {
            int trial[MAXN];
            
            // Construct the trial sequence
            for (int i = 0; i < pos; i++) trial[i] = hold[i];
            trial[pos] = job_to_insert;
            for (int i = pos; i < current_size; i++) trial[i + 1] = hold[i];

            // Evaluate the partial schedule cost
            long long c = evaluate_partial(trial, current_size + 1);

            if (c < best_cost)
            {
                best_cost = c;
                best_pos = pos;
            }
        }

        // Permanently insert at the best position found
        for (int i = current_size; i > best_pos; i--) hold[i] = hold[i - 1];
        hold[best_pos] = job_to_insert;
        current_size++;
    }

    memcpy(result.perm, hold, sizeof(int) * I.n);
    result.cost = evaluate(result.perm);
    return result;
}

Individual neh_heuristic()
{
    Individual result;
    int omega[MAXN];
    for (int i = 0; i < I.n; i++)
        omega[i] = i;

    shuffle(omega, I.n);

    int temp_perm[MAXN];
    temp_perm[0] = omega[0];
    temp_perm[1] = omega[1];
    int alt[2] = {omega[1], omega[0]};
    if (evaluate_partial(alt, 2) < evaluate_partial(temp_perm, 2))
    {
        temp_perm[0] = omega[1];
        temp_perm[1] = omega[0];
    }

    int current_size = 2;

    for (int k = 2; k < I.n; k++)
    {

        int job = omega[k];
        long long min_cost = LLONG_MAX;
        int best_pos = 0;

        for (int pos = 0; pos <= current_size; pos++)
        {

            int trial[MAXN];

            for (int i = 0, j = 0; i <= current_size; i++)
            {
                if (i == pos)
                    trial[i] = job;
                else
                    trial[i] = temp_perm[j++];
            }

            long long c = evaluate_partial(trial, current_size + 1);
            if (c < min_cost)
            {
                min_cost = c;
                best_pos = pos;
            }
        }

        for (int i = current_size; i > best_pos; i--)
            temp_perm[i] = temp_perm[i - 1];

        temp_perm[best_pos] = job;
        current_size++;
    }

    memcpy(result.perm, temp_perm, sizeof(int) * I.n);
    result.cost = evaluate(result.perm);
    return result;
}

/*moore heuristic*/

int compare_due(const void *a, const void *b)
{
    int j1 = *(int *)a;
    int j2 = *(int *)b;
    return I.d[0][j1] - I.d[0][j2];
}

Individual moore_heuristic()
{
    Individual result;
    int jobs[MAXN];

    for (int i = 0; i < I.n; i++)
        jobs[i] = i;

    qsort(jobs, I.n, sizeof(int), compare_due);

    int sigma[MAXN], tau[MAXN];
    int s_size = 0, t_size = 0;
    long long time = 0;

    for (int i = 0; i < I.n; i++)
    {

        sigma[s_size++] = jobs[i];
        time += I.p[jobs[i]];

        if (time > I.d[0][jobs[i]])
        {

            int max_p = -1, idx = -1;
            for (int j = 0; j < s_size; j++)
            {
                if (I.p[sigma[j]] > max_p)
                {
                    max_p = I.p[sigma[j]];
                    idx = j;
                }
            }

            tau[t_size++] = sigma[idx];
            time -= I.p[sigma[idx]];

            for (int j = idx; j < s_size - 1; j++)
                sigma[j] = sigma[j + 1];

            s_size--;
        }
    }

    for (int i = 0; i < s_size; i++)
        result.perm[i] = sigma[i];

    for (int i = 0; i < t_size; i++)
        result.perm[s_size + i] = tau[i];

    result.cost = evaluate(result.perm);
    return result;
}
/* Reverse Moore heuristic - focusing on earliness due dates */

int compare_earliness(const void *a, const void *b)
{
    int j1 = *(int *)a;
    int j2 = *(int *)b;
    if (I.mT < 2) return I.d[0][j1] - I.d[0][j2];
    return I.d[1][j1] - I.d[1][j2];
}

Individual reverse_moore_heuristic()
{
    if (I.mT < 2)
    return moore_heuristic();

    Individual result;
    int jobs[MAXN];

    /* Initialize job list */
    for (int i = 0; i < I.n; i++)
        jobs[i] = i;

    /* Sort jobs by earliness due dates */
    qsort(jobs, I.n, sizeof(int), compare_earliness);

    int sigma[MAXN], tau[MAXN];
    int s_size = 0, t_size = 0;
    long long time = 0;

    /* Traverse from largest earliness due date to smallest */
    for (int i = I.n - 1; i >= 0; i--)
    {

        sigma[s_size++] = jobs[i];
        time += I.p[jobs[i]];

        /* If accumulated time exceeds earliness due date */
        if (time > I.d[1][jobs[i]])
        {

            /* Remove job with largest processing time */
            int max_p = -1, idx = -1;

            for (int j = 0; j < s_size; j++)
            {
                if (I.p[sigma[j]] > max_p)
                {
                    max_p = I.p[sigma[j]];
                    idx = j;
                }
            }
            /* Move removed job to tau set */
            tau[t_size++] = sigma[idx];
            time -= I.p[sigma[idx]];

            /* Shift remaining jobs */
            for (int j = idx; j < s_size - 1; j++)
                sigma[j] = sigma[j + 1];

            s_size--;
        }
    }
    /* Construct final permutation: sigma followed by tau */
    for (int i = 0; i < s_size; i++)
        result.perm[i] = sigma[i];

    for (int i = 0; i < t_size; i++)
        result.perm[s_size + i] = tau[i];

    /* Evaluate using your stepwise ET objective */
    result.cost = evaluate(result.perm);

    return result;
}

/* reading the instance file and adding values to our structure*/

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

/* ================= DIRECTORY PROCESSING ================= */

void process_directory(const char *base, FILE *csv)
{
    WIN32_FIND_DATAA fd;
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
            printf("Processing Instance: %s\n", path);
            if (!read_instance(path))
                continue;

            for (int run = 1; run <= INDEPENDENT_RUNS; run++)
            {

                printf("   Run %d started...\n", run);
                srand((unsigned)(time(NULL) ^ (run * 1234567)));

                clock_t start = clock();

                pop[0] = moore_heuristic();
                pop[1] = reverse_moore_heuristic();

                for (int i = 2; i < 5; i++)
                    pop[i] = neh_heuristic_optimized();

                for (int i = 5; i < POPSIZE; i++)
                {
                    for (int j = 0; j < I.n; j++)
                        pop[i].perm[j] = j;
                    shuffle(pop[i].perm, I.n);
                    pop[i].cost = evaluate(pop[i].perm);
                }

                qsort(pop, POPSIZE, sizeof(Individual), compare_individuals);

                long long best = pop[0].cost;
                long long init_best = pop[0].cost;
                int no_improve_count = 0;

                for (int it = 0; it < MAX_ITERS; it++)
                {

                    Individual p1 = tournament_selection();
                    Individual p2 = tournament_selection();

                    Individual offspring = crossover(&p1, &p2);
                    if (is_duplicate(offspring.perm, offspring.cost))
                        continue;

                    Individual mutant = mutate(&offspring);
                    if (is_duplicate(mutant.perm, mutant.cost))
                        continue;

                    Individual bestNew =
                        (offspring.cost < mutant.cost) ? offspring : mutant;

                    Individual improved = bestNew;
                    if (rand() % 100 < 30)  api_local_search(&improved);
                    replace_worst(&improved);

                    /* Update best solution */
                    if (pop[0].cost < best)
                    {
                        best = pop[0].cost;
                        no_improve_count = 0; // reset stagnation counter
                    }
                    else
                    {
                        no_improve_count++;
                    }

                    /* --- STAGNATION TERMINATION CONDITIONS --- */

                    /* Condition 1: No improvement for 2000 generations */
                    if (no_improve_count >= 2000)
                    {
                        printf("   Terminated due to 2000 stagnation at iteration %d\n", it);
                        break;
                    }

                    /* Condition 2: No improvement for 1000 generations in first 1500 */
                    if (it <= 1500 && no_improve_count >= 1000)
                    {
                        printf("   Early stagnation (1000 in first 1500) at iteration %d\n", it);
                        break;
                    }
                }

                double cpu = (double)(clock() - start) / CLOCKS_PER_SEC;

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

int main()
{
    FILE *csv = fopen("approach2_set2_results.csv", "w");
    fprintf(csv, "InstancePath,RunID,n,InitialBestCost,FinalBestCost,CPUTimeSeconds\n");

    process_directory(BASE_PATH, csv);

    fclose(csv);

    printf("\nAll experiments completed successfully.\n");
    return 0;
}
