#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#define MAX_THREADS 20

int total_hits, hits[MAX_THREADS], num_threads, sample_points, sample_points_per_thread;

void *compute_pi(void *s)
{
    int i, hits;
    double x,y,d;
    int *hit_pointer = (int*) s;
    unsigned int seed = *hit_pointer;
    hits = 0;

    for(i = 0; i < sample_points_per_thread; i++)
    {
        x = (double) ((rand_r(&seed)) / (double) (RAND_MAX));
        y = (double) ((rand_r(&seed)) / (double) (RAND_MAX));

        d = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5);
        if(d < 0.25) hits ++;
        seed *= (i + 1);
    }
    *hit_pointer = hits;
    printf("%d\n", hits);
    pthread_exit(NULL);

}

int main(int argc, char *argv[])
{
    pthread_t p_threads[MAX_THREADS];
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    int i;

    total_hits = 0;

    sample_points_per_thread = sample_points/num_threads;

    for(i = 0; i <num_threads; i++)
    {
        hits[i] = i;
        
        //int pthread_create(pthread_t * thread, const pthread_attr_t * attr, void * (*start_routine)(void *), void *arg);
        pthread_create(&p_threads[i], &attr, compute_pi, (void *) &hits[i]);
        printf("Created thread %d\n", i);

    }
    pthread_attr_destroy(&attr);

    for(i = 0; i <num_threads;i++)
    {
        pthread_join(p_threads[i], NULL);
        total_hits += hits[i];
    }
    printf("%d", total_hits);

    pthread_exit(NULL);
    return 0;
}


