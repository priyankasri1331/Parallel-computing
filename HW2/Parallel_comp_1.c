#include <pthread.h>
#include <stdio.h>
#define MAX_THREADS 8

void *start_routine(void *id)
{
    printf("Howdy from thread %d\n", id);
    pthread_exit(NULL);
}

int main(int argc, char *argv[])
{

    pthread_t threads[MAX_THREADS];
    int status, i;
    for(i = 0; i < MAX_THREADS;i++)
    {
        status = pthread_create(&threads[i],NULL, start_routine, (void *)i);

        if (status)
        {
            printf("Error craeting thread \n");
        }
    }
    pthread_exit(NULL); //main is also a thread. Main creates 8 threads.


    return 0;
}
