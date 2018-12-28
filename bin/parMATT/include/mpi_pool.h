/*
 *  Copyright (c) 2017 Maksim Shegay.
 *  The faculty of Computational Mathematics and Cybernetics of Lomonosov Moscow State University
 *  All rights reserved.
 *
 *  parMatt is licensed under the GNU public license version 2.0.
 */

#ifndef MPI_POOL_H
#define MPI_POOL_H

#include <pthread.h>

#include "mpi_helpers.h"

/* fetch_cb_t is a callback for taskpool_t to fetch next task, sets n to -1 if there is no more tasks 
   fetch frees fetched tasks  */
typedef mpi_task_t *(*fetch_cb_t)(int *n);
/* execute_cb_t is a callback for taskpool_t to execute task, appends return value to the associated processed_t vector */
typedef void *(*execute_cb_t)(mpi_task_t t);

struct node_t {
	mpi_task_t task;
	struct node_t *next;
};

/* taskpool_t is a pool of workers running on nworkers threads plus one thread for fetching the tasks */
typedef struct taskpool_t {
	int len;
	int nworkers;

	struct node_t *head;
	struct node_t *tail;

	pthread_mutex_t *lock;
	pthread_cond_t *not_empty;
	pthread_cond_t *more_tasks;

	fetch_cb_t fetch;
	execute_cb_t execute;
} taskpool_t;

/* procssed_t is a vector with a thread-safe append operation */
typedef struct processed_t {
	int len;
	int cap;

	pthread_mutex_t *lock;
	void **data;
} processed_t;

typedef struct taskpool_ctx_t {
	processed_t *processed;
	taskpool_t *pool;
} taskpool_ctx_t;


taskpool_t *taskpool_init(int n, fetch_cb_t fetch, execute_cb_t execute);
void taskpool_destroy(taskpool_t *p);
void taskpool_run(int nworkers, taskpool_ctx_t ctx, fetch_cb_t fetch, execute_cb_t execute);

processed_t *processed_init();
void processed_destroy(processed_t *p);
void processed_append(processed_t *p, void *x);

#endif
