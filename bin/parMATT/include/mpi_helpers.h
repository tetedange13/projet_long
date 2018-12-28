/*
 *  Copyright (c) 2017 Maksim Shegay.
 *  The faculty of Computational Mathematics and Cybernetics of Lomonosov Moscow State University
 *  All rights reserved.
 *
 *  parMatt is licensed under the GNU public license version 2.0.
 */

#ifndef MPI_HELPERS_H
#define MPI_HELPERS_H

#define MASTER 0

#define MPI_TAG_PAIRSCORES 1
#define MPI_TAG_TASK 2
#define MPI_TAG_TASK_REQ 3
#define MPI_TAG_TASK_INT 4

#define MATT_LOG_LVL "MATT_LOG_LVL"

#include <pthread.h>

static int LOG_LVL;

/* mpi_init_log reads LOG_LVL value from MATT_LOG_LVL env variable
   if variable is undefined than default value will be used.

   LOL_LVL 0 -- log only timings for pairwise alingments and iterative part of the algorithm
   LOG_LVL 1 -- all above plus the information about tasks sent
   LOG_LVL 2 -- all above plus timings for each alignment and barrier wait time (default)
*/
void mpi_init_log();
void mpi_log(const int loglvl, const int master_only, const char *fmt, ...);

typedef struct mpi_task_t {
	/* alignments IDs al2 is ingored in iterative phase */
	int al1;
	int al2;
} mpi_task_t;


typedef mpi_task_t *(*generate_cb_t)(int n_als, int *n_tasks);

/* mpi_send_tasks sends task(s) to a reciever preocess with MPI_TAG_TASK
   for regular tasks and MPI_TAG_TASK_INT for termination tasks. 
   Returns mpi error code for corresponding MPI_Send method */
int mpi_send_tasks(mpi_task_t *t, const int len, const int reciever, const int tag);

/* mpi_recv_tasks recieves tasks from sender, sets len to -1 if
   fetches termination task.
   Returns mpi error code for corresponding MPI_Recv method */
int mpi_recv_tasks(mpi_task_t **t, int *len, const int sender);

mpi_task_t mpi_create_task(const int al1, const int al2);


typedef struct mpi_task_storage_t {
	pthread_mutex_t *lock;
	pthread_cond_t *stop;
	
	mpi_task_t *tasks;

	int size;
	int nprocs;
	int chunk_size;
	int terminated;
	mpi_task_t empty;
} mpi_task_storage_t;

mpi_task_storage_t *mpi_task_storage_create(generate_cb_t generate, const int n_als, const int nprocs, const int chunk_size);
void mpi_task_storage_destroy(mpi_task_storage_t *t);
int mpi_task_storage_size(mpi_task_storage_t *t);
mpi_task_t *mpi_task_storage_get_chunk(mpi_task_storage_t *t, int *len, int *tag);
void mpi_task_storage_wait_empty(mpi_task_storage_t *t);


#endif
