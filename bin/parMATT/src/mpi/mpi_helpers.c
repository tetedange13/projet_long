/*
 *  Copyright (c) 2017 Maksim Shegay.
 *  The faculty of Computational Mathematics and Cybernetics of Lomonosov Moscow State University
 *  All rights reserved.
 *
 *  parMatt is licensed under the GNU public license version 2.0.
 */

#include <assert.h>
#include <errno.h>
#include <mpi.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "MultipleAlignment.h"

#include "mpi_helpers.h"

static void mpi_init_task_type(MPI_Datatype *datatype)
{
	MPI_Type_contiguous(2, MPI_INT, datatype);
	MPI_Type_commit(datatype);
}

static void mpi_free_task_type(MPI_Datatype *datatype)
{
	MPI_Type_free(datatype);
}

int mpi_send_tasks(mpi_task_t *t, const int len, const int reciever, const int tag)
{
	MPI_Datatype datatype;
	mpi_init_task_type(&datatype);
	int ret = MPI_Send(t, len, datatype, reciever, tag, MPI_COMM_WORLD);
	mpi_free_task_type(&datatype);
	return ret;
}

int mpi_recv_tasks(mpi_task_t **t, int *len, const int sender)
{
	MPI_Datatype datatype;
	mpi_init_task_type(&datatype);

	MPI_Status status;
	MPI_Probe(sender, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

	MPI_Get_count(&status, datatype, len);
	*t = calloc(*len, sizeof(mpi_task_t));
	int ret = MPI_Recv(*t, *len, datatype, sender, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	if (status.MPI_TAG == MPI_TAG_TASK_INT) {
		*len = -1;
	}

	mpi_free_task_type(&datatype);
	return ret;
}


mpi_task_t mpi_create_task(const int al1, const int al2)
{
	mpi_task_t t = {
	    .al1 = al1,
	    .al2 = al2,
	};
	return t;
}

void mpi_init_log()
{
	char *v = getenv(MATT_LOG_LVL);
	if (v) {
		errno = 0;
		LOG_LVL = strtol(v, NULL, 10);
		if (errno == ERANGE || errno == EINVAL || LOG_LVL < 0 || LOG_LVL > 2) {
			LOG_LVL = 2;
		}
	} else {
		LOG_LVL = 2;
	}

	mpi_log(0, 0, "loglevel set to %d", LOG_LVL);
}


void mpi_log(const int loglvl, const int master_only, const char *fmt, ...)
{
	if (loglvl > LOG_LVL) {
		return;
	}

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (master_only && rank != MASTER) {
		return;
	}

	char buff[100];
	time_t now = time(0);
	strftime(buff, 100, "%Y-%m-%d %H:%M:%S", localtime(&now));
	if (rank == MASTER) {
		printf("%s master: ", buff);
	} else {
		printf("%s slave %d: ", buff, rank);
	}

	va_list arg;
	va_start(arg, fmt);
	vprintf(fmt, arg);
	va_end(arg);
	printf("\n");
}


mpi_task_storage_t *mpi_task_storage_create(generate_cb_t generate, const int n_als, const int nprocs, const int chunk_size)
{
	mpi_task_storage_t *t = calloc(1, sizeof(mpi_task_storage_t));
	assert(t);

	t->lock = malloc(sizeof(pthread_mutex_t));
	assert(t->lock);
	pthread_mutex_init(t->lock, NULL);

	t->stop = malloc(sizeof(pthread_cond_t));
	assert(t->stop);
	pthread_cond_init(t->stop, NULL);

	t->tasks = generate(n_als, &(t->size));

	t->nprocs = nprocs;
	t->chunk_size = chunk_size;
	return t;
}

void mpi_task_storage_destroy(mpi_task_storage_t *t)
{
	if (t == NULL) {
		return;
	}
	
	pthread_mutex_destroy(t->lock);
	free(t->lock);

	pthread_cond_destroy(t->stop);
	free(t->stop);

	free(t->tasks);
	free(t);
}


int mpi_task_storage_size(mpi_task_storage_t *t)
{
	pthread_mutex_lock(t->lock);
	int size = t->size;
	pthread_mutex_unlock(t->lock);
	return size;
}

mpi_task_t *mpi_task_storage_get_chunk(mpi_task_storage_t *t, int *len, int *tag)
{
	pthread_mutex_lock(t->lock);
	if (t->size == 0) {
		*len = 1;
		*tag = MPI_TAG_TASK_INT;
		t->terminated++;
		if (t->terminated == t->nprocs) {
			pthread_cond_signal(t->stop);
		}
		pthread_mutex_unlock(t->lock);
		return &(t->empty);
	}

	int want = t->size / t->nprocs + (t->size % t->nprocs > 0);
	want = want < t->chunk_size ? want : t->chunk_size;

	if (want > t->size) {
		want = t->size;
	}
	*len = want;
	t->size -= want;

	mpi_task_t *ret = t->tasks + t->size;
	*tag = MPI_TAG_TASK;
	pthread_mutex_unlock(t->lock);
	return ret;
}

void mpi_task_storage_wait_empty(mpi_task_storage_t *t)
{
	pthread_mutex_lock(t->lock);
	while (t->terminated < t->nprocs) {
		pthread_cond_wait(t->stop, t->lock);
	}
	pthread_mutex_unlock(t->lock);
}
