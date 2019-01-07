/*
 *  Copyright (c) 2017 Maksim Shegay.
 *  The faculty of Computational Mathematics and Cybernetics of Lomonosov Moscow State University
 *  All rights reserved.
 *
 *  parMatt is licensed under the GNU public license version 2.0.
 */

#include <assert.h>
#include <math.h>
#include <mpi.h>

#include "MultipleAlignment.h"

#include "mpi_helpers.h"
#include "mpi_pairscores.h"
#include "mpi_pool.h"
#include "mpi_serializer.h"

static MultipleAlignment **mpi_init_als(PDBChain **chains, const int n)
{
	MultipleAlignment **alignments = (MultipleAlignment **)malloc(n * sizeof(MultipleAlignment *));
	assert(alignments);
	for (int i = 0; i < n; i++) {
		alignments[i] = CreateSingleStrandAlignment(Distill(chains[i], i));
	}
	return alignments;
}

static mpi_task_t *mpi_generator_pairwise(int n_als, int *n_tasks)
{
	*n_tasks = (n_als * (n_als - 1)) / 2;
	mpi_task_t *tasks = calloc(*n_tasks, sizeof(mpi_task_t));
	assert(tasks);
	int ind = 0;
	for (int i = 0; i < n_als; i++) {
		for (int j = 0; j < i; j++) {
			tasks[ind++] = mpi_create_task(i, j);
		}
	}
	return tasks;
}

static mpi_task_t *mpi_generator_iterative(int n_als, int *n_tasks)
{
	mpi_task_t *tasks = calloc(n_als, sizeof(mpi_task_t));
	assert(tasks);
	for (int i = 0; i < n_als; i++) {
		tasks[i] = mpi_create_task(i, -1);
	}
	*n_tasks = n_als;
	return tasks;
}

static int mpi_next_chunk_size(const int chunk_size, const int n_tasks, const int nprocs)
{
	int ret = n_tasks / nprocs + (n_tasks % nprocs > 0);
	ret = ret < chunk_size ? ret : chunk_size;
	return ret;
}

struct mpi_master_ctx_t {
	mpi_task_storage_t *task_storage;
	int iteration;
	int nprocs;
	int n_tasks;
};

static void *mpi_master(void *context)
{
	struct mpi_master_ctx_t *ctx = context;
	mpi_task_storage_t *task_storage = ctx->task_storage;
	int iteration = ctx->iteration;
	int nprocs = ctx->nprocs;
	int n_tasks = ctx->n_tasks;

	char dummy;
	int len, tag;
	MPI_Status status;
	// master process recieves tasks not via mpi methods,
	// that's why we want to send terminate tasks only to nprocs-1 nodes
	for (int terminated = 0; terminated < nprocs - 1;) {
		MPI_Recv(&dummy, 1, MPI_CHAR, MPI_ANY_SOURCE, MPI_TAG_TASK_REQ, MPI_COMM_WORLD, &status);
		mpi_task_t *tasks = mpi_task_storage_get_chunk(task_storage, &len, &tag);
		mpi_send_tasks(tasks, len, status.MPI_SOURCE, tag);
		if (tag == MPI_TAG_TASK) {
			int tasks_left = mpi_task_storage_size(task_storage);
			if (iteration == -1) {
				mpi_log(1, 0, "pairwise, sent %d task(s) to %d; (%d/%d)", len, status.MPI_SOURCE, n_tasks - tasks_left, n_tasks);
			} else {
				mpi_log(1, 0, "iteration %d, sent %d task(s) to %d; (%d/%d)", iteration, len, status.MPI_SOURCE, n_tasks - tasks_left, n_tasks);
			}
			continue;
		}
		terminated++;
		if (iteration == -1) {
			mpi_log(1, 0, "pairwise, sent terminate task to %d", status.MPI_SOURCE);
		} else {
			mpi_log(1, 0, "iteration %d, sent terminate task to %d", iteration, status.MPI_SOURCE);
		}
	}
	mpi_task_storage_wait_empty(task_storage);
	return NULL;
}

static mpi_task_t *mpi_fetcher_slave(int *n)
{
	char dummy;
	MPI_Send(&dummy, 1, MPI_CHAR, MASTER, MPI_TAG_TASK_REQ, MPI_COMM_WORLD);
	mpi_task_t *tasks;
	mpi_recv_tasks(&tasks, n, MASTER);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	return tasks;
}

static void mpi_pairwise(MultipleAlignment **als, const int n_als, pairscores_t **scores, processed_t *aligned, int nthreads, int rank, int nprocs)
{
	mpi_task_storage_t *task_storage = NULL;
	pthread_t master;
	int n_tasks = 0;
	if (rank == MASTER) {
		task_storage = mpi_task_storage_create(mpi_generator_pairwise, n_als, nprocs, nthreads);
		n_tasks = mpi_task_storage_size(task_storage);
		struct mpi_master_ctx_t ctx = {
		    .task_storage = task_storage,
		    .iteration = -1,
		    .nprocs = nprocs,
		    .n_tasks = n_tasks,
		};
		pthread_create(&master, NULL, mpi_master, &ctx);
	}

	mpi_task_t *fetcher_master(int *n)
	{
		int tag;
		mpi_task_t *tasks = mpi_task_storage_get_chunk(task_storage, n, &tag);
		mpi_task_t *ret = malloc(*n * sizeof(mpi_task_t));
		assert(ret);
		memcpy(ret, tasks, *n * sizeof(mpi_task_t));
		// races are possible, but this variable is used only for logging
		int tasks_left = mpi_task_storage_size(task_storage);
		if (tag == MPI_TAG_TASK_INT) {
			*n = -1;
			mpi_log(1, 0, "pairwise, sent terminate task to 0");
			return ret;
		}
		mpi_log(1, 0, "pairwise, sent %d task(s) to 0; (%d/%d)", *n, n_tasks - tasks_left, n_tasks);
		return ret;
	}

	fetch_cb_t fetcher;
	if (rank == MASTER) {
		fetcher = fetcher_master;
	} else {
		fetcher = &mpi_fetcher_slave;
	}

	taskpool_ctx_t ctx = {
	    .processed = aligned,
	};

	pthread_mutex_t *lock = malloc(sizeof(pthread_mutex_t));
	assert(lock);
	pthread_mutex_init(lock, NULL);

	void *align(mpi_task_t t)
	{
		double t1 = MPI_Wtime();
		MultipleAlignment *al1 = als[t.al1];
		MultipleAlignment *al2 = als[t.al2];
		MultipleAlignment *ret = AlignAlignments(al1, al2, 0, 0, 0);
		pthread_mutex_lock(lock);
		pairscores_add(scores, t.al1, t.al2, ret->score);
		pthread_mutex_unlock(lock);
		double t2 = MPI_Wtime();
		mpi_log(2, 0, "pairwise, aligned %d %d in %1.5lf seconds", t.al1, t.al2, t2 - t1);
		return ret;
	}

	taskpool_run(nthreads, ctx, fetcher, align);

	pthread_mutex_destroy(lock);
	free(lock);

	if (rank != MASTER) {
		return;
	}

	pthread_join(master, NULL);
	mpi_task_storage_destroy(task_storage);
}

static void mpi_iterative(MultipleAlignment **als, int n_als, processed_t *aligned, MultipleAlignment *ma, pairscores_t **scores, int iteration, int nthreads, int rank, int nprocs)
{
	mpi_task_storage_t *task_storage = NULL;
	pthread_t master;
	int n_tasks;
	if (rank == MASTER) {
		task_storage = mpi_task_storage_create(mpi_generator_iterative, n_als, nprocs, nthreads);
		n_tasks = mpi_task_storage_size(task_storage);
		struct mpi_master_ctx_t ctx = {
		    .task_storage = task_storage,
		    .iteration = iteration,
		    .nprocs = nprocs,
		    .n_tasks = n_tasks,
		};
		pthread_create(&master, NULL, mpi_master, &ctx);
	}

	mpi_task_t *fetcher_master(int *n)
	{
		int tag;
		mpi_task_t *tasks = mpi_task_storage_get_chunk(task_storage, n, &tag);
		mpi_task_t *ret = malloc(*n * sizeof(mpi_task_t));
		assert(ret);
		memcpy(ret, tasks, *n * sizeof(mpi_task_t));
		// races are possible, but this variable is used only for logging
		int tasks_left = mpi_task_storage_size(task_storage);
		if (tag == MPI_TAG_TASK_INT) {
			*n = -1;
			mpi_log(1, 0, "iteration %d, sent terminate task to 0", iteration);
			return ret;
		}
		mpi_log(1, 0, "iteration %d, sent %d task(s) to 0; (%d/%d)", iteration, *n, n_tasks - tasks_left, n_tasks);
		return ret;
	}

	fetch_cb_t fetcher;
	if (rank == MASTER) {
		fetcher = fetcher_master;
	} else {
		fetcher = &mpi_fetcher_slave;
	}

	taskpool_ctx_t ctx = {
	    .processed = aligned,
	};

	void *align(mpi_task_t t)
	{
		double t1 = MPI_Wtime();
		MultipleAlignment *al = als[t.al1];
		double bestScore = 0;
		int best1 = 0, best2 = 0;
		for (int j = 0; j < al->numChains; j++) {
			for (int k = 0; k < ma->numChains; k++) {
				pairscores_t *s = pairscores_find(scores, ma->chains[k]->id, al->chains[j]->id);
				if (s->score > bestScore) {
					bestScore = s->score;
					best1 = j;
					best2 = k;
				}
			}
		}
		MultipleAlignment *ret = AlignAlignments(al, ma, best1, best2, 0);
		double t2 = MPI_Wtime();
		mpi_log(2, 0, "iteration %d, aligned %d in %1.5lf seconds", iteration, t.al1, t2 - t1);
		return ret;
	}

	taskpool_run(nthreads, ctx, fetcher, align);

	if (rank == MASTER) {
		pthread_join(master, NULL);
		mpi_task_storage_destroy(task_storage);
	}
}

static MultipleAlignment *mpi_bcast_best_ma(processed_t *p, int rank)
{
	struct best_t {
		double score;
		int rank;
	};

	struct best_t best = {
	    .score = -INFINITY,
	    .rank = rank,
	};
	int index = 0;
	for (int i = 0; i < p->len; i++) {
		MultipleAlignment *al = p->data[i];
		if (al->score > best.score || i == 0) {
			best.score = al->score;
			index = i;
		}
	}
	struct best_t out;
	MPI_Allreduce(&best, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
	char *buf = NULL;
	int size = 0;
	MultipleAlignment *ret = NULL;
	if (rank == out.rank) {
		buf = mpi_ma_serialize(p->data[index], &size);
		ret = p->data[index];
		p->data[index] = p->data[p->len - 1];
		p->len--;
	}
	MPI_Bcast(&size, 1, MPI_INT, out.rank, MPI_COMM_WORLD);
	if (rank != out.rank) {
		buf = calloc(size, 1);
		assert(buf);
	}
	MPI_Bcast(buf, size, MPI_CHAR, out.rank, MPI_COMM_WORLD);
	if (rank != out.rank) {
		ret = mpi_ma_deserialize(buf);
	}
	free(buf);

	return ret;
}

static int cleanup(MultipleAlignment **als, int n, const int index1, const int index2)
{
	int deleted = 0;
	for (int i = n - 1; i >= 0; i--) {
		for (int j = 0; j < als[i]->numChains; j++) {
			if (als[i]->chains[j]->id != index1 && als[i]->chains[j]->id != index2) {
				continue;
			}

			CleanupAlignment(als[i]);
			n--;
			als[i] = als[n];
			deleted++;
			break;
		}
	}

	return deleted;
}

static void mpi_cleanup_best(MultipleAlignment **als, int *n_als, processed_t *p, MultipleAlignment *ma)
{
	int index1 = ma->chains[0]->id;
	int index2 = ma->chains[ma->numChains - 1]->id;

	*n_als -= cleanup(als, *n_als, index1, index2);

	p->len -= cleanup((MultipleAlignment **)p->data, p->len, index1, index2);
}

static void mpi_barrier_wrapper(const int iteration)
{
	double t1 = MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD);
	double t2 = MPI_Wtime();
	if (iteration == -1) {
		mpi_log(2, 0, "iteration pairwise, waiting on barrier %1.5lf seconds", t2 - t1);
		return;
	}
	mpi_log(2, 0, "iteration %d, waiting on barrier %1.5lf seconds", iteration, t2 - t1);
}

MultipleAlignment *mpi_align(PDBChain **chains, int n_als, int nthreads)
{
	double start_time, end_time;
	start_time = MPI_Wtime();

	int rank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	MultipleAlignment **alignments = mpi_init_als(chains, n_als);

	processed_t *aligned = processed_init();
	pairscores_t *scores = NULL;

	double t1, t2;
	t1 = MPI_Wtime();
	mpi_pairwise(alignments, n_als, &scores, aligned, nthreads, rank, nprocs);
	mpi_barrier_wrapper(-1);
	mpi_pairscores_gather(&scores, rank, nprocs); /* sync */
	t2 = MPI_Wtime();
	mpi_log(0, 1, "pairwise time: %1.5lf seconds", t2 - t1);

	for (int i = 0; n_als > 1; i++) {
		t1 = MPI_Wtime();
		MultipleAlignment *ma = mpi_bcast_best_ma(aligned, rank); /* sync */
		mpi_cleanup_best(alignments, &n_als, aligned, ma);

		mpi_iterative(alignments, n_als, aligned, ma, &scores, i, nthreads, rank, nprocs);

		mpi_barrier_wrapper(i);
		alignments[n_als++] = ma;
		t2 = MPI_Wtime();
		mpi_log(0, 1, "iteration %d time: %1.5lf seconds", i, t2 - t1);
	}
	MultipleAlignment *ret = alignments[0];
	pairscores_free(&scores);
	processed_destroy(aligned);
	free(alignments);
	end_time = MPI_Wtime();
	mpi_log(0, 1, "alignment done in %1.5lf seconds", end_time - start_time);
	if (rank == MASTER) {
		return ret;
	}
	CleanupAlignment(ret);
	return NULL;
}
