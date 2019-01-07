/*
 *  Copyright (c) 2017 Maksim Shegay.
 *  The faculty of Computational Mathematics and Cybernetics of Lomonosov Moscow State University
 *  All rights reserved.
 *
 *  parMatt is licensed under the GNU public license version 2.0.
 */

#include <assert.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#include "mpi_pool.h"

const mpi_task_t TERMINATE_TASK = {
    .al1 = -2,
    .al2 = -2,
};

static void taskpool_push(taskpool_t *p, mpi_task_t task)
{
	struct node_t *next = calloc(1, sizeof(struct node_t));
	assert(next);
	next->task = task;
	if (p->head == NULL) {
		p->head = next;
		p->tail = next;
	} else {
		p->tail->next = next;
		p->tail = next;
	}
	p->len++;
}

static mpi_task_t taskpool_pop(taskpool_t *p)
{
	assert(p->len > 0);
	mpi_task_t ret = p->head->task;
	struct node_t *tmp = p->head;
	p->head = p->head->next;
	free(tmp);
	p->len--;
	return ret;
}

static void *fetcher(void *context)
{
	taskpool_ctx_t *ctx = context;
	taskpool_t *p = ctx->pool;
	for (;;) {
		pthread_mutex_lock(p->lock);
		while (p->len > p->nworkers) {
			pthread_cond_wait(p->more_tasks, p->lock);
		}
		pthread_mutex_unlock(p->lock);
		int n;
		int terminate = 0;
		mpi_task_t *tasks = p->fetch(&n);
		terminate = n == -1;
		if (terminate) {
			free(tasks);
			tasks = calloc(p->nworkers, sizeof(mpi_task_t));
			assert(tasks);
			for (int i = 0; i < p->nworkers; i++) {
				tasks[i] = TERMINATE_TASK;
			}
			n = p->nworkers;
		}
		pthread_mutex_lock(p->lock);
		for (int i = 0; i < n; i++) {
			taskpool_push(p, tasks[i]);
		}
		pthread_mutex_unlock(p->lock);
		pthread_cond_broadcast(p->not_empty);
		free(tasks);
		if (terminate) {
			break;
		}
	}
	return NULL;
}

static void *executor(void *context)
{
	taskpool_ctx_t *ctx = context;
	taskpool_t *p = ctx->pool;
	processed_t *processed = ctx->processed;
	for (;;) {
		pthread_mutex_lock(p->lock);
		while (p->len == 0) {
			pthread_cond_wait(p->not_empty, p->lock);
		}
		mpi_task_t t = taskpool_pop(p);
		pthread_mutex_unlock(p->lock);
		if (t.al1 == TERMINATE_TASK.al1 && t.al2 == TERMINATE_TASK.al2) {
			break;
		}
		void *res = p->execute(t);
		processed_append(processed, res);
		
		pthread_mutex_lock(p->lock);
		if (p->len < p->nworkers) {
			pthread_cond_signal(p->more_tasks);
		}
		pthread_mutex_unlock(p->lock);
	}
	return NULL;
}

taskpool_t *taskpool_init(int n, fetch_cb_t fetch, execute_cb_t execute)
{
	taskpool_t *p = calloc(1, sizeof(taskpool_t));
	assert(p);
	p->len = 0;
	p->nworkers = n;

	p->lock = malloc(sizeof(pthread_mutex_t));
	assert(p->lock);
	pthread_mutex_init(p->lock, NULL);

	p->not_empty = malloc(sizeof(pthread_cond_t));
	assert(p->not_empty);
	pthread_cond_init(p->not_empty, NULL);

	p->more_tasks = malloc(sizeof(pthread_cond_t));
	assert(p->more_tasks);
	pthread_cond_init(p->more_tasks, NULL);

	p->fetch = fetch;
	p->execute = execute;
	return p;
}

void taskpool_destroy(taskpool_t *p)
{
	while (p->len > 0) {
		taskpool_pop(p);
	}
	pthread_mutex_destroy(p->lock);
	free(p->lock);
	pthread_cond_destroy(p->not_empty);
	free(p->not_empty);
	pthread_cond_destroy(p->more_tasks);
	free(p->more_tasks);
	free(p);
}

void taskpool_run(int nworkers, taskpool_ctx_t ctx, fetch_cb_t fetch, execute_cb_t execute)
{
	taskpool_t *p = taskpool_init(nworkers, fetch, execute);
	ctx.pool = p;
	pthread_t executors[nworkers];
	pthread_t fetcher_pthtread;

	for (int i = 0; i < nworkers; i++) {
		pthread_create(executors + i, NULL, executor, &ctx);
	}
	pthread_create(&fetcher_pthtread, NULL, fetcher, &ctx);

	for (int i = 0; i < nworkers; i++) {
		pthread_join(executors[i], NULL);
	}
	pthread_join(fetcher_pthtread, NULL);

	taskpool_destroy(p);
}

processed_t *processed_init()
{
	processed_t *p = calloc(1, sizeof(processed_t));
	assert(p);
	p->len = 0;
	p->cap = 1;
	p->data = calloc(1, sizeof(void *));
	assert(p->data);
	p->lock = malloc(sizeof(pthread_mutex_t));
	assert(p->lock);
	pthread_mutex_init(p->lock, NULL);
	return p;
}

void processed_destroy(processed_t *p)
{
	free(p->data);
	pthread_mutex_destroy(p->lock);
	free(p->lock);
	free(p);
}

void processed_append(processed_t *p, void *x)
{
	pthread_mutex_lock(p->lock);
	if (p->len >= p->cap) {
		p->cap *= 2;
		p->data = realloc(p->data, p->cap * sizeof(void *));
		assert(p->data);
	}
	p->data[p->len++] = x;
	pthread_mutex_unlock(p->lock);
}
