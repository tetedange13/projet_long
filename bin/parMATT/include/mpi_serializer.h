/*
 *  Copyright (c) 2017 Maksim Shegay.
 *  The faculty of Computational Mathematics and Cybernetics of Lomonosov Moscow State University
 *  All rights reserved.
 *
 *  parMatt is licensed under the GNU public license version 2.0.
 */

#ifndef SERIALIZER_H
#define SERIALIZER_H

#include "MultipleAlignment.h"
#include "AssemblyOrder.h"

void *mpi_ma_serialize(MultipleAlignment *ma, int *size);

MultipleAlignment *mpi_ma_deserialize(char *buf);

#endif
