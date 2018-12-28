/*
 *  Copyright (c) 2017 Maksim Shegay.
 *  The faculty of Computational Mathematics and Cybernetics of Lomonosov Moscow State University
 *  All rights reserved.
 *
 *  parMatt is licensed under the GNU public license version 2.0.
 */

#include <assert.h>

#include "mpi_serializer.h"

static void memcpy_wrapper(char *dst, int *offset, void *src, size_t size)
{
	memcpy(dst + *offset, src, size);
	*offset += size;
}

static void *memdup_wrapper(char *src, int *offset, size_t size)
{
	void *out = malloc(size);
	assert(out || !size);
	memcpy(out, src + *offset, size);
	*offset += size;
	return out;
}

static char deserialize_char_value(char *src, int *offset)
{
	char ret = *(src + *offset);
	*offset += sizeof(char);
	return ret;
}

static int deserialize_int_value(char *src, int *offset)
{
	int ret = *((int *)(src + *offset));
	*offset += sizeof(int);
	return ret;
}

static double deserialize_double_value(char *src, int *offset)
{
	double ret = *((double *)(src + *offset));
	*offset += sizeof(double);
	return ret;
}

static size_t pdb_chain_size(PDBChain *ma)
{
	size_t size = 14 * sizeof(int) + 2;
	size += ma->numCisPeps * sizeof(CisPep);
	size += ma->numSeqAdvs * sizeof(SeqAdv);
	size += ma->numDbrefs * sizeof(DBRef);
	size += ma->numSSbonds * sizeof(SSbond);

	for (int i = 0; i < ma->numBetaSheets; i++) {
		size += 4 + sizeof(int);
		size += ma->betaSheets->numStrands * sizeof(BetaSheetPair);
	}

	size += ma->numAlphaHelices * sizeof(AlphaHelix);
	size += ma->numHbonds * sizeof(HydrogenBond);
	size += ma->numBetaPairs * sizeof(BetaPair);
	size += ma->numAtoms * sizeof(Atom);
	size += ma->numTempResidues * sizeof(TempResidue);
	size += ma->length * sizeof(Residue);
	size += strlen(ma->idString) + 1;
	size += strlen(ma->seq->name) + 1;
	size += ma->seq->length;
	return size;
}

static int mpi_get_buf_size(MultipleAlignment *ma)
{
	int size = 3 * sizeof(double) + 3 * sizeof(int);

	size += ma->numBlocks * sizeof(AlignedBlock);
	size += ma->numChains * (ma->numResidues * sizeof(ResiduePosition) + sizeof(double));
	for (int i = 0; i < ma->numChains; i++) {
		size += 2 * sizeof(int);
		size += ma->chains[i]->length * sizeof(ResiduePosition);
		size += pdb_chain_size(ma->chains[i]->pdb);
	}

	size += sizeof(int) + 3 * ma->order->nodes * sizeof(int);
	return size;
}

static void assembly_order_serialize(char *buf, int *offset, AssemblyOrder *order)
{
	size_t size = sizeof(int) + 3 * order->nodes * sizeof(int);
	int *tree = malloc(size);
	assert(tree);
	tree[0] = order->nodes;

	/* AssemblyOrder contains all the nodes in a continguous memory area */
	/* located right after AssemblyOrder */
	AssemblyNode *root = (AssemblyNode *)(order + 1);
	for (int i = 0; i < order->nodes; i++) {
		AssemblyNode *node = (AssemblyNode *)(order->root + i);
		tree[3 * i + 1] = node->id;
		tree[3 * i + 2] = node->left ? node->left - root : 0;
		tree[3 * i + 3] = node->right ? node->right - root : 0;
	}
	memcpy_wrapper(buf, offset, tree, size);
	free(tree);
}

static void assembly_order_deserialize(MultipleAlignment *ma, char *buf, int *offset)
{
	int nodes = deserialize_int_value(buf, offset);
	AssemblyOrder *order = malloc(sizeof(AssemblyOrder) + nodes * sizeof(AssemblyNode));
	assert(order);
	order->nodes = nodes;
	order->root = (AssemblyNode *)(order + 1);
	for (int i = 0; i < nodes; i++) {
		AssemblyNode *node = (AssemblyNode *)(order->root + i);

		int id = deserialize_int_value(buf, offset);
		int left = deserialize_int_value(buf, offset);
		int right = deserialize_int_value(buf, offset);

		node->id = id;
		node->left = left ? order->root + left : NULL;
		node->right = right ? order->root + right : NULL;
	}
	ma->order = order;
}

static void weighted_residue_positions_serialize(char *buf, int *offset, MultipleAlignment *ma)
{
	memcpy_wrapper(buf, offset, &ma->numResidues, sizeof(int));
	if (ma->numResidues == 0)
		return;
	for (int i = 0; i < ma->numChains; i++) {
		memcpy_wrapper(buf, offset, ma->residues[i].res, ma->numResidues * sizeof(ResiduePosition));
		memcpy_wrapper(buf, offset, &ma->residues[i].weight, sizeof(double));
	}
}

static void weighted_residue_positions_deserialize(MultipleAlignment *ma, char *buf, int *offset, int numChains)
{
	ma->numResidues = deserialize_int_value(buf, offset);
	if (ma->numResidues == 0){
		ma->residues = NULL;
		return;
	}
	ma->residues = malloc(numChains * sizeof(WeightedResiduePositions));
	assert(ma->residues || !numChains);
	for (int i = 0; i < numChains; i++) {
		ma->residues[i].res = memdup_wrapper(buf, offset, ma->numResidues * sizeof(ResiduePosition));
		ma->residues[i].weight = deserialize_double_value(buf, offset);
	}
}

static void pdb_chain_serialize(char *buf, int *offset, PDBChain *chain)
{
	memcpy_wrapper(buf, offset, &chain->numCisPeps, sizeof(int));
	memcpy_wrapper(buf, offset, chain->cisPeps, chain->numCisPeps * sizeof(CisPep));

	memcpy_wrapper(buf, offset, &chain->numSeqAdvs, sizeof(int));
	memcpy_wrapper(buf, offset, chain->seqAdvs, chain->numSeqAdvs * sizeof(SeqAdv));

	memcpy_wrapper(buf, offset, &chain->numDbrefs, sizeof(int));
	memcpy_wrapper(buf, offset, chain->dbrefs, chain->numDbrefs * sizeof(DBRef));

	memcpy_wrapper(buf, offset, &chain->numSSbonds, sizeof(int));
	memcpy_wrapper(buf, offset, chain->ssbonds, chain->numSSbonds * sizeof(SSbond));

	memcpy_wrapper(buf, offset, &chain->numBetaSheets, sizeof(int));

	for (int i = 0; i < chain->numBetaSheets; i++) {
		memcpy_wrapper(buf, offset, chain->betaSheets->id, 4);
		memcpy_wrapper(buf, offset, &chain->betaSheets->numStrands, sizeof(int));
		memcpy_wrapper(buf, offset, chain->betaSheets->strands, chain->betaSheets->numStrands * sizeof(BetaSheetPair));
	}

	memcpy_wrapper(buf, offset, &chain->numAlphaHelices, sizeof(int));
	memcpy_wrapper(buf, offset, chain->alphaHelices, chain->numAlphaHelices * sizeof(AlphaHelix));

	memcpy_wrapper(buf, offset, &chain->numHbonds, sizeof(int));
	memcpy_wrapper(buf, offset, chain->hbonds, chain->numHbonds * sizeof(HydrogenBond));

	memcpy_wrapper(buf, offset, &chain->numBetaPairs, sizeof(int));
	memcpy_wrapper(buf, offset, chain->betaPairs, chain->numBetaPairs * sizeof(BetaPair));

	memcpy_wrapper(buf, offset, &chain->length, sizeof(int));

	memcpy_wrapper(buf, offset, &chain->numAtoms, sizeof(int));
	memcpy_wrapper(buf, offset, chain->atoms, chain->numAtoms * sizeof(Atom));

	memcpy_wrapper(buf, offset, &chain->chainName, sizeof(char));
	memcpy_wrapper(buf, offset, &chain->tempAtoms, sizeof(char));

	memcpy_wrapper(buf, offset, &chain->terminated, sizeof(int));
	memcpy_wrapper(buf, offset, &chain->numTempResidues, sizeof(int));
	memcpy_wrapper(buf, offset, chain->tempResidues, chain->numTempResidues * sizeof(TempResidue));

	memcpy_wrapper(buf, offset, &chain->secondaryCalculated, sizeof(int));

	memcpy_wrapper(buf, offset, chain->residues, chain->length * sizeof(Residue));

	size_t size = strlen(chain->idString) + 1;
	memcpy_wrapper(buf, offset, chain->idString, size);

	size = strlen(chain->seq->name) + 1;
	memcpy_wrapper(buf, offset, chain->seq->name, size);
	memcpy_wrapper(buf, offset, &chain->seq->length, sizeof(int));
	memcpy_wrapper(buf, offset, chain->seq->seq, chain->seq->length);
}

static PDBChain *pdb_chain_deserialize(char *buf, int *offset)
{
	PDBChain *chain = malloc(sizeof(PDBChain));
	assert(chain);
	chain->numCisPeps = deserialize_int_value(buf, offset);
	chain->cisPeps = memdup_wrapper(buf, offset, chain->numCisPeps * sizeof(CisPep));

	chain->numSeqAdvs = deserialize_int_value(buf, offset);
	chain->seqAdvs = memdup_wrapper(buf, offset, chain->numSeqAdvs * sizeof(SeqAdv));

	chain->numDbrefs = deserialize_int_value(buf, offset);
	chain->dbrefs = memdup_wrapper(buf, offset, chain->numDbrefs * sizeof(DBRef));

	chain->numSSbonds = deserialize_int_value(buf, offset);
	chain->ssbonds = memdup_wrapper(buf, offset, chain->numSSbonds * sizeof(SSbond));

	chain->numBetaSheets = deserialize_int_value(buf, offset);

	chain->betaSheets = malloc(chain->numBetaSheets * sizeof(BetaSheet));
	assert(chain->betaSheets || !chain->numBetaSheets);
	for (int i = 0; i < chain->numBetaSheets; i++) {
		memcpy(chain->betaSheets[i].id, buf + *offset, 4);
		*offset += 4;
		chain->betaSheets[i].numStrands = deserialize_int_value(buf, offset);
		chain->betaSheets[i].strands = memdup_wrapper(buf, offset, chain->betaSheets[i].numStrands * sizeof(BetaSheetPair));
	}

	chain->numAlphaHelices = deserialize_int_value(buf, offset);
	chain->alphaHelices = memdup_wrapper(buf, offset, chain->numAlphaHelices * sizeof(AlphaHelix));

	chain->numHbonds = deserialize_int_value(buf, offset);
	chain->hbonds = memdup_wrapper(buf, offset, chain->numHbonds * sizeof(HydrogenBond));

	chain->numBetaPairs = deserialize_int_value(buf, offset);
	chain->betaPairs = memdup_wrapper(buf, offset, chain->numBetaPairs * sizeof(BetaPair));

	chain->length = deserialize_int_value(buf, offset);

	chain->numAtoms = deserialize_int_value(buf, offset);
	chain->atoms = memdup_wrapper(buf, offset, chain->numAtoms * sizeof(Atom));

	chain->chainName = deserialize_char_value(buf, offset);
	chain->tempAtoms = deserialize_char_value(buf, offset);

	chain->terminated = deserialize_int_value(buf, offset);
	chain->numTempResidues = deserialize_int_value(buf, offset);
	chain->tempResidues = memdup_wrapper(buf, offset, chain->numTempResidues * sizeof(TempResidue));

	chain->secondaryCalculated = deserialize_int_value(buf, offset);

	chain->residues = memdup_wrapper(buf, offset, chain->length * sizeof(Residue));

	size_t size = strlen(buf + *offset) + 1;
	chain->idString = memdup_wrapper(buf, offset, size);

	chain->seq = malloc(sizeof(Sequence));
	assert(chain->seq);

	size = strlen(buf + *offset) + 1;
	chain->seq->name = memdup_wrapper(buf, offset, size);
	chain->seq->length = deserialize_int_value(buf, offset);
	chain->seq->seq = memdup_wrapper(buf, offset, chain->seq->length);

	return chain;
}

static void residue_positions_serialize(char *buf, int *offset, MultipleAlignment *ma)
{
	for (int i = 0; i < ma->numChains; i++) {
		memcpy_wrapper(buf, offset, &ma->chains[i]->length, sizeof(int));
		memcpy_wrapper(buf, offset, &ma->chains[i]->id, sizeof(int));
		memcpy_wrapper(buf, offset, ma->chains[i]->res, sizeof(ResiduePosition) * ma->chains[i]->length);
		pdb_chain_serialize(buf, offset, ma->chains[i]->pdb);
	}
}

static void residue_positions_deserialize(MultipleAlignment *ma, char *buf, int *offset, int numChains)
{
	ma->chains = malloc(numChains * sizeof(ResiduePositions *));
	assert(ma->chains || !numChains);
	for (int i = 0; i < numChains; i++) {
		ma->chains[i] = malloc(sizeof(ResiduePositions));
		assert(ma->chains[i]);
		int length = deserialize_int_value(buf, offset);
		ma->chains[i]->length = length;

		ma->chains[i]->id = deserialize_int_value(buf, offset);

		ma->chains[i]->res = memdup_wrapper(buf, offset, length * sizeof(ResiduePosition));
		ma->chains[i]->pdb = pdb_chain_deserialize(buf, offset);
	}
}

void *mpi_ma_serialize(MultipleAlignment *ma, int *size)
{
	*size = mpi_get_buf_size(ma);
	int offset = 0;
	char *buf = malloc(*size);
	assert(buf || !(*size));

	memcpy_wrapper(buf, &offset, &ma->score, sizeof(double));
	memcpy_wrapper(buf, &offset, &ma->rmsd, sizeof(double));
	memcpy_wrapper(buf, &offset, &ma->pvalue, sizeof(double));
	memcpy_wrapper(buf, &offset, &ma->numChains, sizeof(int));
	memcpy_wrapper(buf, &offset, &ma->numBlocks, sizeof(int));
	memcpy_wrapper(buf, &offset, ma->blocks, ma->numBlocks * sizeof(AlignedBlock));

	weighted_residue_positions_serialize(buf, &offset, ma);
	residue_positions_serialize(buf, &offset, ma);
	assembly_order_serialize(buf, &offset, ma->order);
	return buf;
}

MultipleAlignment *mpi_ma_deserialize(char *buf)
{
	int offset = 0;
	MultipleAlignment *ma = malloc(sizeof(MultipleAlignment));
	assert(ma);

	ma->score = deserialize_double_value(buf, &offset);
	ma->rmsd = deserialize_double_value(buf, &offset);
	ma->pvalue = deserialize_double_value(buf, &offset);
	ma->numChains = deserialize_int_value(buf, &offset);
	ma->numBlocks = deserialize_int_value(buf, &offset);

	ma->blocks = memdup_wrapper(buf, &offset, ma->numBlocks * sizeof(AlignedBlock));

	weighted_residue_positions_deserialize(ma, buf, &offset, ma->numChains);
	residue_positions_deserialize(ma, buf, &offset, ma->numChains);
	assembly_order_deserialize(ma, buf, &offset);

	/* following fields will be recalculated upon alignment */
	ma->conflictMap = NULL;
	ma->averages = NULL;
	return ma;
}
