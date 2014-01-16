/*
 *      Wapiti - A linear-chain CRF tool
 *
 * Copyright (c) 2009-2011  CNRS
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <inttypes.h>
#include <float.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include "wapiti.h"
#include "gradient.h"
#include "model.h"
#include "quark.h"
#include "reader.h"
#include "sequence.h"
#include "thread.h"
#include "tools.h"
#include "decoder.h"
#include "vmath.h"
#include "dual.h"

#include <math.h>

/******************************************************************************
 * Sequence tagging
 *
 *   This module implement sequence tagging using a trained model and model
 *   evaluation on devlopment set.
 *
 *   The viterbi can be quite intensive on the stack if you push in it long
 *   sequence and use large labels set. It's less a problem than in gradient
 *   computations but it can show up in particular cases. The fix is to call it
 *   through the mth_spawn function and request enough stack space, this will be
 *   fixed in next version.
 ******************************************************************************/

/* tag_expsc:
 *   Compute the score lattice for classical Viterbi decoding. This is the same
 *   as for the first step of the gradient computation with the exception that
 *   we don't need to take the exponential of the scores as the Viterbi decoding
 *   works in log-space.
 */
static int tag_expsc(mdl_t *mdl, const seq_t *seq, double *vpsi,dual_t* dual) {
	const double  *x = mdl->theta;
	const uint32_t Y = mdl->nlbl;
	const uint32_t T = seq->len;
	double (*psi)[T][Y][Y] = (void *)vpsi;


	// We first have to compute the Ψ_t(y',y,x_t) weights defined as
	//   Ψ_t(y',y,x_t) = \exp( ∑_k θ_k f_k(y',y,x_t) )
	// So at position 't' in the sequence, for each couple (y',y) we have
	// to sum weights of all features.
	// This is the same than what we do for computing the gradient but, as
	// the viterbi algorithm also work in the logarithmic space, we can
	// remove the exponential.
	//
	// Only the observations present at this position will have a non-nul
	// weight so we can sum only on thoses.
	//
	// As we use only two kind of features: unigram and bigram, we can
	// rewrite this as
	//   ∑_k μ_k(y, x_t) f_k(y, x_t) + ∑_k λ_k(y', y, x_t) f_k(y', y, x_t)
	// Where the first sum is over the unigrams features and the second is
	// over bigrams ones.
	//
	// This allow us to compute Ψ efficiently in two steps
	//   1/ we sum the unigrams features weights by looping over actives
	//        unigrams observations. (we compute this sum once and use it
	//        for each value of y')
	//   2/ we add the bigrams features weights by looping over actives
	//        bigrams observations (we don't have to do this for t=0 since
	//        there is no bigrams here)
	for (uint32_t t = 0; t < T; t++) {
		const pos_t *pos = &(seq->pos[t]);
		for (uint32_t y = 0; y < Y; y++) {
			double sum = 0.0;
			for (uint32_t n = 0; n < pos->ucnt; n++) {
				const uint64_t o = pos->uobs[n];
				sum += x[mdl->uoff[o] + y];
			}
                        sum += dual_get_unary_penalty(dual,t,y);
                        //printf("(%d,%d):%f\n",t,y,dual_get_penalty(d,t,y));
			for (uint32_t yp = 0; yp < Y; yp++)
                          (*psi)[t][yp][y] = sum;

		}
	}
	for (uint32_t t = 1; t < T; t++) {
		const pos_t *pos = &(seq->pos[t]);
		for (uint32_t yp = 0, d = 0; yp < Y; yp++) {
			for (uint32_t y = 0; y < Y; y++, d++) {
				double sum = 0.0;

                                /* double tmp = dual_get_binary_penalty(dual,t,yp,y); */
                                /* if (tmp != 0.0) */
                                /* { */
                                /*   fprintf(stderr, "penalty detected: %f\n", tmp); */
                                /*   fprintf(stderr, "previous tag: %s\n", qrk_id2str(mdl->reader->lbl, yp)); */
                                /*   fprintf(stderr, "current tag: %s\n", qrk_id2str(mdl->reader->lbl, y)); */
                                /* } */

                                sum += dual_get_binary_penalty(dual,t,yp,y);

				for (uint32_t n = 0; n < pos->bcnt; n++) {
					const uint64_t o = pos->bobs[n];
					sum += x[mdl->boff[o] + d];
                                        //sum += x[mdl->boff[o] + d] + dual_get_binary_penalty(dual,t,yp,y);
				}
				(*psi)[t][yp][y] += sum;
			}
		}
	}

	return 0;
}

/* tag_postsc:
 *   This function compute score lattice with posteriors. This generally result
 *   in a slightly best labelling and allow to output normalized score for the
 *   sequence and for each labels but this is more costly as we have to perform
 *   a full forward backward instead of just the forward pass.
 */
static int tag_postsc(mdl_t *mdl, const seq_t *seq, double *vpsi) {
	const uint32_t Y = mdl->nlbl;
	const uint32_t T = seq->len;
	double (*psi)[T][Y][Y] = (void *)vpsi;
	grd_st_t *grd_st = grd_stnew(mdl, NULL);
	grd_st->first = 0;
	grd_st->last  = T - 1;
	grd_stcheck(grd_st, seq->len);
	if (mdl->opt->sparse) {
		grd_spdopsi(grd_st, seq);
		grd_spfwdbwd(grd_st, seq);
	} else {
		grd_fldopsi(grd_st, seq);
		grd_flfwdbwd(grd_st, seq);
	}
	double (*alpha)[T][Y] = (void *)grd_st->alpha;
	double (*beta )[T][Y] = (void *)grd_st->beta;
	double  *unorm        =         grd_st->unorm;
	for (uint32_t t = 0; t < T; t++) {
		for (uint32_t y = 0; y < Y; y++) {
			double e = (*alpha)[t][y] * (*beta)[t][y] * unorm[t];
			for (uint32_t yp = 0; yp < Y; yp++)
				(*psi)[t][yp][y] = e;
		}
	}
	grd_stfree(grd_st);
	return 1;
}

/* tag_viterbi:
 *   This function implement the Viterbi algorithm in order to decode the most
 *   probable sequence of labels according to the model. Some part of this code
 *   is very similar to the computation of the gradient as expected.
 *
 *   And like for the gradient, the caller is responsible to ensure there is
 *   enough stack space.
 */
void tag_viterbi(mdl_t *mdl, const seq_t *seq,
	         uint32_t out[], double *score, double psc[], dual_t* d)
{
	const uint32_t Y = mdl->nlbl;
	const uint32_t T = seq->len;

	double   *vpsi  = xmalloc(sizeof(double  ) * T * Y * Y);
	uint32_t *vback = xmalloc(sizeof(uint32_t) * T * Y);
	double   (*psi) [T][Y][Y] = (void *)vpsi;
	uint32_t (*back)[T][Y]    = (void *)vback;

	double *cur = xmalloc(sizeof(double) * Y);
	double *old = xmalloc(sizeof(double) * Y);
	// We first compute the scores for each transitions in the lattice of
	// labels.
	int op;
	op = tag_expsc(mdl, seq, vpsi, d);


        //fprintf(stderr, "T = %d\n", T);

        mdl->opt = &opt_defaults;
        //mdl->opt->lblpost = true;
        //mdl->opt->sparse = true;

        /* fprintf(stderr, "Before tag_postsc\n"); */
        /* op = tag_postsc(mdl,seq,vpsi); */
        /* fprintf(stderr, "After tag_postsc\n"); */

	// Now we can do the Viterbi algorithm. This is very similar to the
	// forward pass
	//   | α_1(y) = Ψ_1(y,x_1)
	//   | α_t(y) = max_{y'} α_{t-1}(y') + Ψ_t(y',y,x_t)
	// We just replace the sum by a max and as we do the computation in the
	// logarithmic space the product become a sum. (this also mean that we
	// don't have to worry about numerical problems)
	//
	// Next we have to walk backward over the α in order to find the best
	// path. In order to do this efficiently, we keep in the 'back' array
	// the indice of the y value selected by the max. This also mean that
	// we only need the current and previous value of the α vectors, not
	// the full matrix.
	for (uint32_t y = 0; y < Y; y++)
          if (op)
            //cur[y] = log((*psi)[0][0][y]);
            cur[y] = (*psi)[0][0][y];
          else
            cur[y] = (*psi)[0][0][y];
	for (uint32_t t = 1; t < T; t++) {
		for (uint32_t y = 0; y < Y; y++)
			old[y] = cur[y];
		for (uint32_t y = 0; y < Y; y++) {
			double   bst = -DBL_MAX;
			uint32_t idx = 0;
			for (uint32_t yp = 0; yp < Y; yp++) {
				double val = old[yp];
				if (op)
                                {
                                  //fprintf(stderr, "%f\n", (*psi)[t][yp][y]);
                                  val *= (*psi)[t][yp][y];
                                  //val += dual_get_binary_penalty(d,t,yp,y);
                                }
				else
                                {
                                  //fprintf(stderr, "%f\n", (*psi)[t][yp][y]);
                                  val += (*psi)[t][yp][y];
                                  //val += dual_get_binary_penalty(d,t,yp,y);
                                }
				if (val > bst) {
					bst = val;
					idx = yp;
				}
			}
			(*back)[t][y] = idx;
			cur[y]        = bst;
		}
	}
	// We can now build the sequence of labels predicted by the model. For
	// this we search in the last α vector the best value. Using this index
	// as a starting point in the back-pointer array we finally can decode
	// the best sequence.
	uint32_t bst = 0;
	for (uint32_t y = 1; y < Y; y++)
		if (cur[y] > cur[bst])
                  bst = y;

        //        fprintf(stderr, "bst: %f\n", cur[bst]);

	if (score != NULL)
		*score = cur[bst];
	for (uint32_t t = T; t > 0; t--) {
		const uint32_t yp = (t != 1) ? (*back)[t - 1][bst] : 0;
		const uint32_t y  = bst;
		out[t - 1] = y;
		if (psc != NULL)
			psc[t - 1] = (*psi)[t - 1][yp][y];
		bst = yp;
	}
	free(old);
	free(cur);
	free(vback);
	free(vpsi);
}

/* tag_nbviterbi:
 *   This function implement the Viterbi algorithm in order to decode the N-most
 *   probable sequences of labels according to the model. It can be used to
 *   compute only the best one and will return the same sequence than the
 *   previous function but will be slower to do it.
 */
void tag_nbviterbi(mdl_t *mdl, const seq_t *seq, uint32_t N,
                   uint32_t out[][N], double sc[], double psc[][N]) {
	const uint32_t Y = mdl->nlbl;
	const uint32_t T = seq->len;
	double   *vpsi  = xmalloc(sizeof(double  ) * T * Y * Y);
	uint32_t *vback = xmalloc(sizeof(uint32_t) * T * Y * N);
	double   (*psi) [T][Y    ][Y] = (void *)vpsi;
	uint32_t (*back)[T][Y * N]    = (void *)vback;
	double *cur = xmalloc(sizeof(double) * Y * N);
	double *old = xmalloc(sizeof(double) * Y * N);
	// We first compute the scores for each transitions in the lattice of
	// labels.
	int op;
	if (mdl->opt->lblpost)
		op = tag_postsc(mdl, seq, (double *)psi);
	else
		op = tag_expsc(mdl, seq, (double *)psi,NULL); /// WARNING SHOULD NOT WORK!!
	// Here also, it's classical but we have to keep the N best paths
	// leading to each nodes of the lattice instead of only the best one.
	// This mean that code is less trivial and the current implementation is
	// not the most efficient way to do this but it works well and is good
	// enough for the moment.
	// We first build the list of all incoming arcs from all paths from all
	// N-best nodes and next select the N-best one. There is a lot of room
	// here for later optimisations if needed.
	for (uint32_t y = 0, d = 0; y < Y; y++) {
		cur[d++] = (*psi)[0][0][y];
		for (uint32_t n = 1; n < N; n++)
			cur[d++] = -DBL_MAX;
	}
	for (uint32_t t = 1; t < T; t++) {
		for (uint32_t d = 0; d < Y * N; d++)
			old[d] = cur[d];
		for (uint32_t y = 0; y < Y; y++) {
			// 1st, build the list of all incoming
			double lst[Y * N];
			for (uint32_t yp = 0, d = 0; yp < Y; yp++) {
				for (uint32_t n = 0; n < N; n++, d++) {
					lst[d] = old[d];
					if (op)
						lst[d] *= (*psi)[t][yp][y];
					else
						lst[d] += (*psi)[t][yp][y];
				}
			}
			// 2nd, init the back with the N first
			uint32_t *bk = &(*back)[t][y * N];
			for (uint32_t n = 0; n < N; n++)
				bk[n] = n;
			// 3rd, search the N highest values
			for (uint32_t i = N; i < N * Y; i++) {
				// Search the smallest current value
				uint32_t idx = 0;
				for (uint32_t n = 1; n < N; n++)
					if (lst[bk[n]] < lst[bk[idx]])
						idx = n;
				// And replace it if needed
				if (lst[i] > lst[bk[idx]])
					bk[idx] = i;
			}
			// 4th, get the new scores
			for (uint32_t n = 0; n < N; n++)
				cur[y * N + n] = lst[bk[n]];
		}
	}
	// Retrieving the best paths is similar to classical Viterbi except that
	// we have to search for the N bet ones and there is N time more
	// possibles starts.
	for (uint32_t n = 0; n < N; n++) {
		uint32_t bst = 0;
		for (uint32_t d = 1; d < Y * N; d++)
			if (cur[d] > cur[bst])
				bst = d;
		if (sc != NULL)
			sc[n] = cur[bst];
		cur[bst] = -DBL_MAX;
		for (uint32_t t = T; t > 0; t--) {
			const uint32_t yp = (t != 1) ? (*back)[t - 1][bst] / N: 0;
			const uint32_t y  = bst / N;
			out[t - 1][n] = y;
			if (psc != NULL)
				psc[t - 1][n] = (*psi)[t - 1][yp][y];
			bst = (*back)[t - 1][bst];
		}
	}
	free(old);
	free(cur);
	free(vback);
	free(vpsi);
}





/* tag_label:
 *   Label a data file using the current model. This output an almost exact copy
 *   of the input file with an additional column with the predicted label. If
 *   the check option is specified, the input file must be labelled and the
 *   predicted labels will be checked against the provided ones. This will
 *   output error rates during the labelling and detailed statistics per label
 *   at the end.
 */

void tag_label(mdl_t *mdl, FILE *fin, FILE *fout) {

	qrk_t *lbls = mdl->reader->lbl;
	const uint32_t Y = mdl->nlbl;

	uint64_t scnt = 0;
	uint64_t stat[3][Y];
	for (uint32_t y = 0; y < Y; y++)
		stat[0][y] = stat[1][y] = stat[2][y] = 0;
	// Next read the input file sequence by sequence and label them, we have
	// to take care of not discarding the raw input as we want to send it
	// back to the output with the additional predicted labels.

        /* WE READ AND TAG ONLY FIRST SENTENCE */

        //while (!feof(fin)) {
		// So, first read an input sequence keeping the raw_t object
		// available, and label it with Viterbi.
		raw_t *raw = rdr_readraw(mdl->reader, fin);
                dual_t* d = dual_init(5,Y);


                //printf("P+D=%d\n",qrk_str2id(lbls, "P+D")); // print -1 if P+D not found!
                //printf("ADP=%d\n",qrk_str2id(lbls, "ADP"));
                //printf("NOUN=%d\n",qrk_str2id(lbls, "NOUN"));

                //sentence of length 5
                // "NPP", "V", "D", "NC", "PUNCT" (gold)
                // "PRO", "V", "P+D", "NC", "PUNCT" (predicted)
                //dual_set_in_z(d,0,qrk_str2id(lbls, "PRO"));
                //dual_set_in_z(d,1,qrk_str2id(lbls, "V"));
                //dual_set_in_z(d,2,qrk_str2id(lbls, "ADP"));
                //dual_set_in_z(d,3,qrk_str2id(lbls, "NOUN"));
                //dual_set_in_z(d,4,qrk_str2id(lbls, "."));

                dual_set_unary_penalty(d,0,qrk_str2id(lbls, "PRO"),-1);
                dual_set_unary_penalty(d,0,qrk_str2id(lbls, "V"),-1);
                dual_set_unary_penalty(d,0,qrk_str2id(lbls, "DET"),-1);
                dual_set_unary_penalty(d,2,qrk_str2id(lbls, "ADP"),-10);
                //dual_set_in_u(d,3,qrk_str2id(lbls, "NOUN"),-1);
                //dual_set_in_u(d,4,qrk_str2id(lbls, "."),-1);

                //if (raw == NULL)
		//	break;

                seq_t *seq = rdr_raw2seq(mdl->reader, raw, false);
		const uint32_t T = seq->len;
		uint32_t *out = xmalloc(sizeof(uint32_t) * T);
		double   *psc = xmalloc(sizeof(double  ) * T);
		double   *scs = xmalloc(sizeof(double  ));
		tag_viterbi(mdl, seq, (uint32_t*)out, scs, (double*)psc,d);
		// Next we output the raw sequence with an aditional column for
		// the predicted labels
                      for (uint32_t t = 0; t < T; t++) {
				fprintf(fout, "%s\t", raw->lines[t]);
				uint32_t lbl = out[t];
				const char *lblstr = qrk_id2str(lbls, lbl);
				fprintf(fout, "%s", lblstr);

				fprintf(fout, "\n");
			}
			fprintf(fout, "\n");
		fflush(fout);

		// Cleanup memory used for this sequence
		free(scs);
		free(psc);
		free(out);
		rdr_freeseq(seq);
		rdr_freeraw(raw);

		// And report our progress, at regular interval we display how
		// much sequence are labelled and if possible the current tokens
		// and sequence error rates.
		if (++scnt % 1000 == 0) {
			info("%10"PRIu64" sequences labeled", scnt);
			info("\n");
		}
	//}
}

/* eval_t:
 *   This a state tracker used to communicate between the main eval function and
 *   its workers threads, the <mdl> and <dat> fields are used to transmit to the
 *   workers informations needed to make the computation, the other fields are
 *   for returning the partial results.
 */
typedef struct eval_s eval_t;
struct eval_s {
	mdl_t    *mdl;
	dat_t    *dat;
	uint64_t  tcnt;  // Processed tokens count
	uint64_t  terr;  // Tokens error found
	uint64_t  scnt;  // Processes sequences count
	uint64_t  serr;  // Sequence error found
};

/* tag_evalsub:
 *   This is where the real evaluation is done by the workers, we process data
 *   by batch and for each batch do a simple Viterbi and scan the result to find
 *   errors.
 */
static void tag_evalsub(job_t *job, uint32_t id, uint32_t cnt, eval_t *eval) {
	unused(id && cnt);
	mdl_t *mdl = eval->mdl;
	dat_t *dat = eval->dat;
	eval->tcnt = 0;
	eval->terr = 0;
	eval->scnt = 0;
	eval->serr = 0;
	// We just get a job a process all the squence in it.
	uint32_t count, pos;
	while (mth_getjob(job, &count, &pos)) {
		for (uint32_t s = pos; s < pos + count; s++) {
			// Tag the sequence with the viterbi
			const seq_t *seq = dat->seq[s];
			const uint32_t T = seq->len;
			uint32_t out[T];
			tag_viterbi(mdl, seq, out, NULL, NULL,NULL);
			// And check for eventual (probable ?) errors
			bool err = false;
			for (uint32_t t = 0; t < T; t++)
				if (seq->pos[t].lbl != out[t])
					eval->terr++, err = true;
			eval->tcnt += T;
			eval->scnt += 1;
			eval->serr += err;
		}
	}
}

/* tag_eval:
 *   Compute the token error rate and sequence error rate over the devel set (or
 *   taining set if not available).
 */
void tag_eval(mdl_t *mdl, double *te, double *se) {
	const uint32_t W = mdl->opt->nthread;
	dat_t *dat = (mdl->devel == NULL) ? mdl->train : mdl->devel;
	// First we prepare the eval state for all the workers threads, we just
	// have to give them the model and dataset to use. This state will be
	// used to retrieve partial result they computed.
	eval_t *eval[W];
	for (uint32_t w = 0; w < W; w++) {
		eval[w] = xmalloc(sizeof(eval_t));
		eval[w]->mdl = mdl;
		eval[w]->dat = dat;
	}
	// And next, we call the workers to do the job and reduce the partial
	// result by summing them and computing the final error rates.
	mth_spawn((func_t *)tag_evalsub, W, (void *)eval, dat->nseq,
		mdl->opt->jobsize);
	uint64_t tcnt = 0, terr = 0;
	uint64_t scnt = 0, serr = 0;
	for (uint32_t w = 0; w < W; w++) {
		tcnt += eval[w]->tcnt;
		terr += eval[w]->terr;
		scnt += eval[w]->scnt;
		serr += eval[w]->serr;
		free(eval[w]);
	}
	*te = (double)terr / tcnt * 100.0;
	*se = (double)serr / scnt * 100.0;
}
