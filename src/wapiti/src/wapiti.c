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
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "decoder.h"
#include "model.h"
#include "options.h"
#include "progress.h"
#include "quark.h"
#include "reader.h"
#include "sequence.h"
#include "tools.h"
#include "trainers.h"
#include "wapiti.h"

/*******************************************************************************
 * Training
 ******************************************************************************/
static void trn_auto(mdl_t *mdl) {
	const uint32_t maxiter = mdl->opt->maxiter;
	mdl->opt->maxiter = 3;
	trn_sgdl1(mdl);
	mdl->opt->maxiter = maxiter;
	trn_lbfgs(mdl);
}

static const char *typ_lst[] = {
	"maxent",
	"memm",
	"crf"
};
static const uint32_t typ_cnt = sizeof(typ_lst) / sizeof(typ_lst[0]);

static const struct {
	char *name;
	void (* train)(mdl_t *mdl);
} trn_lst[] = {
	{"l-bfgs", trn_lbfgs},
	{"sgd-l1", trn_sgdl1},
	{"bcd",    trn_bcd  },
	{"rprop",  trn_rprop},
	{"rprop+", trn_rprop},
	{"rprop-", trn_rprop},
	{"auto",   trn_auto }
};
static const uint32_t trn_cnt = sizeof(trn_lst) / sizeof(trn_lst[0]);


/*******************************************************************************
 * Labeling
 ******************************************************************************/



static void dolabel(mdl_t *mdl, char *model, char *input, char *output) {
	// First, load the model provided by the user. This is mandatory to
	// label new datas ;-)
	if (model == NULL)
		fatal("you must specify a model");
	info("* Load model\n");
	FILE *file = fopen(model, "r");
	if (file == NULL)
		pfatal("cannot open input model file");
	mdl_load(mdl, file);
	// Open input and output files
	FILE *fin = stdin, *fout = stdout;
	if (input != NULL) {
		fin = fopen(input, "r");
		if (fin == NULL)
			pfatal("cannot open input data file");
	}
	if (output != NULL) {
		fout = fopen(output, "w");
		if (fout == NULL)
			pfatal("cannot open output data file");
	}
	// Do the labelling
	info("* Label sequences\n");
	tag_label(mdl, fin, fout);
	info("* Done\n");
	// And close files
	if (input != NULL)
		fclose(fin);
	if (output != NULL)
		fclose(fout);
}

/*******************************************************************************
 * Entry point
 ******************************************************************************/
int main(int argc, char *argv[argc]) {
       	mdl_t *mdl = mdl_new(rdr_new(false));
        dolabel(mdl,"model","test.txt","out.txt");
	mdl_free(mdl);
	return EXIT_SUCCESS;
}

