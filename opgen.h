#ifndef _OPGEN_H
#define _OPGEN_H

#include <string.h>
#include <errno.h>
#include <search.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "pcflib.h"
#include "opdefs.h"

const char * skip_to_colon(const char * line);

const char * assert_token(const char * line, char * buf, char * bitr, const char * token);

const char * read_token(const char * line, char * bitr);

PCFOP * read_label(const char * line, struct PCFState * st, uint32_t iptr);

PCFOP * read_initbase(const char * line);

PCFOP * read_gate(const char * line);

PCFOP * read_copy(const char * line);

PCFOP * read_arith(const char * line, void (*op)(struct PCFState *, struct PCFOP *));

PCFOP * read_add(const char * line);
PCFOP * read_mul(const char * line);

PCFOP * read_copy_indir(const char * line);

PCFOP * read_indir_copy(const char * line);

PCFOP * read_const(const char * line);

int count_tokens_to_close_paren(const char * line);

PCFOP * read_bits(const char * line);

PCFOP * read_join(const char * line);

PCFOP * read_mkptr(const char * line);

PCFOP * read_call(const char * line);

PCFOP * read_branch(const char * line);

PCFOP * read_ret(const char * line);

PCFOP * read_clear(const char * line); 

PCFOP * read_instr(struct PCFState * st, const char * line, uint32_t iptr);


#endif
