#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <errno.h>
#include <search.h>
#include <assert.h>

#include "pcflib.h"
#include "opdefs.h"
#include "opflows.h"

void check_alloc(void * ptr)
{
 if(ptr == 0)
    {
      fprintf(stderr, "Failed to allocate memory: %s\n", strerror(errno));
      exit(-1);
    }
}

struct label {
  char * str;
};

const char * skip_to_colon(const char * line)
{
  while(line[0] != ':')
    {
      assert(line[0] != '\0');
      line++;
    }
  return line;
}

const char * assert_token(const char * line, char * buf, char * bitr, const char * token)
{
  while(line[0] != ' ')
    {
      assert(line[0] != '\0');
      bitr[0] = line[0];
      line++;
      bitr++;
    }
  bitr[0] = '\0';

  assert(strcmp(buf, token) == 0);

  return line;
}

const char * read_token(const char * line, char * bitr)
{
  while(line[0] == ' ')
    {
      assert(line[0] != '\0');
      line++;
    }
  while(line[0] != ' ')
    {
      assert(line[0] != '\0');
      bitr[0] = line[0];
      bitr++;
      line++;
    }
  bitr[0] = '\0';
  return line;
}

PCFOP * read_label(const char * line, struct PCFState * st, uint32_t iptr)
{
  ENTRY * newent, * r;
  char buf[LINE_MAX], *bitr;
  PCFOP * ret = (PCFOP *)malloc(sizeof(struct PCFOP));
  check_alloc(ret);

  ret->op = nop;
  ret->type = LABEL_OP;

  bitr = buf;
  
  line = skip_to_colon(line);

  // Skip over the ':'
  line++;

  line = assert_token(line, buf, bitr, "STR");

  while((line[0] == ' ') || (line[0] == '"'))
    {
      assert(line[0] != '\0');
      line++;
    }

  bitr = buf;

  while((line[0] != ' ') && (line[0] != '"') && (line[0] != ')'))
    {
      assert(line[0] != '\0');
      bitr[0] = line[0];
      line++;
      bitr++;
    }
  bitr[0] = '\0';

  newent = (ENTRY*)malloc(sizeof(ENTRY));
  check_alloc(newent);

  newent->key = (char*)malloc(strlen(buf)+1);
  check_alloc(newent->key);
  strcpy(newent->key, buf);

  newent->data = (char*)malloc(sizeof(uint32_t));
  check_alloc(newent->data);
  *((uint32_t*)newent->data) = iptr;

  //hsearch_r(*newent, ENTER, &r, st->labels);
#ifdef __APPLE__
  if((r=hsearch(*newent, ENTER)) == 0)
#else
  if(hsearch_r(*newent, ENTER, &r, st->labels) == 0)
#endif
    {
      fprintf(stderr, "Problem inserting hash table for %s %d: %s\n", newent->key, *((uint32_t*)newent->data), strerror(errno));
      abort();
    }

  return ret;
}

PCFOP * read_initbase(const char * line)
{
  char buf[LINE_MAX], *bitr;
  PCFOP * ret = (PCFOP *)malloc(sizeof(struct PCFOP));
  check_alloc(ret);

  ret->op = initbase_op;
  ret->type = INITBASE_OP;

  uint32_t base;
  bitr = buf;
  bitr[0] = '\0';
  line = skip_to_colon(line);
  line++;

  line = assert_token(line, buf, bitr, "BASE");

  sscanf(line, "%d\n", &base);
  assert(base >= 0);

  ret->data = malloc(sizeof(uint32_t));
  check_alloc(ret->data);

  *((uint32_t*)ret->data) = base;

  return ret;
}

PCFOP * read_gate(const char * line)
{
  char buf[LINE_MAX], *bitr;
  PCFOP * ret = (PCFOP *)malloc(sizeof(PCFOP));
  struct PCFGate * data = (PCFGate *) malloc(sizeof(struct PCFGate));
  uint32_t i = 0;
  check_alloc(ret);
  check_alloc(data);
  bitr = buf;
  ret->op = gate_op;
  ret->type = GATE_OP;
  ret->data = data;

  data->tag = TAG_INTERNAL;

  bitr[0] = '\0';
  for(i = 0; i < 4; i++)
    {
      line = skip_to_colon(line);
      line++;

      bitr = buf;
      line = read_token(line, bitr);
      if(strcmp(buf, "DEST") == 0)
        {
          bitr = buf;
          line = read_token(line, bitr);
          assert(sscanf(buf, "%d", &data->reswire) == 1);
        }
      else if(strcmp(buf, "OP1") == 0)
        {
          bitr = buf;
          line = read_token(line, bitr);
          assert(sscanf(buf, "%d", &data->wire1) == 1);
        }

      else if(strcmp(buf, "OP2") == 0)
        {
          bitr = buf;
          line = read_token(line, bitr);
          assert(sscanf(buf, "%d", &data->wire2) == 1);
        }

      else if(strcmp(buf, "TRUTH-TABLE") == 0)
        {
          bitr = buf;
          line = read_token(line, bitr);
          assert(buf[0] == '#');
          assert(buf[1] == '*');
          data->truth_table = 
            (buf[2] == '1' ? 1 : 0) |
            (buf[3] == '1' ? 2 : 0) |
            (buf[4] == '1' ? 4 : 0) |
            (buf[5] == '1' ? 8 : 0);
        }
      else assert(0);
    }
  return ret;
}

PCFOP * read_copy(const char * line)
{
  char buf[LINE_MAX], *bitr;
  PCFOP * ret =(PCFOP *) malloc(sizeof(PCFOP));
  struct copy_op_data * data = (copy_op_data *)malloc(sizeof(struct copy_op_data));
  uint32_t i = 0;
  check_alloc(ret);
  check_alloc(data);
  bitr = buf;
  ret->op = copy_op;
  ret->type = COPY_OP;
  ret->data = data;
  
  bitr[0] = '\0';
  for(i = 0; i < 3; i++)
    {
      line = skip_to_colon(line);
      line++;

      bitr = buf;
      line = read_token(line, bitr);
      if(strcmp(buf, "DEST") == 0)
        {
          bitr = buf;
          line = read_token(line, bitr);
          assert(sscanf(buf, "%d", &data->dest) == 1);
        }
      else if(strcmp(buf, "OP1") == 0)
        {
          bitr = buf;
          line = read_token(line, bitr);
          assert(sscanf(buf, "%d", &data->source) == 1);
        }

      else if(strcmp(buf, "OP2") == 0)
        {
          bitr = buf;
          line = read_token(line, bitr);
          assert(sscanf(buf, "%d", &data->width) == 1);
        }
      else assert(0);
    }
  return ret;
}

PCFOP * read_arith(const char * line, void (*op)(struct PCFState *, struct PCFOP *))
{
  char buf[LINE_MAX], *bitr;
  PCFOP * ret = (PCFOP *)malloc(sizeof(PCFOP));
  struct arith_op_data * data = (arith_op_data *)malloc(sizeof(struct arith_op_data));
  uint32_t i = 0;
  check_alloc(ret);
  check_alloc(data);
  bitr = buf;
  ret->op = op;
  ret->data = data;

  bitr[0] = '\0';
  for(i = 0; i < 3; i++)
    {
      line = skip_to_colon(line);
      line++;

      bitr = buf;
      line = read_token(line, bitr);
      if(strcmp(buf, "DEST") == 0)
        {
          bitr = buf;
          line = read_token(line, bitr);
          assert(sscanf(buf, "%d", &data->dest) == 1);
        }
      else if(strcmp(buf, "OP1") == 0)
        {
          bitr = buf;
          line = read_token(line, bitr);
          assert(sscanf(buf, "%d", &data->op1) == 1);
        }

      else if(strcmp(buf, "OP2") == 0)
        {
          bitr = buf;
          line = read_token(line, bitr);
          assert(sscanf(buf, "%d", &data->op2) == 1);
        }
      else assert(0);
    }

  return ret;
}

PCFOP * read_add(const char * line)
{
  return read_arith(line, add_op);
}

PCFOP * read_mul(const char * line)
{
  return read_arith(line, mul_op);
}

PCFOP * read_copy_indir(const char * line)
{
  char buf[LINE_MAX], *bitr;
  PCFOP * ret = ( PCFOP *) malloc(sizeof(PCFOP));
  struct copy_op_data * data = (copy_op_data *)malloc(sizeof(struct copy_op_data));
  uint32_t i = 0;
  check_alloc(ret);
  check_alloc(data);
  bitr = buf;
  ret->op = copy_indir_op;
  ret->type = COPY_INDIR_OP;
  ret->data = data;

  bitr[0] = '\0';
  for(i = 0; i < 3; i++)
    {
      line = skip_to_colon(line);
      line++;

      bitr = buf;
      line = read_token(line, bitr);
      if(strcmp(buf, "DEST") == 0)
        {
          bitr = buf;
          line = read_token(line, bitr);
          assert(sscanf(buf, "%d", &data->dest) == 1);
        }
      else if(strcmp(buf, "OP1") == 0)
        {
          bitr = buf;
          line = read_token(line, bitr);
          assert(sscanf(buf, "%d", &data->source) == 1);
        }

      else if(strcmp(buf, "OP2") == 0)
        {
          bitr = buf;
          line = read_token(line, bitr);
          assert(sscanf(buf, "%d", &data->width) == 1);
        }
      else assert(0);
    }
  return ret;
}


PCFOP * read_indir_copy(const char * line)
{
  char buf[LINE_MAX], *bitr;
  PCFOP * ret = (PCFOP *)malloc(sizeof(PCFOP));
  struct copy_op_data * data = (copy_op_data *)malloc(sizeof(struct copy_op_data));
  uint32_t i = 0;
  check_alloc(ret);
  check_alloc(data);
  bitr = buf;
  ret->op = indir_copy_op;
  ret->type = INDIR_COPY_OP;
  ret->data = data;

  bitr[0] = '\0';
  for(i = 0; i < 3; i++)
    {
      line = skip_to_colon(line);
      line++;

      bitr = buf;
      line = read_token(line, bitr);
      if(strcmp(buf, "DEST") == 0)
        {
          bitr = buf;
          line = read_token(line, bitr);
          assert(sscanf(buf, "%d", &data->dest) == 1);
        }
      else if(strcmp(buf, "OP1") == 0)
        {
          bitr = buf;
          line = read_token(line, bitr);
          assert(sscanf(buf, "%d", &data->source) == 1);
        }

      else if(strcmp(buf, "OP2") == 0)
        {
          bitr = buf;
          line = read_token(line, bitr);
          assert(sscanf(buf, "%d", &data->width) == 1);
        }
      else assert(0);
    }
  return ret;
}

PCFOP * read_const(const char * line)
{
  char buf[LINE_MAX], *bitr;
  PCFOP * ret = ( PCFOP *) malloc(sizeof(PCFOP));
  struct const_op_data * data = ( const_op_data *)malloc(sizeof(struct const_op_data));

  check_alloc(ret);
  check_alloc(data);


  bitr = buf;

  ret->op = const_op;
  ret->type = CONST_OP;
  ret->data = data;

  bitr[0] = '\0';
  line = skip_to_colon(line);
  line++;

  bitr = buf;
  line = read_token(line, bitr);
  if(strcmp(buf, "DEST") == 0)
    {
      bitr = buf;
      line = read_token(line, bitr);
      assert(sscanf(buf, "%d", &data->dest) == 1);
    }
  else if(strcmp(buf, "OP1") == 0)
    {
      bitr = buf;
      line = read_token(line, bitr);
      assert(sscanf(buf, "%d", &data->value) == 1);
    }
  else assert(0);

  line = skip_to_colon(line);
  line++;

  bitr = buf;
  line = read_token(line, bitr);
  if(strcmp(buf, "DEST") == 0)
    {
      bitr = buf;
      line = read_token(line, bitr);
      assert(sscanf(buf, "%d", &data->dest) == 1);
    }
  else if(strcmp(buf, "OP1") == 0)
    {
      bitr = buf;
      line = read_token(line, bitr);
      assert(sscanf(buf, "%d", &data->value) == 1);

    }
  else assert(0);

  return ret;
}

int count_tokens_to_close_paren(const char * line)
{
  uint32_t cnt = 0;
  while(line[0] == ' ') line++;
  assert(line[0] == '(');
  line++;
  while(line[0] != ')')
    {
      assert(line[0] != '\0');
      cnt++;
      while((line[0] != ' ') && (line[0] != ')')) line++;
      while((line[0] == ' ') && (line[0] != ')')) line++;
    }
  return cnt;
}

PCFOP * read_bits(const char * line)
{
  char buf[LINE_MAX], *bitr;
  PCFOP * ret = (PCFOP*)malloc(sizeof(PCFOP));
  uint32_t i = 0;

  struct bits_op_data * data = (struct bits_op_data*)malloc(sizeof(struct bits_op_data));

  check_alloc(ret);
  check_alloc(data);

  ret->op = bits_op;
  ret->type = BITS_OP;
  ret->data = data;
  bitr = buf;

  for(i = 0; i < 2; i++)
    {
      line = skip_to_colon(line);
      line++;

      bitr = buf;
      line = read_token(line, bitr);

      if(strcmp(buf, "DEST") == 0)
        {
          int cnt = count_tokens_to_close_paren(line);
          data->ndests = cnt;
          char buf2[10];
          int i, j, k;

          data->dests = (uint32_t*)malloc(cnt * sizeof(uint32_t));
          check_alloc(data->dests);
          i = 0;
          j = 0;
          k = 0;
          while(line[i] == '(') i++;
          while(line[i] == ' ') i++;
          while(line[i] == '(') i++;
          while(line[i] != ')')
            {
              assert(line[i] != '\0');
              while((line[i] != ' ') && (line[i] != ')'))
                {
                  assert(line[i] != '\0');
                  buf2[j] = line[i];
                  i++;
                  j++;
                }
              buf2[j] = '\0';
              assert(sscanf(buf2, "%d", &data->dests[k]) == 1);
              k++;
              while(line[i] == ' ') 
                {
                  assert(line[i] != '\0');
                  i++;
                }
              j = 0;
            }
          assert(k == cnt);
          line += i;
        }
      else if(strcmp(buf, "OP1") == 0)
        {
          bitr = buf;
          line = read_token(line, bitr);
          assert(sscanf(buf, "%d", &data->source) == 1);
        }
      else
        {
          assert(0);
        }
    }
  return ret;
}

PCFOP * read_join(const char * line)
{
  char buf[LINE_MAX], *bitr;
  PCFOP * ret = (PCFOP*)malloc(sizeof(PCFOP));
  uint32_t i = 0;

  struct join_op_data * data = (struct join_op_data*)malloc(sizeof(struct join_op_data));

  check_alloc(ret);
  check_alloc(data);

  ret->op = join_op;
  ret->type = JOIN_OP;
  ret->data = data;
  bitr = buf;

  for(i = 0; i < 2; i++)
    {
      line = skip_to_colon(line);
      line++;

      bitr = buf;
      line = read_token(line, bitr);

      if(strcmp(buf, "OP1") == 0)
        {
          int cnt = count_tokens_to_close_paren(line);
          data->nsources = cnt;
          char buf2[10];
          int i, j, k;

          data->sources = (uint32_t*)malloc(cnt * sizeof(uint32_t));
          check_alloc(data->sources);
          i = 0;
          j = 0;
          k = 0;
          while(line[i] == '(') i++;
          while(line[i] == ' ') i++;
          while(line[i] == '(') i++;
          while(line[i] != ')')
            {
              assert(line[i] != '\0');
              while((line[i] != ' ') && (line[i] != ')'))
                {
                  assert(line[i] != '\0');
                  buf2[j] = line[i];
                  i++;
                  j++;
                }
              buf2[j] = '\0';
              assert(sscanf(buf2, "%d", &data->sources[k]) == 1);
              k++;
              while(line[i] == ' ') 
                {
                  assert(line[i] != '\0');
                  i++;
                }
              j = 0;
            }
          assert(k == cnt);
          line += i;
        }
      else if(strcmp(buf, "DEST") == 0)
        {
          bitr = buf;
          line = read_token(line, bitr);
          assert(sscanf(buf, "%d", &data->dest) == 1);
        }
      else
        {
          assert(0);
        }
    }
  return ret;
}


PCFOP * read_mkptr(const char * line)
{
  char buf[LINE_MAX], *bitr;
  PCFOP * ret = ( PCFOP * )malloc(sizeof(PCFOP));
  uint32_t * data = (uint32_t * )malloc(sizeof(uint32_t));

  check_alloc(ret);
  check_alloc(data);

  ret->op = mkptr_op;
  ret->type = MKPTR_OP;
  ret->data = data;

  bitr = buf;

  line = skip_to_colon(line);
  line++;

  bitr = buf;
  line = read_token(line, bitr);

  if(strcmp(buf, "DEST") == 0)
    {
      bitr = buf;
      line = read_token(line, bitr);
      assert(sscanf(buf, "%d", data) == 1);
    }
  else
    assert(0);

  return ret;
}

PCFOP * read_call(const char * line)
{
  char buf[LINE_MAX], *bitr;
  PCFOP * ret = (PCFOP*)malloc(sizeof(PCFOP));
  uint32_t i = 0;

  struct call_op_data * data = (struct call_op_data*)malloc(sizeof(struct call_op_data));

  check_alloc(ret);
  check_alloc(data);

  ret->op = call_op;
  ret->type = CALL_OP;
  ret->data = data;
  bitr = buf;

  for(i = 0; i < 2; i++)
    {
      line = skip_to_colon(line);
      line++;

      bitr = buf;
      line = read_token(line, bitr);

      if(strcmp(buf, "NEWBASE") == 0)
        {
          bitr = buf;
          line = read_token(line, bitr);
          assert(sscanf(buf, "%d", &data->newbase) == 1);
        }
      else if(strcmp(buf, "FNAME") == 0)
        {
          ENTRY * ent = (ENTRY*)malloc(sizeof(ENTRY));
          check_alloc(ent);

          while((line[0] == ' ') || (line[0] == '"'))
            {
              assert(line[0] != '\0');
              line++;
            }

          bitr = buf;
          line = read_token(line, bitr);

          assert(buf[strlen(buf)-1] == '"');
          buf[strlen(buf)-1] = '\0';

          ent->key = (char*)malloc(strlen(buf)+1);
          strcpy(ent->key, buf);
          data->target = ent;
        }
      else
        {
          assert(0);
        }
    }
  
  return ret;
}

PCFOP * read_branch(const char * line)
{
  char buf[LINE_MAX], *bitr;
  PCFOP * ret = (PCFOP*)malloc(sizeof(PCFOP));
  uint32_t i = 0;

  struct branch_op_data * data = (struct branch_op_data*)malloc(sizeof(struct call_op_data));

  check_alloc(ret);
  check_alloc(data);

  ret->op = branch_op;
  ret->type = BRANCH_OP;
  ret->data = data;
  bitr = buf;

  for(i = 0; i < 2; i++)
    {
      line = skip_to_colon(line);
      line++;

      bitr = buf;
      line = read_token(line, bitr);

      if(strcmp(buf, "CND") == 0)
        {
          bitr = buf;
          line = read_token(line, bitr);
          assert(sscanf(buf, "%d", &data->cnd_wire) == 1);
        }
      else if(strcmp(buf, "TARG") == 0)
        {
          ENTRY * ent = (ENTRY*)malloc(sizeof(ENTRY));
          check_alloc(ent);

          while((line[0] == ' ') || (line[0] == '"'))
            {
              assert(line[0] != '\0');
              line++;
            }

          bitr = buf;
          line = read_token(line, bitr);

          assert(buf[strlen(buf)-1] == '"');
          buf[strlen(buf)-1] = '\0';

          ent->key = (char*)malloc(strlen(buf)+1);
          strcpy(ent->key, buf);
          data->target = ent;
        }
      else
        {
          assert(0);
        }
    }
  
  return ret;
}

PCFOP * read_ret(const char * line)
{
  PCFOP * ret = (PCFOP*)malloc(sizeof(PCFOP));
  ret->op = ret_op;
  ret->type = RET_OP;
  return ret;
}

PCFOP * read_clear(const char * line)
{
  char buf[LINE_MAX], *bitr;
  PCFOP * ret = (PCFOP*)malloc(sizeof(PCFOP));
  struct clear_op_data * data = (struct clear_op_data*)malloc(sizeof(struct clear_op_data));

  check_alloc(ret);
  check_alloc(data);

  ret->op = clear_op;
  ret->type = CLEAR_OP;

  bitr = buf;
  bitr[0] = '\0';
  line = skip_to_colon(line);
  line++;

  line = assert_token(line, buf, bitr, "LOCALSIZE");

  sscanf(line, "%d\n", &data->localsize);
  assert(data->localsize >= 0);

  ret->data = data;
  return ret;
}

PCFOP * read_instr(struct PCFState * st, const char * line, uint32_t iptr)
{
  char buf[LINE_MAX], *bitr;
  buf[0] = '\0';
  bitr = buf;

  assert(line[0] == '(');
  line++;

  // step through the line, copying the parts that we actually want
  while((line[0] != ' ') && (line[0] != ')'))
    {
      bitr[0] = line[0];
      line++;
      bitr++;
    }
  bitr[0] = '\0';

  if(strcmp(buf, "LABEL") == 0)
    return read_label(line, st, iptr);
  else if(strcmp(buf, "INITBASE") == 0)
    return read_initbase(line);
  else if(strcmp(buf, "CONST") == 0)
    return read_const(line);
  else if(strcmp(buf, "GATE") == 0)
    return read_gate(line);
  else if(strcmp(buf, "BITS") == 0)
    return read_bits(line);
  else if(strcmp(buf, "MKPTR") == 0)
    return read_mkptr(line);
  else if(strcmp(buf, "COPY") == 0)
    return read_copy(line);
  else if(strcmp(buf, "COPY-INDIR") == 0)
    return read_copy_indir(line);
  else if(strcmp(buf, "INDIR-COPY") == 0)
    return read_indir_copy(line);
  else if(strcmp(buf, "CALL") == 0)
    return read_call(line);
  else if(strcmp(buf, "RET") == 0)
    return read_ret(line);
  else if(strcmp(buf, "BRANCH") == 0)
    return read_branch(line);
  else if(strcmp(buf, "CLEAR") == 0)
    return read_clear(line);
  else if(strcmp(buf, "JOIN") == 0)
    return read_join(line);
  else if(strcmp(buf, "ADD") == 0)
    return read_add(line);
  else if(strcmp(buf, "MUL") == 0)
    return read_mul(line);
  assert(0);
}


/**
   this function returns a PCFState object
   it accepts a filename, two keys, and a function used to copy keys
   when this funciton is called, the two keys 
   (for some reason yet to be discerned) are offset by 1
*/ 

//PCFState * load_pcf_file(const char * fname, void * key0, void * key1, void *(*copy_key)(void*))
PCFState * load_pcf_file(const char * fname, void * key0, void * key1, void (*copy_key)(void*,void*))
{
  FILE * input;
  PCFState * ret; 
  char line[LINE_MAX];
  uint32_t icount = 0;
  uint32_t i = 0;

  uint32_t num_wires = NUM_WIRES;

  ret = (PCFState*)malloc(sizeof(struct PCFState));
  check_alloc(ret);

  //ret->alice_outputs = 0;
  //ret->bob_outputs = 0;
  ret->inp_i = 0;
  ret->inp_idx = 0; // should not be strictly necessary because will be set before first use

  ret->accum=0.0;
  ret->accum2=0.0;

  fprintf(stdout,"allocate const keys:\n");
  
  ret->constant_keys[0] =(void*) malloc(4*sizeof(uint32_t));
  ret->constant_keys[1] =(void*) malloc(4*sizeof(uint32_t));
  fprintf(stdout,"set up const keys:\n");

  copy_key(key0, ret->constant_keys[0]);
  fprintf(stdout,"set const key 0\n");
  copy_key(key1, ret->constant_keys[1]);
  fprintf(stdout,"set const key 1\n");
  ret->copy_key = copy_key;
  fprintf(stdout,"set copy key function\n");
  ret->call_stack = 0;
  ret->done = 0;
#ifdef __APPLE__
// TODO:fix this ????
#else
  ret->labels = (struct hsearch_data *)malloc(sizeof(struct hsearch_data));
  check_alloc(ret->labels);
  memset(ret->labels, 0, sizeof(struct hsearch_data));
#endif

  ret->wires = (struct wire *)malloc(num_wires * sizeof(struct wire)); // !note here the limit on size of the wire table
  check_alloc(ret->wires);
  
  fprintf(stdout,"allocate keys\n");
  for(i = 0; i < num_wires; i++)
    {
      //fprintf(stdout,"%i",i);
      ret->wires[i].flags = KNOWN_WIRE;
      ret->wires[i].value = 0;
      ret->wires[i].keydata = (void *)malloc(4*sizeof(uint32_t));
      check_alloc(ret->wires[i].keydata);
      
      //      copy_key(key0,ret->wires[i].keydata);
      //ret->wires[i].keydata = copy_key(key0);
    }
  fprintf(stdout,"keys allocated\n");
  

  
  ret->done = 0;
  ret->base = 1;
  ret->PC = 0;

  fprintf(stderr, "%s\n", fname);
  input = fopen(fname, "r");
  if(input == 0)
    {
      fprintf(stderr, "%s: %s\n", fname, strerror(errno));
      assert(0);
    }

  while(!feof(input))
    {
      fgets(line, LINE_MAX-1, input);
      icount++;
    }

#ifdef __APPLE__
  if(hcreate(icount) == 0)
#else
  if(hcreate_r(icount, ret->labels) == 0)
#endif
    {
      fprintf(stderr, "Unable to allocate hash table: %s\n", strerror(errno));
      abort();
      //      exit(-1);
    }

  ret->icount = icount;
  ret->ops = (PCFOP*)malloc(icount * sizeof(PCFOP));
  check_alloc(ret->ops);

  assert(fseek(input, 0, SEEK_SET) == 0);

  icount = 0;


  // read all of the ops into memory
  while(!feof(input))
    {
      PCFOP * op;
      fgets(line, LINE_MAX-1, input);
      op = read_instr(ret, line, icount);
      op->idx = icount;

      ret->ops[icount] = *op;
      
      //free(op);
      // is this a good idea after setting the pointer to it?
      // read_instr uses malloc, and then ret sets the shallow copy
      // but I really think free is a mistake
      
      icount++;
    }

  fclose(input);

  ret->wires[0].value = 1;

  ret->copy_key(ret->constant_keys[1],ret->wires[0].keydata);
  //  ret->wires[0].keydata = ret->copy_key(ret->constant_keys[1]);
  ret->wires[0].flags = KNOWN_WIRE;

  
  // at this point, we can construct a graph of the circuit using dependency analysis
  


  return ret;
}

void finalize(PCFState * st)
{

  //  fprintf(stderr, "finalize\n");
  
  // uint32_t i = 0;
  
  // still need some way to delete the keys and free memory
  //for(i = 0; i < 200000; i++)
  //  {
  //    if(st->wires[i].keydata != 0)
  //      st->delete_key(st->wires[i].keydata);
  //  } 
  uint32_t i = 0;
  for(i=0;i< NUM_WIRES;i++){
    free(st->wires[i].keydata);
  }
  //  free(st->wires);
  free(st);
  
  //  fprintf(stderr,"done finalize\n");
  
  fprintf(stderr,"accumulator: %f\n",st->accum);
  fprintf(stderr,"accumulator2: %f\n",st->accum2);

}

void apply_flow(struct PCFState *st, struct PCFOP * op, uint32_t * table){
  switch (op->type){
  case GATE_OP:
    gate_flow(st,op,table);
    break;
  case BITS_OP:
    bits_flow(st,op,table);
    break;
  case CONST_OP:
    const_flow(st,op,table);
    break;
  case ADD_OP:
  case SUB_OP:
  case MUL_OP:
    arith_flow(st,op,table);
    break;
  case INITBASE_OP:
    initbase_flow(st,op,table);
    break;
  case CLEAR_OP:
    clear_flow(st,op,table);
    break;
  case MKPTR_OP:
    mkprt_flow(st,op,table);
    break;
  case COPY_INDIR_OP:
    copy_indir_flow(st,op,table);
    break;
  case INDIR_COPY_OP:
    indir_copy_flow(st,op,table);
    break;
  case CALL_OP:
    call_flow(st,op,table);
    break;
  case RET_OP:
    ret_flow(st,op,table);
    break;
  case BRANCH_OP:
    branch_flow(st,op,table);
    break;
  case LABEL_OP:
    label_flow(st,op,table);
    break;
  case JOIN_OP:
    join_flow(st,op,table);
    break;
  default:
    fprintf(stdout,"error determining op type!");
    break;
  }
  
}


PCFState * build_tree(struct PCFState *st){
  uint32_t i;

  // table good for the whole circuit, in the worst case
  // must be cleared for each function though
  uint32_t * p = (uint32_t *)malloc(sizeof(uint32_t*)*st->icount);
  for(uint32_t j = 0; j < st->icount; j++){
    p[j]=0;
  }

  // run through the list of ops
  // building dependencies
  for(i=0; i< st->icount; i++){
    apply_flow(st,&st->ops[i],p);
  }

  return st;
}


void evaluate_circuit(struct PCFState *st){

  // will need here a couple of threads
  // because no internet, use a threading abstraction for now
  
  // execute circuit on main thread,
  // dispatching gates in parallel
  

#ifndef __APPLE__
  clock_gettime(CLOCK_REALTIME, &(st->requestStart)); 
#endif

  uint32_t i = 0;
  i++;
  while(st->done != 0)
    {
      if(st->curgate ==0){
        //dispatch
      } else{
        // take the top off of the ready queue
        // execute it
        // add successors to ready queue, if ready
        // continue

        // in fact, I think this can be evaluated in parallel as well

        /*
         PCFOP op = readyQueue.pop();
         // this pop should be atomic and blocking
         // so only one thread at a time can get it
         
         if(op != 0){
           //
           

         } else{
           // nothing on the ready queue
           // either an error or we just need to wait
           // figure that out later
         }
        */

      }
      // note that first I am implementing the parallel garbler,
      // then I can make adaptations to communicate properly with the evaluator
      // but want to do first construct the parallel infrastructure to run in sequence
    }



#ifndef __APPLE__
  clock_gettime(CLOCK_REALTIME, &(st->requestEnd));
  st->accum += ( st->requestEnd.tv_sec - st->requestStart.tv_sec )
    + ( st->requestEnd.tv_nsec - st->requestStart.tv_nsec )
    / BILLION;
#endif

}

// this function will change for parallel garbler
// and really the whole execution paradigm
struct PCFGate * get_next_gate(struct PCFState * st)
{
  
  //fprintf(stderr, "get next gate");
  //  std::cout << "get next gate" << std::endl;
#ifndef __APPLE__
  clock_gettime(CLOCK_REALTIME, &(st->requestStart)); 
#endif

  st->curgate = 0;
  //fprintf(stderr,"program counter: %u\n",st->PC);
  while((st->curgate == 0) && (st->done == 0))
  {
    // if curgate is 0, why are we executing things?
    // that check is an indicator (a null pointer)
    // that we are on a non-gate instruction


      // call the instruction's op function
      st->ops[st->PC].op(st, &st->ops[st->PC]);
      // then increment the program counter
      st->PC++;
      // and verify that the PC has not eclipsed the instruction count
      //assert((st->PC < st->icount));
    }

  //  if((st->curgate == 0) || (st->done != 0))
  // this seems redundant, since st->curgate should only be 0 after the above loop if st->done is nonzero
  // might want to look this up in an old version of PCF - June 3 2016
 
#ifndef __APPLE__
  clock_gettime(CLOCK_REALTIME, &(st->requestEnd));
  st->accum += ( st->requestEnd.tv_sec - st->requestStart.tv_sec )
    + ( st->requestEnd.tv_nsec - st->requestStart.tv_nsec )
    / BILLION;
#endif

  if(st->done != 0)
  {
      // no more gates and the PCFState thinks it's done
    //      fprintf(stderr,"enter finalize\n");
    finalize(st);
    return 0;
  }
  else {
    // return the next gate
    return st->curgate;
  }
}


void * get_wire_key(struct PCFState * st, uint32_t idx)
{
  return st->wires[idx].keydata;
}
/*
void set_wire_key(struct PCFState * st, uint32_t idx, void * kd)
{
  setWireKey(&st->wires[idx], kd);
}
*/

/*
wire * getWire(struct PCFState * st, uint32_t idx)
{
  return &st->wires[idx];
}
*/

void set_external_circuit(struct PCFState * st, void * ec)
{
  st->external_circuit = ec;
}

void * get_external_circuit(struct PCFState * st)
{
  return st->external_circuit;
}

void set_callback(struct PCFState * st, void* (*callback)(struct PCFState *,struct PCFGate*))
{
  st->callback = callback;
}

//void set_key_copy_function(struct PCFState * st, void *(*f)(void*))
void set_key_copy_function(struct PCFState * st, void (*f)(void*,void*))
{
  st->copy_key = f;
}

void set_key_delete_function(struct PCFState * st, void (*f)(void*))
{
  st->delete_key = f;
}

uint32_t read_alice_length(const char * fname)
{
  uint32_t alice_in = -1;
  char line[LINE_MAX];
  FILE * cnf;
  cnf = fopen(fname, "r");
  // Assume that the inputs are all on one line -- probably a bad assumption, but whatever
  fgets(line, LINE_MAX-1, cnf);
  fclose(cnf);
  alice_in = (strlen(line) - 1) * 4;
  fprintf(stderr,"alice length: %u\n",alice_in);
  return alice_in;
}

uint32_t read_bob_length(const char * fname)
{
  uint32_t bob_in = 32;
  FILE * cnf;
  char line[LINE_MAX];
  cnf = fopen(fname, "r");
  // Bob's input is on the second line
  //
  // The input file format should be like this:
  //
  // 0xALICEINPUTSINHEX
  // 0x0000000000000000
  //
  // or
  //
  // 0x0000000000000000
  // 0xBOBINPUTSINHEX00
  fgets(line, LINE_MAX-1, cnf);
  fgets(line, LINE_MAX-1, cnf);
  fclose(cnf);
  bob_in = (strlen(line) - 1) * 4;
  fprintf(stderr,"bob length: %u\n",bob_in);
  return bob_in;
}

char* get_alice_input(uint32_t length, const char* fname){
  FILE * cnf;
  char * alice_in = (char*) malloc(sizeof(char)*length);
  cnf = fopen(fname,"r");
  // alice's input is on the first line.
  // read length characters into alice_in
  fgets(alice_in, length, cnf);
  fclose(cnf);
  return alice_in;
}

char* get_bob_input(uint32_t length, const char* fname){
  FILE * cnf;
  char * bob_in = (char*)malloc(sizeof(char)*length);
  char discard[LINE_MAX];
  cnf = fopen(fname,"r");
  // bob's input is on the second line
  // read length characters into bob_in
  fgets(discard, LINE_MAX-1,cnf);
  fgets(bob_in, length, cnf);
  fclose(cnf);
  return bob_in;
}
