#include "opgen.h"



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


  //ret->num_succs = 1;
  //ret->num_preds = 2; // is it always 2?

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


  //ret->num_succs = 1;
  //ret->num_preds = 1;

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


  //ret->num_preds = 2;
  //ret->num_succs =1;

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

  // ret-> num_preds = &data->width;

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

