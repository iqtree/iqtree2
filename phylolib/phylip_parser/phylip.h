#ifndef PHYLIP_H
#define PHYLIP_H

#define PHYLIP_SEQUENTIAL       0x00000001
#define PHYLIP_INTERLEAVED      0x00000002

#define DNA_DATA                0x00000001
#define PROT_DATA               0x00000002

struct phylip_data
 {
   int          taxa;
   int          seqlen;
   char      ** label;
   char      ** seq;
 };

struct phylip_data * alloc_phylip_struct (int taxa, int seqlen);
struct phylip_data * pl_phylip_parse (const char *, int);
void pl_phylip_subst (struct phylip_data *, int);
void free_phylip_struct (struct phylip_data * pd);

#endif
