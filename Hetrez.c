/*******************************************************************************************
 *
 *  Hetrez: Analyze and compare het-mers between two k-mer tables to see if they
 *    represent the same sample, different sample but same species, or are from
 *    two different species.
 *
 *  Author:  Gene Myers
 *  Date  :  May, 2021
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>
#include <pthread.h>

#undef  DEBUG_PARTITION
#define DEBUG_THREADS
#define REPORT

#include "libfastk.h"

static char *Usage = " [-e<int>:<int>] <source1>[.ktab] <source2>[.ktab]";

static int   ETHRESH_1, ETHRESH_2;   //  Count threshold for reliable k-mers
static int   VERBOSE;
static int   NTHREADS;

/****************************************************************************************
 *
 *  Print & compare utilities
 *
 *****************************************************************************************/

#define  COUNT_PTR(p)   ((uint16 *) (p+kbyte))

static char dna[4] = { 'a', 'c', 'g', 't' };

static char *fmer[256], _fmer[1280];

static void setup_fmer_table()
{ char *t;
  int   i, l3, l2, l1, l0;

  i = 0;
  t = _fmer;
  for (l3 = 0; l3 < 4; l3++)
   for (l2 = 0; l2 < 4; l2++)
    for (l1 = 0; l1 < 4; l1++)
     for (l0 = 0; l0 < 4; l0++)
       { fmer[i] = t;
         *t++ = dna[l3];
         *t++ = dna[l2];
         *t++ = dna[l1];
         *t++ = dna[l0];
         *t++ = 0;
         i += 1;
       }
}

static void print_seq(uint8 *seq, int len)
{ int i, b, k;

  b = len >> 2;
  for (i = 0; i < b; i++)
    printf("%s",fmer[seq[i]]);
  k = 6;
  for (i = b << 2; i < len; i++)
    { printf("%c",dna[(seq[b] >> k) & 0x3]);
      k -= 2;
    }
}

static void print_hap(uint8 *seq, int len, int half)
{ int i, b, h, k;

  h = half >> 2;
  b = len >> 2;
  for (i = 0; i < h; i++)
    printf("%s",fmer[seq[i]]);
  k = 6;
  for (i = h << 2; k >= 0; i++)
    { if (i == half)
        printf("%c",dna[(seq[h] >> k) & 0x3]-32);
      else
        printf("%c",dna[(seq[h] >> k) & 0x3]);
      k -= 2;
    }
  for (i = h+1; i < b; i++)
    printf("%s",fmer[seq[i]]);
  k = 6;
  for (i = b << 2; i < len; i++)
    { printf("%c",dna[(seq[b] >> k) & 0x3]);
      k -= 2;
    }
}

static inline int mycmp(uint8 *a, uint8 *b, int n)
{ while (n-- > 0)
    { if (*a++ != *b++)
        return (a[-1] < b[-1] ? -1 : 1);
    }
  return (0);
}

static inline int mypref(uint8 *a, uint8 *b, int n)
{ int   i;
  uint8 x, y;
 
  for (i = 0; i <= n; i += 4)
    { if (*a != *b)
        { x = *a;
          y = *b;
          if ((x & 0xf0) != (y & 0xf0))
            if ((x & 0xc0) != (y & 0xc0))
              return (i);
            else
              return (i + 1);
          else
            if ((x & 0xfc) != (y & 0xfc))
              return (i + 2);
            else
              return (i + 3);
        }
      a += 1;
      b += 1;
    }
  return (n+1);
}


/****************************************************************************************
 *
 *  Find Haplotype Pairs
 *
 *****************************************************************************************/

typedef struct
  { int          tid;
    Kmer_Stream *T;
    Kmer_Stream *U;
    int64       *begs;
    int64       *ends;
    int64        uhc1, uhc2;
    int64        hpne, hpeq;
  } TP;

static void *analysis_thread(void *args)
{ TP *parm = (TP *) args;
  int           tid   = parm->tid;
  Kmer_Stream  *T     = parm->T;
  Kmer_Stream  *U     = parm->U;
  int64        *begs  = parm->begs;
  int64         Tend  = parm->ends[0];
  int64         Uend  = parm->ends[1];

   uint8  prefs[] = { 0x3f, 0x0f, 0x03, 0x00 };

  uint8 *cache1, *cptr1, *ctop1, *nptr1;
  uint8 *cache2, *cptr2, *ctop2, *nptr2;

  int    e, f, i;
  int    index[8];
  uint8 *finger[8];
  uint8 *flimit[8];

  int    a, advn[8];
  int    c, good[8];
  int    mc, hc;
  uint8 *mr, *hr;

  int64  uhc1, uhc2;
  int64  hpne, hpeq;

  (void) print_seq;
  (void) print_hap;

  setup_fmer_table();

  int kmer  = T->kmer;
  int tbyte = T->tbyte;
  int kbyte = T->kbyte;

  int khalf = kmer/2;
  int mask  = prefs[khalf&0x3];
  int offs  = (khalf >> 2) + 1;
  int rem   = ((kmer+3) >> 2) - offs;

#ifdef DEBUG_PARTITION
  printf("Extension = K[%d]&%02x . K[%d..%d)\n",offs-1,mask,offs,offs+rem);
#endif

  cache1 = Malloc(4097*tbyte,"Allocating entry buffer");
  ctop1  = cache1 + 4096*tbyte;

  cache2 = Malloc(4097*tbyte,"Allocating entry buffer");
  ctop2  = cache2 + 4096*tbyte;

  uhc1  = 0;
  uhc2  = 0;
  hpne  = 0;
  hpeq  = 0;

#ifdef DEBUG_THREADS
  printf("Doing %d:",tid);
  for (c = 0; c < 2; c++)
    printf(" [%lld-%lld]",begs[c],parm->ends[c]);
  printf("\n");
#endif

  if (tid != 0)
    { T = Clone_Kmer_Stream(T);
      U = Clone_Kmer_Stream(U);
    }

  GoTo_Kmer_Index(T,begs[0]);
  GoTo_Kmer_Index(U,begs[1]);

  while (1)
    { if (T->cidx < Tend)
        {  Current_Entry(T,cache1);
           if (U->cidx < Uend)
             { Current_Entry(U,cache2);
               e = mycmp(cache1,cache2,kbyte);
             }
           else
             e = -1;
        }
      else
         { if (U->cidx < Tend)
             { Current_Entry(U,cache2);
               e = 1;
             }
           else
             break;
        }

      if (e <= 0)
        { cptr1 = cache1;
          nptr1 = cache1 + tbyte;
          Next_Kmer_Entry(T);
        }
      else
        { cptr1 = cache2;
          nptr1 = cache1;
        }
      if (e >= 0)
        { cptr2 = cache2;
          nptr2 = cache2 + tbyte;
          Next_Kmer_Entry(U);
        }
      else
        { cptr2 = cache1;
          nptr2 = cache2;
        }

      f = 0;
      index[f++] = 0;
      for ( ; T->cidx < Tend; Next_Kmer_Entry(T))
        { int x = mypref(cptr1,Current_Entry(T,nptr1),khalf);
          if (x < khalf)
            break;
          if (x == khalf)
            index[f++] = nptr1-cache1;
          if (nptr1 >= ctop1)
            { int64 cidx = ctop1-cache1;
              int64 cmax = ((cidx*14)/(10*tbyte) + 2048)*tbyte;
              cache1 = Realloc(cache1,cmax+tbyte,"Reallocting entry buffer");
              ctop1 = cache1 + cmax;
              nptr1 = cache1 + cidx;
            }
          cptr1 = nptr1;
          nptr1 = cptr1+tbyte;
        }

      e = f;
      index[f++] = 0;
      for ( ; U->cidx < Uend; Next_Kmer_Entry(U))
        { int x = mypref(cptr2,Current_Entry(U,nptr2),khalf);
          if (x < khalf)
            break;
          if (x == khalf)
            index[f++] = nptr2-cache2;
          if (nptr2 >= ctop2)
            { int64 cidx = ctop2-cache2;
              int64 cmax = ((cidx*14)/(10*tbyte) + 2048)*tbyte;
              cache2 = Realloc(cache2,cmax+tbyte,"Reallocting entry buffer");
              ctop2 = cache2 + cmax;
              nptr2 = cache2 + cidx;
            }
          cptr2 = nptr2;
          nptr2 = cptr2+tbyte;
        }

#ifdef DEBUG_PARTITION
      printf("part %d:",f);
      for (i = 0; i < e; i++)
        printf(" %d",index[i]/tbyte);
      printf(" /");
      for (i = e; i < f; i++)
        printf(" %d",index[i]/tbyte);
      printf(" == %ld %ld\n",(nptr1-cache1)/tbyte,(nptr2-cache2)/tbyte);
#endif

      if (e <= 1 && f-e <= 1)
        continue;

      for (i = 0; i < f; i++)
        if (i < e)
          finger[i] = cache1 + index[i];
        else
          finger[i] = cache2 + index[i];
      for (i = 1; i < f; i++)
        flimit[i-1] = finger[i];
      flimit[e-1] = nptr1;
      flimit[f-1] = nptr2;

#define ADD(i)                                  \
{ int cn = *COUNT_PTR(finger[i]);               \
  advn[a++] = i;                                \
  if (i < e)                                    \
    { if (ETHRESH_1 <= cn) 			\
        good[c++] = i;                          \
    }                                           \
  else                                          \
    { if (ETHRESH_2 <= cn) 			\
        good[c++] = i;                          \
    }                                           \
}

#define SET(i)  \
{ mc = hc;      \
  mr = hr;      \
  a = c = 0;    \
  ADD(i);       \
}

      while (1)
        { for (i = 0; i < f; i++)
            if (finger[i] < flimit[i])
              break;
          if (i >= f)
            break;
          hr = finger[i]+offs;
          hc = hr[-1] & mask;
          SET(i);
          for (i++; i < f; i++)
            if (finger[i] < flimit[i])
              { hr = finger[i]+offs;
                hc = hr[-1] & mask;
                if (hc == mc)
                  { int v = mycmp(hr,mr,rem);
                    if (v == 0)
                      ADD(i)
                    else if (v < 0)
                      SET(i)
                  }
                else if (hc < mc)
                  SET(i)
              }

#ifdef DEBUG_PARTITION
          { int j;

            j = 0;
            for (i = 0; i < a; i++)
              { printf(" %d",advn[i]);
                if (j < c && advn[i] == good[j])
                  { printf("+");
                    j += 1;
                  }
              }
            printf("\n");
          }
#endif

          if (c > 1)
            { int j, d1, d2;

/*
              j  = 0;
              for (i = 0; i < a; i++)
                { x = advn[i];
                  s = ((finger[x]+offs)>>shift) & 0x3;
                  if (x < e)
                    { table1[s].fing = finger[x]; 
                      s = (1<<s);
                      if (j < c && x == good[j])
                        { table1[s].good = 1;
                          gset |= (1 << s);
                          j  += 1;
                        }
                      aset |= (1 << s);
                    }
                  else
                    { table2[s].fing = finger[x]; 
                      s = (1<<(s+4));
                      if (j < c && x == good[j])
                        { table2[s].good = 1;
                          gset |= (1 << s);
                          j  += 1;
                        }
                    }
                  aset |= (1 << s);
                }

              for (i = 0; i < a; i++)
                { x = advn[i];
                  s = ((finger[x]+offs)>>shift) & 0x3;
                  if (x < e)
                    s = (1<<s);
                  else
                    s = (1<<(s+4));
                  if (j < c && x == good[j])
                    { gset |= (1 << s);
                      j += 1;
                    }
                  aset |= (1 << s);
                }
              gset = complete[gset]
              if ((gset & aset) == gset)
                salvage
*/
                 
              d1 = d2 = 0;
              for (i = 0; i < c; i++)
                if (good[i] < e)
                  d1 += 1;
                else
                  d2 += 1;
              if (d1 > 1 || d2 > 1)
                {
#ifdef REPORT
                  printf("\nVariant %d %d\n",d1,d2);
                  j = 0;
                  for (i = 0; i < a; i++)
                    { if (j < c && advn[i] == good[j])
                        { if (good[j] < e)
                            printf("  D1: ");
                          else
                            printf("  D2: ");
                          j += 1;
                        }
                      else
                        printf("   *: ");
                      print_hap(finger[advn[i]],kmer,khalf);
                      printf(" %d\n",*COUNT_PTR(finger[advn[i]]));
                    }
#endif
                  if (d1 == 0)
                    { uhc2 += 1;        // unique het-context to d2
#ifdef REPORT
                      printf("    Het unique to D2\n");
#endif
                     }
                  else if (d2 == 0)
                    { uhc1 += 1;        // unique het-context to d2
#ifdef REPORT
                      printf("    Het unique to D1\n");
#endif
                     }
                  else if (d1 != d2)
                    { hpne += 1;        // het-pairs not equal
#ifdef REPORT
                      printf("    Hets not equal (a)\n");
#endif
                     }
                  else
                    { for (j = 0; j < d1; j++)
                        if (mycmp(finger[good[j]],finger[good[d1+j]],kbyte) != 0)
                          break;
                      if (j == d1)
                        { hpeq += 1;        // het-pairs equal
#ifdef REPORT
                          printf("    Hets are equal\n");
#endif
                        }
                      else
                        { hpne += 1;        // het-pairs not equal
#ifdef REPORT
                          printf("    Hets not equal (b)\n");
#endif
                        }
                    }
                }
            }

          for (i = 0; i < a; i++)
            finger[advn[i]] += tbyte;

#ifdef DEBUG_PARTITION
          for (i = 0; i < e; i++)
            printf(" %ld",(finger[i]-cache1)/tbyte);
          printf(" /");
          for (i = e; i < f; i++)
            printf(" %ld",(finger[i]-cache2)/tbyte);
          printf("\n");
#endif
        }
    }

  if (tid != 0)
    { Free_Kmer_Stream(T);
      Free_Kmer_Stream(U);
    }

  parm->uhc1 = uhc1;
  parm->uhc2 = uhc2;
  parm->hpne = hpne;
  parm->hpeq = hpeq;
  return (NULL);
}


/****************************************************************************************
 *
 *  Main
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ Kmer_Stream *T;
  Kmer_Stream *U;

  { int    i, j, k;
    int    flags[128];
    char  *eptr, *fptr;

    ARG_INIT("Hetrez");

    ETHRESH_1 = 1;
    ETHRESH_2 = 1;
    NTHREADS  = 4;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'e':
            ETHRESH_1 = strtol(argv[i]+2,&eptr,10);
            if (eptr > argv[i]+2)
              { if (ETHRESH_1 < 1 || ETHRESH_1 > 0x7fff)
                  { fprintf(stderr,"%s: Error count threshold %d is out of range\n",
                                   Prog_Name,ETHRESH_1);
                    exit (1);
                  }
                if (*eptr == ':')
                  { ETHRESH_2 = strtol(eptr+1,&fptr,10);
                    if (fptr > eptr+1 && *fptr == '\0')
                      { if (ETHRESH_2 < 1 || ETHRESH_2 > 0x7fff)
                          { fprintf(stderr,"%s: Error count threshold %d is out of range\n",
                                           Prog_Name,ETHRESH_2);
                            exit (1);
                          }
                        break;
                      }
                  }
              }
            fprintf(stderr,"%s: Syntax of -e option invalid -e<int>:<int>\n",
                           Prog_Name);
            exit (1);
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    if (argc != 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode\n");
        fprintf(stderr,"      -e: Accept only k-mers with counts >= threshold\n");
        fprintf(stderr,"      -T: number of threads to use\n");
        exit (1);
      }
  }

  T = Open_Kmer_Stream(argv[1]);
  U = Open_Kmer_Stream(argv[2]);

  if (T == NULL)
    { fprintf(stderr,"%s: Cannot open table %s\n",Prog_Name,argv[1]);
      exit (1);
    }
  if (U == NULL)
    { fprintf(stderr,"%s: Cannot open table %s\n",Prog_Name,argv[2]);
      exit (1);
    }
  if (T->kmer != U->kmer)
    { fprintf(stderr,"%s: Tables have different k-mer lengths, %d and %d, respectively\n",
                     Prog_Name,T->kmer,U->kmer);
      exit (1);
    }

  //  Analyze the pair of tables in thread blocks

  { Kmer_Stream  *T, *U;
    int64     range[NTHREADS+1][2];
#ifndef DEBUG_THREADS
    pthread_t threads[NTHREADS];
#endif
    TP        parm[NTHREADS];
    char     *seq;
    uint8    *ent;
    int       t, i;
    int64     p;
  
    T = Open_Kmer_Stream(argv[1]);
    U = Open_Kmer_Stream(argv[2]);

    range[0][0] = 0;
    range[NTHREADS][0] = T->nels;
    range[0][1] = 0;
    range[NTHREADS][1] = U->nels;

    seq = Current_Kmer(T,NULL);
    ent = Current_Entry(T,NULL);
    for (t = 1; t < NTHREADS; t++)
      { p = (T->nels*t)/NTHREADS; 
        GoTo_Kmer_Index(T,p);
#ifdef DEBUG
        printf("\n%d: %0*x\n",t,2*T->ibyte,T->cpre);
        printf(" %lld: %s\n",p,Current_Kmer(T,seq));
#endif
        ent = Current_Entry(T,ent);                //  Break at prefix boundaries
        for (i = T->ibyte; i < T->kbyte; i++)
          ent[i] = 0;
        GoTo_Kmer_Entry(T,ent);
#ifdef DEBUG
        printf(" %lld: %s\n",T->cidx,Current_Kmer(T,seq));
#endif
        range[t][0] = T->cidx;
        GoTo_Kmer_Entry(U,ent);
#ifdef DEBUG
        printf(" %lld: %s\n",U->cidx,Current_Kmer(U,seq));
#endif
        range[t][1] = U->cidx;
      }
    free(seq);

    for (t = 0; t < NTHREADS; t++)
      { parm[t].tid   = t;
        parm[t].T     = T;
        parm[t].U     = U;
        parm[t].begs  = range[t];
        parm[t].ends  = range[t+1];
      }

#ifdef DEBUG_THREADS
    for (t = 0; t < NTHREADS; t++)
      analysis_thread(parm+t);
#else
    for (t = 1; t < NTHREADS; t++)
      pthread_create(threads+t,NULL,analysis_thread,parm+t);
    analysis_thread(parm);
    for (t = 1; t < NTHREADS; t++)
      pthread_join(threads[t],NULL);
#endif

    for (t = 1; t < NTHREADS; t++)
      { parm[0].uhc1 += parm[t].uhc1;
        parm[0].uhc2 += parm[t].uhc2;
        parm[0].hpne += parm[t].hpne;
        parm[0].hpeq += parm[t].hpeq;
      }

    printf("\nHet-contexts unique to source1: ");
    Print_Number(parm[0].uhc1,10,stdout);
    printf("\nHet-contexts unique to source2: ");
    Print_Number(parm[0].uhc2,10,stdout);
    printf("\nDistinct Het-pairs: ");
    Print_Number(parm[0].hpne,10,stdout);
    printf("\nShared Het-pairs:   ");
    Print_Number(parm[0].hpeq,10,stdout);
    printf("\n");

    Free_Kmer_Stream(T);
    Free_Kmer_Stream(U);
  }

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
