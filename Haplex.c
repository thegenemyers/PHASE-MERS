/*********************************************************************************************\
 *
 *  Code for produce all potential haplotype k-mers with a single SNP in the center of
 *    the k-mers
 *
 *  Author:  Gene Myers
 *  Date  :  October, 2020
 *
 *********************************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>

#undef  DEBUG_PARTITION
#define DEBUG_TABLE

#include "libfastk.h"

static char *Usage = " [-h<int>:<int>] <source1>[.ktab] <source2>[.ktab]";

static int   HAPLO_LOW, HAPLO_HGH;   //  Count range for acceptable het k-mers
static int   VERBOSE;

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
        printf("*");
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


/****************************************************************************************
 *
 *  Find Haplotype Pairs
 *
 *****************************************************************************************/

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

static void Analyze_Het_Mers(Kmer_Stream *T, Kmer_Stream *U)
{ uint8  prefs[] = { 0x3f, 0x0f, 0x03, 0x00 };

  uint8 *cache1, *cptr1, *ctop1, *nptr1;
  uint8 *cache2, *cptr2, *ctop2, *nptr2;
  uint8 *start;

  int    e, f, i;
  int    index[8];
  uint8 *finger[8];
  uint8 *flimit[8];

  int    a, advn[4];
  int    c, good[4];
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

  First_Kmer_Entry(T);
  First_Kmer_Entry(U);
  while (T->csuf != NULL && U->csuf != NULL)
    { f = 0;
      index[f++] = 0;

      cptr1 = nptr1 = cache1;
      Current_Entry(T,cptr1);

      cptr2 = nptr2 = cache2;
      Current_Entry(U,cptr2);

      if (mycmp(cptr1,cptr2,kbyte) < 0)
        start = cptr1;
      else
        start = cptr2;

      for ( ; T->csuf != NULL; Next_Kmer_Entry(T))
        { int x = mypref(start,Current_Entry(T,nptr1),khalf); 
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

      for ( ; U->csuf != NULL; Next_Kmer_Entry(U))
        { int x = mypref(start,Current_Entry(T,nptr2),khalf); 
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
      printf("part %d"",f);
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

#define ADD(i)					\
{ int cn = *COUNT_PTR(finger[i]);		\
  advn[a++] = i;				\
  if (HAPLO_LOW <= cn && cn <= HAPLO_HGH)	\
    good[c++] = i;				\
}

#define SET(i)	\
{ mc = hc;	\
  mr = hr;	\
  a = c = 0;	\
  ADD(i);	\
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

              d1 = d2 = 0;
              for (i = 0; i < c; i++)
                if (good[i] < e)
                  d1 += 1;
                else
                  d2 += 1;
              if (d1 > 1 || d2 > 1)
                { if (d1 == 0)
                    uhc2 += 1;        // unique het-context to d2
                  else if (d2 == 0)
                    uhc1 += 1;        // unique het-context to d2
                  else if (d1 != d2)
                    hpne += 1;        // het-pairs not equal
                  else
                    { for (j = 0; j < d1; j++)
                        if (mycmp(finger[good[j]],finger[good[d1+j]],kbyte) != 0)
                          break;
                      if (j == d1)
                        hpeq += 1;        // het-pairs equal
                      else
                        hpne += 1;        // het-pairs not equal
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

  if (VERBOSE)
    { printf("\nHet-contexts unique to source1: ");
      Print_Number(uhc1,10,stdout);
      printf("\nHet-contexts unique to source2: ");
      Print_Number(uhc2,10,stdout);
      printf("\nDistinct Het-pairs: ");
      Print_Number(hpne,10,stdout);
      printf("\nShared Het-pairs:   ");
      Print_Number(hpeq,10,stdout);
      printf("\n");
    }
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

    ARG_INIT("Haplex");

    HAPLO_LOW = 1;
    HAPLO_HGH = 0x7fff;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'h':
            HAPLO_LOW = strtol(argv[i]+2,&eptr,10);
            if (eptr > argv[i]+2)
              { if (HAPLO_LOW < 1 || HAPLO_LOW > 0x7fff)
                  { fprintf(stderr,"%s: Haplo minimum count %d is out of range\n",
                                   Prog_Name,HAPLO_LOW);
                    exit (1);
                  }
                if (*eptr == ':')
                  { HAPLO_HGH = strtol(eptr+1,&fptr,10);
                    if (fptr > eptr+1 && *fptr == '\0')
                      { if (HAPLO_HGH < 1 || HAPLO_HGH > 0x7fff)
                          { fprintf(stderr,"%s: Haplo maximum count %d is out of range\n",
                                           Prog_Name,HAPLO_HGH);
                            exit (1);
                          }
                        if (HAPLO_LOW > HAPLO_HGH)
                          { fprintf(stderr,"%s: Haplo count range is invalid\n",Prog_Name);
                            exit (1);
                          }
                        break;
                      }
                  }
              }
            fprintf(stderr,"%s: Syntax of -h option invalid -h<int>:<int>\n",Prog_Name);
            exit (1);
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode\n");
        fprintf(stderr,"      -h: Accept only het-mers with count in given range (inclusive).\n");
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

  Analyze_Het_Mers(T,U);

  Free_Kmer_Stream(U);
  Free_Kmer_Stream(T);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
