/*********************************************************************************************\
 *
 *  Example code for opening and fetching compressed profiles produced by FastK
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

#include "libfastk.h"
#include "reader.h"

#undef DEBUG

static int BC_PREFIX;

static char *Usage[] = { " [-10X] <lower>:[.prof] <upper>[.prof]",
                         "    <source>[.cram|.[bs]am|.db|.dam|.f[ast][aq][.gz] ..."
                       };

static char DNA[4] = { 'A', 'C', 'G', 'T' };

static uint8 code[128] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

static uint8 comp[128] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 3, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 3, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

static int is_minimal(char *seq, int len)
{ int j, k;
  int x, y;

  for (k = 0, j = len-1; k < j; k++, j--)
    { x = code[(int) seq[k]];
      y = comp[(int) seq[j]];
      if (x < y)
        return (1);
      if (x > y)
        return (0);
    }
  if (k <= j)
    { x = code[(int) seq[k]];
      if (x < 2)
        return (1);
      else
        return (0);
    }
  else
    return (1);
}

static int64 Scan_File(Profile_Index *LOW, Profile_Index *HGH, char *input, int64 beg) 
{ char      *rseq;
  uint16    *lowprof, *hghprof;
  int        pmax, plen;
  int        nreads, kmer, km1;
  int        i, id;

  Fetcher    fetch;
  Entry     *read;

  pmax    = 20000;
  lowprof = Malloc(2*pmax*sizeof(uint16),"Profile array");
  hghprof = lowprof + pmax;

  kmer    = LOW->kmer;
  km1     = kmer - 1;
  nreads  = LOW->nbase[LOW->nparts-1];

  fetch = Open_File(input);
  if (fetch == NULL)
    { fprintf(stderr,"%s: Cannot recognize input file %s\n",Prog_Name,input);
      exit (1);
    }

  for (id = beg; id < nreads; id++) 
    { read = fetch();
      if (read == NULL)
        return (id);
      rseq = read->seq;

      plen = Fetch_Profile(LOW,(int64) id,pmax,lowprof);
      if (plen > pmax)
        { pmax    = 1.2*plen + 1000;
          lowprof = Realloc(lowprof,2*pmax*sizeof(uint16),"Profile array");
          hghprof = lowprof + pmax;
          Fetch_Profile(LOW,(int64) id,pmax,lowprof);
        }
      if (Fetch_Profile(HGH,(int64) id,pmax,hghprof) != plen)
        { fprintf(stderr,"%s: Low & Hgh profiles for id %d are not the same length !?\n",
                         Prog_Name,id+1);
          exit (1);
        }
      if ((int) read->len != plen+km1)
        { fprintf(stderr,"%s: Sequence and profile not the same length?\n",
                         Prog_Name);
          exit (1);
        }

#ifdef DEBUG
      printf("\nRead %d:\n",id+1);
      for (i = 0; i < BC_PREFIX+km1; i++)
        { if (i == 16)
            printf(" ---------\n");
          printf(" %5d: %c\n",i-(BC_PREFIX+km1),rseq[i]);
        }
      for (i = BC_PREFIX; i < plen; i++)
        { uint32 hgh = hghprof[i];
          uint32 low = lowprof[i];
          uint32 code = (hgh << 16 | low);

          if (code > 0)
            { if ((code & 0x7) == 0)
                printf(" %5d: %c %10d -",i,rseq[i+km1],code>>3);
              else
                printf(" %5d: %c %10d %c",i,rseq[i+km1],code>>3,DNA[(code&0x7)-1]);
              if (is_minimal(rseq+i,kmer))
                printf("\t+\n");
              else
                printf("\t-\n");
            }
          else
            printf(" %5d: %c 0",i,rseq[i+km1]);
          printf("\n");
        }
#else

      printf("R\t%d\t",id+1);
      if (BC_PREFIX == 0)
        printf("-\n");
      else
        { for (i = 0; i < 16; i++)
            putchar(rseq[i]);
          putchar('\n');
        }

      for (i = BC_PREFIX; i < plen; i++)
        { uint32 hgh = hghprof[i];
          uint32 low = lowprof[i];
          uint32 code = (hgh << 16 | low);

          if (code > 0)
            { if ((code & 0x7) == 0)
                printf("H\t%d\t%10d",i+kmer,code>>3);
              else
                printf("%c\t%d\t%10d",DNA[(code&0x7)-1],i+kmer,code>>3);
              if (is_minimal(rseq+i,kmer))
                printf("\t+\n");
              else
                printf("\t-\n");
            }
        }
#endif
    }

  if (fetch() != NULL)
    { fprintf(stderr,"%s: More reads in input then in profiles??\n",Prog_Name);
      exit (1);
    }

  free(lowprof);

  return (-1);
}


/****************************************************************************************
 *
 *  Main
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ Profile_Index *LOW;
  Profile_Index *HGH;

  { int    i, j, k;
    int    flags[128];

    (void) flags;

    ARG_INIT("Phasemap");

    BC_PREFIX = 0;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
          case '1':
            if (argv[i][2] != '0' || argv[i][3] != 'X')
              { fprintf(stderr,"\n%s: -%s is not a legal optional argument\n",Prog_Name,argv[i]);
                exit (1);
              }
            BC_PREFIX = 23;
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc < 4)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"\n");
        fprintf(stderr,"    -10X: Each read begins with a 10X bar code and linker\n");
        exit (1);
      }
  }

  LOW = Open_Profiles(argv[1]);
  if (LOW == NULL)
    { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,argv[1]);
      exit (1);
    }
  HGH = Open_Profiles(argv[2]);
  if (HGH == NULL)
    { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,argv[2]);
      exit (1);
    }

  { int   c;
    int64 beg;

    beg = 0;
    for (c = 3; c < argc; c++)
      { if (beg < 0)
          { fprintf(stderr,"%s: More reads in input then in profiles??\n",Prog_Name);
            exit (1);
          }
        beg = Scan_File(LOW,HGH,argv[c],beg);
      }
    if (beg >= 0)
      { fprintf(stderr,"%s: More reads in profiles then in input??\n",Prog_Name);
        exit (1);
      }
  }

  Free_Profiles(HGH);
  Free_Profiles(LOW);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
