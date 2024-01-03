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

static char *Usage = " [-T<int(4)>] [-h<int>:<int>::<int>:<int>] <source1>[.ktab] <source2>[.ktab]";

static int   HAPLO_LOW1, HAPLO_HGH1;   //  Count ranges for acceptable het k-mers
static int   HAPLO_LOW2, HAPLO_HGH2;
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
    int64        sx, sy, sx2, sy2, sxy, sn;
    int64        hx, hy, hx2, hy2, hxy, hn;
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

  int    e, f;
  int    index[8];
  uint8 *finger[8];
  uint8 *flimit[8];

  int64  sx, sy, sx2, sy2, sxy, sn;
  int64  hx, hy, hx2, hy2, hxy, hn;
#ifdef REPORT
  double scor, hcor;
#endif

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

  sx2 = sx = sy2 = sy = sn = sxy = 0;
  hx2 = hx = hy2 = hy = hn = hxy = 0;

#ifdef DEBUG_THREADS
  printf("Doing %d:",tid);
  printf(" [%lld-%lld]",begs[0],parm->ends[0]);
  printf(" [%lld-%lld]",begs[1],parm->ends[1]);
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
      { int i;

        printf("part %d:",f);
        for (i = 0; i < e; i++)
          printf(" %d",index[i]/tbyte);
        printf(" /");
        for (i = e; i < f; i++)
          printf(" %d",index[i]/tbyte);
        printf(" == %ld %ld\n",(nptr1-cache1)/tbyte,(nptr2-cache2)/tbyte);
      }
#endif

      if (e <= 1 && f-e <= 1)
        continue;

      { int i;

        for (i = 0; i < f; i++)
          if (i < e)
            finger[i] = cache1 + index[i];
          else
            finger[i] = cache2 + index[i];
        for (i = 1; i < f; i++)
          flimit[i-1] = finger[i];
        flimit[e-1] = nptr1;
        flimit[f-1] = nptr2;
      }

      //  Accumulate count correlation of all k-mers good in either table

      { uint8 *t, *et;
        uint8 *u, *eu;
        int    x, c, d;

        et = flimit[e-1];
        eu = flimit[f-1];
        t  = finger[0];
        u  = finger[e];
        while (1)
          { if (t >= et)
              { if (u >= eu)
                  break;
                x = 1;
              }
            else if (u >= eu)
              x = -1;
            else
              x = mycmp(t,u,kbyte);
            if (x < 0)
              { c = *COUNT_PTR(t);
                if (c >= HAPLO_LOW1 && c <= HAPLO_HGH1)
                  { sx2 += c*c;
                    sx  += c;
                    sn  += 1;
                  }
                t += tbyte;
              }
            else if (x > 0)
              { d = *COUNT_PTR(u);
                if (d >= HAPLO_LOW2 && d <= HAPLO_HGH2)
                  { sy2 += d*d;
                    sy  += d;
                    sn  += 1;
                  }
                u += tbyte;
              }
            else
	      { c = *COUNT_PTR(t);
                d = *COUNT_PTR(u);
                if ((c >= HAPLO_LOW1 && c <= HAPLO_HGH1) || (d >= HAPLO_LOW2 && d <= HAPLO_HGH2))
                  { sx2 += c*c;
                    sx  += c;
                    sy2 += d*d;
                    sy  += d;
                    sxy += c*d;
                    sn  += 1;
                  }
                t += tbyte;
                u += tbyte;
              }
          }
#ifdef REPORT
        scor = (sn*sxy - sx*sy) / (sqrt(sn*sx2 - sx*sx) * sqrt(sn*sy2 - sy*sy));
        printf(" %lld %lld %lld %lld %lld %lld = %.4f\n",sn,sx,sy,sx2,sy2,sxy,scor);
#endif
      }

      //  Accumulate count correlation of all het k-mers good in either table

#define ADD(i)                                  	\
{ int cn = *COUNT_PTR(finger[i]);               	\
  advn[a] = i; 		                        	\
  if (i < e)                                    	\
    good[a++] = (cn >= HAPLO_LOW1 && cn <= HAPLO_HGH1);	\
  else                                          	\
    good[a++] = (cn >= HAPLO_LOW2 && cn <= HAPLO_HGH2);	\
}

#define SET(i)  \
{ mc = hc;      \
  mr = hr;      \
  a = 0;        \
  ADD(i);       \
}

      { int    a, b, advn[8], good[8];
        int    mc, hc;
        uint8 *mr, *hr;
        int    d1, d2;
        int    i, j, x;
        int    c, d;

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
            for (i = 0; i < a; i++)
              { printf(" %d",advn[i]);
                if (good[i])
                  printf("+");
              }
            printf("\n");
#endif
  
            d1 = d2 = 0;
            for (i = 0; i < a; i++)
              { if (advn[i] >= e)
                  break;
                d1 += good[i];
              }
            for (b = i; i < a; i++)
              d2 += good[i];
 
            if (d1 > 1 || d2 > 1)
              {
#ifdef REPORT
                printf("\nVariant %d %d\n",d1,d2);
                for (i = 0; i < a; i++)
                  { if (good[i])
                      { if (i < b)
                          printf("  D1: ");
                        else
                          printf("  D2: ");
                      }
                    else
                      printf("   *: ");
                    print_hap(finger[advn[i]],kmer,khalf);
                    printf(" %d\n",*COUNT_PTR(finger[advn[i]]));
                  }
#endif

                i = 0;
                j = b;
                while (1)
                  { if (i >= b)
                      { if (j >= a)
                          break;
                        x = 1;
                      }
                    else if (j >= a)
                      x = -1;
                    else
                      x = mycmp(finger[advn[i]],finger[advn[j]],kbyte);
                    if (x < 0)
                      { if (good[i])
                          { c = *COUNT_PTR(finger[advn[i]]);
                            hx2 += c*c;
                            hx  += c;
                            hn  += 1;
                          }
printf(" <%d",i);
                        i += 1;
                      }
                    else if (x > 0)
                      { if (good[j])
                          { d = *COUNT_PTR(finger[advn[j]]);
                            hy2 += d*d;
                            hy  += d;
                            hn  += 1;
                          }
printf(" >%d",j);
                        j += 1;
                      }
                    else
                      { if (good[i] || good[j])
                          { c = *COUNT_PTR(finger[advn[i]]);
                            d = *COUNT_PTR(finger[advn[j]]);
                            hx2 += c*c;
                            hx  += c;
                            hy2 += d*d;
                            hy  += d;
                            hxy += c*d;
                            hn  += 1;
                          }
printf(" %d+%d",i,j);
                        i += 1;
                        j += 1;
                      }
                  }
printf("\n");
#ifdef REPORT
                hcor = (hn*hxy - hx*hy) / (sqrt(hn*hx2 - hx*hx) * sqrt(hn*hy2 - hy*hy));
                printf(" %lld %lld %lld %lld %lld %lld = %.4f\n",hn,hx,hy,hx2,hy2,hxy,hcor);
#endif
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
    }

  if (tid != 0)
    { Free_Kmer_Stream(T);
      Free_Kmer_Stream(U);
    }

  parm->sn  = sn;
  parm->sx  = sx;
  parm->sy  = sy;
  parm->sx2 = sx2;
  parm->sy2 = sy2;
  parm->sxy = sxy;

  parm->hn  = hn;
  parm->hx  = hx;
  parm->hy  = hy;
  parm->hx2 = hx2;
  parm->hy2 = hy2;
  parm->hxy = hxy;

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

    ARG_INIT("Hetcor");

    HAPLO_LOW1 = 1;
    HAPLO_HGH1 = 0x7fff;
    HAPLO_LOW2 = 1;
    HAPLO_HGH2 = 0x7fff;
    NTHREADS  = 4;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'h':
            HAPLO_LOW1 = strtol(argv[i]+2,&eptr,10);
            if (eptr > argv[i]+2)
             { if (HAPLO_LOW1 < 1 || HAPLO_LOW1 > 0x7fff)
                { fprintf(stderr,"%s: First haplo minimum count %d is out of range\n",
                                 Prog_Name,HAPLO_LOW1);
                  exit (1);
                }
               if (*eptr == ':')
                { HAPLO_HGH1 = strtol(eptr+1,&fptr,10);
                  if (fptr > eptr+1 && *fptr == ':')
                   { if (HAPLO_HGH1 < 1 || HAPLO_HGH1 > 0x7fff)
                      { fprintf(stderr,"%s: First haplo maximum count %d is out of range\n",
                                       Prog_Name,HAPLO_HGH1);
                        exit (1);
                      }
                     if (HAPLO_LOW1 > HAPLO_HGH1)
                      { fprintf(stderr,"%s: First haplo count %d:%d range is invalid\n",
                                       Prog_Name,HAPLO_LOW1,HAPLO_HGH1);
                        exit (1);
                      }
                     if (*fptr != ':' || fptr[1] != ':')
                      { fprintf(stderr,"%s: Expecting :: between ranges, not '%c%c'\n",
                                       Prog_Name,fptr[0],fptr[1]);
                        exit (1);
                      }
                     HAPLO_LOW2 = strtol(fptr+2,&eptr,10);
                     if (eptr > fptr+2)
                      { if (HAPLO_LOW2 < 1 || HAPLO_LOW2 > 0x7fff)
                         { fprintf(stderr,"%s: Second haplo minimum count %d is out of range\n",
                                          Prog_Name,HAPLO_LOW2);
                           exit (1);
                         }
                        if (*eptr == ':')
                         { HAPLO_HGH2 = strtol(eptr+1,&fptr,10);
                           if (fptr > eptr+1 && *fptr == '\0')
                            { if (HAPLO_HGH2 < 1 || HAPLO_HGH2 > 0x7fff)
                               { fprintf(stderr,"%s: Second haplo maximum count",Prog_Name);
                                 fprintf(stderr," %d is out of range\n",HAPLO_HGH2);
                                 exit (1);
                               }
                              if (HAPLO_LOW2 > HAPLO_HGH2)
                               { fprintf(stderr,"%s: Second haplo count %d:%d range is invalid\n",
                                                Prog_Name,HAPLO_LOW2,HAPLO_HGH2);
                                 exit (1);
                               }
                              break;
                            }
                         }
                      }
                   }
                }
              }
            fprintf(stderr,"%s: Syntax of -h option invalid -h<int>:<int>::<int>:<int>\n",
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
        fprintf(stderr,"      -h: Accept only k-mers with count in given range (inclusive).\n");
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

    { double  sx, sy, sx2, sy2, sxy, sn;
      double  hx, hy, hx2, hy2, hxy, hn;
      double  scor, hcor;

      sx = sy = sx2 = sy2 = sxy = sn = 0.;
      hx = hy = hx2 = hy2 = hxy = hn = 0.;
      for (t = 0; t < NTHREADS; t++)
        { sn  += parm[t].sn;
          sx  += parm[t].sx;
          sy  += parm[t].sy;
          sx2 += parm[t].sx2;
          sy2 += parm[t].sy2;
          sxy += parm[t].sxy;

          hn  += parm[t].hn;
          hx  += parm[t].hx;
          hy  += parm[t].hy;
          hx2 += parm[t].hx2;
          hy2 += parm[t].hy2;
          hxy += parm[t].hxy;
        }

      scor = (sn*sxy - sx*sy) / (sqrt(sn*sx2 - sx*sx) * sqrt(sn*sy2 - sy*sy));
      hcor = (hn*hxy - hx*hy) / (sqrt(hn*hx2 - hx*hx) * sqrt(hn*hy2 - hy*hy));

      printf("\nCorrelations:\n");
      printf("  Entire genome %.4f\n",scor);
      printf("  Good het-mers %.4f\n",hcor);
    } 

    Free_Kmer_Stream(T);
    Free_Kmer_Stream(U);
  }

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
