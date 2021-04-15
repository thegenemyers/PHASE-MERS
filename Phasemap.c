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
#include "DB.h"

#undef DEBUG

static char *Usage = "[-D<dazz_db>] <lower>:[.prof] <upper>[.prof]";

static char DNA[4] = { 'a', 'c', 'g', 't' };

static void Scan_All(Profile_Index *LOW, Profile_Index *HGH, DAZZ_DB *DBP, Profile_Index *DPR)
{ DAZZ_READ *read;
  uint8     *rseq, *radj;
  uint16    *lowprof, *hghprof, *dbsprof;
  int        pmax, plen;
  int        nreads, kmer, km1;
  int        i, id;

  pmax    = 20000;
  lowprof = Malloc(3*pmax*sizeof(uint16),"Profile array");
  hghprof = lowprof + pmax;
  dbsprof = hghprof + pmax;

  kmer    = LOW->kmer;
  km1     = kmer - 1;
  nreads  = LOW->nbase[LOW->nparts-1];
  if (DBP != NULL)
    { rseq = (uint8 *) New_Read_Buffer(DBP);
      radj = rseq+km1;
      read = DBP->reads;
    }

  for (id = 0; id < nreads; id++) 
    { plen = Fetch_Profile(LOW,(int64) id,pmax,lowprof);
      if (plen > pmax)
        { pmax    = 1.2*plen + 1000;
          lowprof = Realloc(lowprof,3*pmax*sizeof(uint16),"Profile array");
          hghprof = lowprof + pmax;
          dbsprof = hghprof + pmax;
          Fetch_Profile(LOW,(int64) id,pmax,lowprof);
        }
      if (Fetch_Profile(HGH,(int64) id,pmax,hghprof) != plen)
        { fprintf(stderr,"%s: Low & Hgh profiles for id %d are not the same length !?\n",
                         Prog_Name,id+1);
          exit (1);
        }
      if (DBP != NULL)
        { Load_Read(DBP,id,(char *) rseq,0);
          if (read[id].rlen != plen + km1)
            { fprintf(stderr,"%s: Profiles and DB sequence lengths do not correspond\n",
                             Prog_Name);
              exit (1);
            }
          if (DPR != NULL)
            { if (Fetch_Profile(DPR,(int64) id,pmax,dbsprof) != plen)
                { fprintf(stderr,"%s: Het-mer and DB profiles for id %d are not the same length?\n",
                                 Prog_Name,id+1);
                  exit (1);
                }
            }
        }

#ifdef DEBUG
      printf("\nRead %d:\n",id+1);
      if (DBP != NULL)
        { for (i = 0; i < km1; i++)
            printf(" %5d: %c\n",i-km1,DNA[rseq[i]]);
          for (i = 0; i < plen; i++)
            { uint32 hgh = hghprof[i];
              uint32 low = lowprof[i];
              uint32 code = (hgh << 16 | low);

              if (code > 0)
                { if ((code & 0x7) == 0)
                    printf(" %5d: %c %10d -",i,DNA[radj[i]],code>>3);
                  else
                    printf(" %5d: %c %10d %c",i,DNA[radj[i]],code>>3,DNA[(code&0x7)-1]);
                }
              else
                printf(" %5d: %c 0",i,DNA[radj[i]]);
              if (DPR != NULL)
                printf(" [%5d]",dbsprof[i]);
              printf("\n");
            }
        }
      else
        for (i = 0; i < plen; i++)
          { uint32 hgh = hghprof[i];
            uint32 low = lowprof[i];
            uint32 code = (hgh << 16 | low);

            if (code > 0)
              { if ((code & 0x7) == 0)
                  printf(" %5d: %10d -\n",i,code>>3);
                else
                  printf(" %5d: %10d %c\n",i,code>>3,DNA[(code&0x7)-1]);
              }
            else
              printf(" %5d: 0\n",i);
          }
#else
      printf("\nRead\t%d\n",id+1);
      for (i = 0; i < plen; i++)
        { uint32 hgh = hghprof[i];
          uint32 low = lowprof[i];
          uint32 code = (hgh << 16 | low);

          if (code > 0)
            { if ((code & 0x7) == 0)
                printf("\t%d\t%d\t-",i+kmer,code>>3);
              else
                printf("\t%d\t%d\t%c",i+kmer,code>>3,DNA[(code&0x7)-1]);
              if (DPR != NULL)
                printf("\t%5d",dbsprof[i]);
              printf("\n");
            }
        }
#endif
    }

  free(lowprof);
}


/****************************************************************************************
 *
 *  Main
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ Profile_Index *LOW;
  Profile_Index *HGH;
  Profile_Index *DPR;
  DAZZ_DB        DB, *DBP;

  { int    i, j, k;
    int    flags[128];

    (void) flags;

    ARG_INIT("Phaselink");

    DBP = NULL;
    DPR = NULL;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
          case 'D':
            if (Open_DB(argv[i]+2,&DB) < 0)
              exit (1);
            DBP = &DB;
            DPR = Open_Profiles(argv[i]+2);
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc != 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
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

  Scan_All(LOW,HGH,DBP,DPR);

  Free_Profiles(HGH);
  Free_Profiles(LOW);
  if (DBP != NULL)
    Close_DB(DBP);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
