/********************************************************************************************
 *
 *  General sequence file reader
 *
 *  Author:  Gene Myers
 *  Date  :  April, 2021
 *
 *********************************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <inttypes.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>
#include <pthread.h>

#include <zlib.h>

#include "libfastk.h"
#include "sam.h"
#include "reader.h"


/*******************************************************************************************
 *
 *  Read fasta or fastq sequences
 *
 ********************************************************************************************/

static FILE   *fid;
static gzFile  gid;
static char *(*read_line)(int64 *);

//  Read next line into a buffer and return a pointer to the buffer and set *plen
//    the length of the line.  NB: replaces '\n' with '\0'.

static char *Read_Line(int64 *plen)
{ static char *buffer;
  static int   bmax = 0;
  int len;
  
  if (bmax == 0)
    { bmax = 500;
      buffer = (char *) Malloc(bmax,"Allocating line buffer");
      if (buffer == NULL)
        exit (1);
    }
  
  if (fgets(buffer,bmax,fid) == NULL)
    return (NULL);
  
  len = strlen(buffer);
  while (buffer[len-1] != '\n') 
    { bmax = ((int) (1.4*bmax)) + 100;
      buffer = (char *) Realloc(buffer,bmax,"Reallocating line buffer");
      if (buffer == NULL)
        exit (1);
      if (fgets(buffer+len,bmax-len,fid) == NULL)
        { fprintf(stderr,"%s: Last line of file does not end with new-line\n",Prog_Name);
          exit (1);
        } 
      len += strlen(buffer+len);
    }
  buffer[--len] = '\0';
  
  *plen = len;
  return (buffer);
}

//  Read next line into a buffer and return a pointer to the buffer and set *plen
//    the length of the line.  NB: replaces '\n' with '\0'.

static char *Read_Zip_Line(int64 *plen)
{ static char *buffer;
  static int   bmax = 0;
  int len;
  
  if (bmax == 0)
    { bmax = 500;
      buffer = (char *) Malloc(bmax,"Allocating line buffer");
      if (buffer == NULL)
        exit (1);
    }
  
  if (gzgets(gid,buffer,bmax) == NULL)
    return (NULL);
  
  len = strlen(buffer);
  while (buffer[len-1] != '\n') 
    { bmax = ((int) (1.4*bmax)) + 100;
      buffer = (char *) Realloc(buffer,bmax,"Reallocating line buffer");
      if (buffer == NULL)
        exit (1);
      if (gzgets(gid,buffer+len,bmax-len) == NULL)
        { fprintf(stderr,"%s: Last line of file does not end with new-line\n",Prog_Name);
          exit (1);
        } 
      len += strlen(buffer+len);
    }
  buffer[--len] = '\0';
  
  *plen = len;
  return (buffer);
}

//  Read next fastq entry returning it in Entry record

static Entry *Get_Fastq_Entry()
{ static Entry   entry;
  static int64   smax = 0;
  static char   *seq = NULL;

  int64   slen, qlen;
  char   *line;

  line = read_line(&qlen);
  if (line == NULL)
     return (NULL);

  if (line[0] != '@')
    { fprintf(stderr,"%s: Entry header does not start with an @-sign\n",Prog_Name);
      exit (1);
    }

  line = read_line(&slen);
  if (line == NULL)
    { fprintf(stderr,"%s: Incomplete fastq entry\n",Prog_Name);
      exit (1);
    }
  if (slen > smax)
    { smax = 1.2*slen + 10000; 
      seq = Realloc(seq,smax+1,"Allocating sequence buffer");
      if (seq == NULL)
        exit (1);
    }
  memcpy(seq,line,slen+1);

  line = read_line(&qlen);
  if (line == NULL)
    { fprintf(stderr,"scaff2contig: Incomplete fastq entry\n");
      exit (1);
    }
  if (line[0] != '+')
    { fprintf(stderr,"scaff2contig: Divider line does not start with a +-sign\n");
      exit (1);
    }

  line = read_line(&qlen);
  if (line == NULL)
    { fprintf(stderr,"scaff2contig: Incomplete fastq entry\n");
      exit (1);
    }
  if (slen != qlen)
    { fprintf(stderr,"scaff2contig: QV line does not have the same length as sequence line\n");
      exit (1);
    }

  entry.seq    = seq;
  entry.len    = slen;
  return (&entry);
}

//  Read next fasta entry returning it in Entry record

static Entry *Get_Fasta_Entry()
{ static Entry  entry;
  static int64  smax = 0;
  static char  *seq = NULL;
  static char  *next = NULL;
  static int    empty = 1;

  int64 slen, m;
  char *line;

  if (empty)
    { empty = 0;
      next = read_line(&slen);
      if (next == NULL)
        { fprintf(stderr,"%s: fasta is empty!\n",Prog_Name);
          exit (1);
        }
      if (next[0] != '>')
        { fprintf(stderr,"%s: Entry header does not start with a >-sign\n",Prog_Name);
          exit (1);
        }
    }
  else
    { if (next == NULL)
        return (NULL);
    }

  slen = 0;
  while (1)
    { line = read_line(&m);
      if (line == NULL || line[0] == '>')
        break;

      if (slen+m > smax)
        { smax = 1.2*(slen+m) + 1000000; 
          seq = Realloc(seq,smax+1,"Allocating sequence buffer");
          if (seq == NULL)
            exit (1);
        }

      memcpy(seq+slen,line,m);
      slen += m;
    }
  if (slen == 0)
    { if (line == NULL)
        fprintf(stderr,"%s: Sequence missing after last header\n",Prog_Name);
      else
        fprintf(stderr,"%s: Sequence missing between two headers\n",Prog_Name);
      exit (1);
    }

  next = line;
  entry.seq    = seq;
  entry.len    = slen;
  return (&entry);
}


/*******************************************************************************************
 *
 *  Read BAM/SAM file
 *
 ********************************************************************************************/

static samFile *sid;

static Entry *Get_Sam_Entry()
{ static Entry entry;

  samRecord *same;

  if (sam_eof(sid))
    { sam_close(sid);
      return (NULL);
    }

  same = sam_record_extract(sid,0);

  entry.seq = same->seq;
  entry.len = same->len;
  return (&entry);
}



/*******************************************************************************************
 *
 *  Read Cram file
 *
 ********************************************************************************************/

/*
static cram_fd *cid;

static Entry *Get_Cram_Entry()
{ static Entry entry;

  cram_record *rec;

  rec = cram_get_seq(cid);
  if (rec == NULL)
    { cram_close(cid);
      return (NULL);
    }

  entry.seq = (char *) (rec->s->seqs_blk->data+rec->seq);
  entry.len = rec->len;
  return (&entry);
}
*/


/*******************************************************************************************
 *
 *  Read Dazzler Database
 *
 ********************************************************************************************/

static DAZZ_DB db;

static Entry *Get_DB_Entry()
{ static Entry   entry;
  static int     rnum = 0;

  if (rnum == 0)
    entry.seq = New_Read_Buffer(&db);

  if (rnum >= db.nreads)
    { free(entry.seq);
      rnum = 0;
      Close_DB(&db);
      return (NULL);
    }

  Load_Read(&db,rnum,entry.seq,1);
  entry.len = db.reads[rnum].rlen;

  rnum += 1;
  return (&entry);
}


/*******************************************************************************************
 *
 *  Open a file and return it fid and type
 *
 ********************************************************************************************/

  //  Open and get info about each input file

Fetcher Open_File(char *arg)

{ static char *suffix[] = { ".cram", ".bam", ".sam", ".db", ".dam",
                            ".fastq", ".fasta", ".fq", ".fa",
                            ".fastq.gz", ".fasta.gz", ".fastq", ".fasta",
                            ".fq.gz",  ".fa.gz", ".fq", ".fa" };
  static char *extend[] = { ".cram", ".bam", ".sam", ".db", ".dam",
                            ".fastq", ".fasta", ".fq", ".fa",
                            ".fastq.gz", ".fasta.gz", ".fastq.gz", ".fasta.gz",
                            ".fq.gz",  ".fa.gz", ".fq.gz", ".fa.gz" };

  static DAZZ_DB db;

  char  *pwd, *root, *path;
  int    i;

  pwd = PathTo(arg);
  for (i = 0; i < 17; i++)
    { root  = Root(arg,suffix[i]);
      path  = Catenate(pwd,"/",root,extend[i]);
      fid   = fopen(path,"r");
      if (fid != NULL) break;
      free(root);
    }
  if (fid < 0)
    { fprintf(stderr,"\n%s: Cannot open %s as a .cram|[bs]am|f{ast}[aq][.gz]|db|dam file\n",
                     Prog_Name,arg);
      exit (1);
    }

  switch (i)
  { case 0:
      fclose(fid);
      return (NULL);
/*
      cid = cram_open(path);
      if (cid == NULL)
        {
        }
      return (Get_Cram_Entry);
*/

    case 1:
    case 2:
      fclose(fid);
      sid = sam_open(path);
      if (sid == NULL)
        { fprintf(stderr,"%s: Cannot open %s as a sam/bam file\n",Prog_Name,path);
          exit (1);
        }
      sam_header_process(sid,0);
      return (Get_Sam_Entry);

    case 3:
    case 4:
      fclose(fid);
      if (Open_DB(path,&db) < 0)
        { fprintf(stderr,"%s: Cannot open %s as a Dazzler DB or DAM\n",Prog_Name,path);
          exit (1);
        }
      return (Get_DB_Entry);

    case 5:
    case 7:
      read_line = Read_Line;
      return (&Get_Fastq_Entry);

    case 6:
    case 8:
      read_line = Read_Line;
      return (Get_Fasta_Entry);

    case 9:
    case 11:
    case 13:
    case 15:
      fclose(fid);
      read_line = Read_Zip_Line;
      gid = gzopen(path,"r");
      return (Get_Fastq_Entry);

    case 10:
    case 12:
    case 14:
    case 16:
      fclose(fid);
      gid = gzopen(path,"r");
      read_line = Read_Zip_Line;
      return (Get_Fasta_Entry);
  }

  return (NULL);
}
