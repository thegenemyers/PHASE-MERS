# Phasmer: A FastK K-mer phase chaining preprocessor  
<font size ="4">**_Author:  Haynes Heaton, Richard Durbin, Gene Myers_**<br>
**_First:   April 5, 2021_**<br>

This is a repository under initial development so nothing is guaranteed to work and may be
broken at any time :-)

&nbsp;

  - [Phasemer](#phasemer): produce a table or listing of het-kmers and optional homo-kers
  
  - [Phasemap](#phasemap): maps the phase-mers to locations along a read data set 
  
  - [An Example](#example): a complete example on a 10X data set of 8 gzipped fastq files. 

&nbsp;


<a name="phasemer"></a>
```
1. Phasemer [-L[s]] [-h<int>:<int>] [-m<read:%> [-d<int>:<int>]]
                        [-N<path_name>] <source>[.ktab]
```

In a scan of the FastK table \<source>, Phasemer identifies all bundles of 2&#8209;4 k&#8209;mers that differ only in their middle base, i.e. the &lfloor;k/2&rfloor;<sup>th</sup> base.
These k&#8209;mers we call **het&#8209;mers**, all others we call **hom&#8209;mers**
(pronounced \`homers\`), and any combination are referred to as **phase&#8209;mers**.

If the &#8209;h option is given then only k&#8209;mers whose count is in the given
range (inclusive) are considered legitimate het&#8209;mers, otherwise all are accepted.
If the &#8209;d option is given then only k&#8209;mers whose count is in the given
range (inclusive) are considered legitimate hom&#8209;mers, otherwise all are accepted.

Phasemer either produces tables or prints to stdout depending on the setting of
its' options as described below.

### Listing Options:

With the &#8209;L option set Phasemer sends a printout of each het-mer bundle to stdout (in no
particular order) along with the count of each.  For example,  

```
    ...
    
    1374344: cgataatcgacagtcaaaaaAgggtcaaacataaggcccc 23
             cgataatcgacagtcaaaaaGgggtcaaacataaggcccc 23

    1374345: cgataatcgcaaggattctaGcatttacgggtagtcgccc 27
             cgataatcgcaaggattctaTcatttacgggtagtcgccc 21

    1374346: cgataatcgcaccgattctaCcttttcccggtactgaccc 17
             cgataatcgcaccgattctaTcttttcccggtactgaccc 12

    1374347: cgataatcgcacggattctaCcatttcccggtactcaatc 12
             cgataatcgcacggattctaTcatttcccggtactcaatc 35

    1374348: cgataatcgcacggattctaAccatttacgggtagtcgcc 33
             cgataatcgcacggattctaCccatttacgggtagtcgcc 11

    1374349: cgataatcgcacggattctaAtttttcccggtactcgccc 11
             cgataatcgcacggattctaCtttttcccggtactcgccc 8

    ...
```
in response to `Phasemer -L -h8:40 CBS.ktab` where CBS is a 50X HiFi data set of
Cabernet Sauvignon.

With the &#8209;Ls option set, Phasemer produces an alphabetical listing of each het&#8209;mer.
In this order the het&#8209;mers are no longer necessarily organized into bundles as above.  So each
het&#8209;mer is preceded by the **id** assigned to the het-mer bundle and its variant base.
For example, in this listing:

```
        ...
        272/g: aaaaaaaaaaaaaaaaaaaaGtcaaaagtcatgaatagtg
        273/g: aaaaaaaaaaaaaaaaaaaaGtcaaatataattaaaatta
        277/g: aaaaaaaaaaaaaaaaaaaaGtcatagggtttctcatatc
        286/g: aaaaaaaaaaaaaaaaaaaaGtctactaagttgaaaacct
        289/g: aaaaaaaaaaaaaaaaaaaaGtctgctaagttgaaaacct
        294/g: aaaaaaaaaaaaaaaaaaaaGtgaaagatttgacttatgg
        297/g: aaaaaaaaaaaaaaaaaaaaGtgattggagaacacgtgtc
        307/g: aaaaaaaaaaaaaaaaaaaaGttgtgtttaggaactaaat
        309/g: aaaaaaaaaaaaaaaaaaaaGttttttttttttttttttt
          1/t: aaaaaaaaaaaaaaaaaaaaTaaaaaaaaaaataaaaata
         92/t: aaaaaaaaaaaaaaaaaaaaTaagcttctattctttcctt
        101/t: aaaaaaaaaaaaaaaaaaaaTacaaaaaaacaaaacaaaa
        104/t: aaaaaaaaaaaaaaaaaaaaTacaaagagcatttccatat
        125/t: aaaaaaaaaaaaaaaaaaaaTagaaaagaaaagaaaagaa
        128/t: aaaaaaaaaaaaaaaaaaaaTagaaatgatgacttgtatg
    ...
```
`272/g` indicates the given k-mer comes from the 272'nd het-bundle and is the
variant that has a 'g' as its middle base.

### Table Building Options:

When the &#8209;L option is not set, then Phasemer produces FastK tables as output, optionally using
the path name of the &#8209;N option as the root name for the table.  Ideally Phasemer would produce
a single table, say `HETS`, of relevant k-mers each associated with a 4 byte integer that encodes
a k&#8209;mer's type and id (more on this in a moment).  But FastK limits the associated values to 2&#8209;byte counts, so to work around this, Phasemer produces **two** tables, `HETS.U` and `HETS.L`
that contain exactly the same k-mers but the "U"&#8209;table counts give the upper 2 bytes of the
desired integer and the "L"&#8209;table counts give the lower 2 bytes.  By reading both tables in
sync or producing relative profiles for both and scanning both in sync, one can combine the
pair of 2&#8209;byte counts to construct the desired integer on the fly.

If the &#8209;N option is set, then the table names use the specified path as the root name.
If absent, then the root name of the source is used.  The option has no effect for the
list options involving &#8209;L.

With the &#8209;m option *not* set, Phasemer produces tables for the legitimate het&#8209;mers,
and if &#8209;m is set then it produces tables that contain both the legitimate het-mers
and &#8209;m percent of the available, legitimate hom&#8209;mers using the "modimizer" concept.

The integer associated with each k-mer gives a 29&#8209;bit id, and a 3&#8209;bit type.
Every het **bundle** gets an id and every modimizer-selected hom-mer individually gets an id,
where the id's begin at 0 and are consecutively allocated to bundles and hom&#8209;mers.  The
id is in the high order 29-bits of a k&#8209;-mers' integer.  The lower order 3 bits contain 0
if the k&#8209;mer is a hom&#8209;mer, and 1, 2, 3, or 4 if it is a het&#8209;mer, where the
codes correspond to a, c, g, t as the variant base, respectively.  In this way, without having
to consult the underly squence of the k&#8209;mer, one immediately knows which het&#8209;mer
of a bundle is being addressed.

&nbsp;

<a name="phasemap"></a>
```2. Phasemap [-10X] <lower>[.prof] <upper>[.prof]```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```<source>[.cram|.[bs]am|.db|.dam|.f[ast][aq][.gz] ...```


Given a pair of upper & lower phase-mer tables produced by Phasemer, say ```Mers.U``` & ```Mers.L```, one can then build relative upper & lower profiles of the phase-mers therein against *any* input data set, potentially in a series of individual files.  For example:

```
FastK -p:Mers.U -NMap.U X.fa.gz Y.fq Z.bam  // Build relative profiles Map.U of X,Y,Z w.r.t Mers.U
FastK -p:Mers.L -NMap.L X.fa.gz Y.fq Z.bam  // Build relative profiles Map.L of X,Y,Z w.r.t Mers.L
```

One can then give the profile pair and the input series profiled to Phasemap.  Continuing the
example:

```
Phasemap Map.L Map.U X.fa.gz Y.fq Z.bam  // Map the profile locations to X, Y, & Z
```

Phasemap outputs a series of lines to stdout that contain tab-separated values.  For each read,
a line that has 'R' in the first field is followed by a sequence of lines giving the location
of phase-mers found along that read in left-to-right order along the read.  The first field of
a hom-mer hit is 'H' and for het-mer hits it is 'A', 'C', 'G', or 'T' indicating the base found
at the variant position of the het-mer.  For both types of phase-mers the next three fields give
the position along the read, then the id of the hom-mer or het-mer bundle, and lastly a '+' or '-'
denoting the orientation of the k-mer as it is found in the read.

Often the input is paired 10X data implying that half the input reads have bar-codes and linkers.
If the -10x option is specified, then Phasemap ignores the first 23bp of a sequence, and outputs
the first 16bp in the third field of the read's 'R' line.  Otherwise a '-' is the third field.
The second field of a read line is always the ordinal index of the read within the input stream.

<a name="example"></a>
## An Example:

We start with a directory that contains 4 pairs of 10X gzip'd fastq files totalling 36.8Gbp
of data stored in 25.6GB:

```
> ls
Fish_L1_R1.fastq.gz    Fish_L3_R1.fastq.gz
Fish_L1_R2.fastq.gz    Fish_L3_R2.fastq.gz
Fish_L2_R1.fastq.gz    Fish_L4_R1.fastq.gz
Fish_L2_R2.fastq.gz    Fish_L4_R2.fastq.gz
```

Each 'R1' and 'R2' file is of paired reads where the 'R1' file contains the forward reads with
the 23bp bar-code and linker.  We want a table of the all 21-mers occuring 6 or more times,
**excluding** the bar-code and linker.  We do so as follows:

```
FastK -v -k21 -t6 -bc23 -T8 *R1.fastq.gz -NForward
FastK -v -k21 -t6 -T8 *R1.fastq.gz -NReverse
Logex -h -T8 'Fish=A|+B' Forward Reverse
Fastrm Forward Reverse
```

The commands above first build tables Forward and Reverse of the appropriate k-mers in the set
of forward and reverse reads.  Note carefully, that the option -bc23 in the first call to FastK
ensures that k-mers involving the bar-codes and linkers are not in the table ```Forward```.
Then Logex produces the union of the two tables in Fish along with a histogram, and the
final command removes the
tables Forward and Reverse as they are no longer needed.  On my laptop the above took 11:18 minutes
and the ```Fish``` table occupies 6.24GB.   At this point, we have:

```
> ls
Fish_L1_R1.fastq.gz    Fish_L3_R1.fastq.gz    Fish.hist
Fish_L1_R2.fastq.gz    Fish_L3_R2.fastq.gz    Fish.ktab
Fish_L2_R1.fastq.gz    Fish_L4_R1.fastq.gz
Fish_L2_R2.fastq.gz    Fish_L4_R2.fastq.gz
```

Now one can produce a pair of phase-mer tables, called ```Phase``` here.  This took 35 seconds
on my laptop and the two tables occupy only about 395 MB.

```
Phasemer -h6:40 -m2.5 -d6:100 -NPhase Fish
Fastrm Fish
> ls
Fish_L1_R1.fastq.gz    Fish_L3_R1.fastq.gz    Phase.L.ktab
Fish_L1_R2.fastq.gz    Fish_L3_R2.fastq.gz    Phase.U.ktab
Fish_L2_R1.fastq.gz    Fish_L4_R1.fastq.gz
Fish_L2_R2.fastq.gz    Fish_L4_R2.fastq.gz
```

The ```Phase``` tables can be used to map the phase-mers against *any* data set but for the
purposes of thise example, we do so for the original 10X data set.  This step is compute intensive
and took about 31 minutes on my laptop with the 4 profiles occupying 13.6GB of disk space.

```
FastK -v -k21 -p:Phase.L -T8 *R1.fastq.gz -NForward.L
FastK -v -k21 -p:Phase.U -T8 *R1.fastq.gz -NForward.U
FastK -v -k21 -p:Phase.L -T8 *R2.fastq.gz -NReverse.L
FastK -v -k21 -p:Phase.L -T8 *R2.fastq.gz -NReverse.U
> ls
Fish_L1_R1.fastq.gz    Fish_L3_R1.fastq.gz    Phase.L.ktab      Reverse.L.prof
Fish_L1_R2.fastq.gz    Fish_L3_R2.fastq.gz    Phase.U.ktab      Reverse.U.prof
Fish_L2_R1.fastq.gz    Fish_L4_R1.fastq.gz    Forward.L.prof
Fish_L2_R2.fastq.gz    Fish_L4_R2.fastq.gz    Forward.U.prof
```
Note we produced separate profile pairs for the forward and reverse reads, but this
time without ignoring the prefix of the forward reads so that Phasemap can output
barcodes.  In the last step we, then call Phasemap to produce maps for the forward
and reverse reads:

```
Phasemap -10X Forward.L Forward.U *R1.fastq.gz >Fmap.txt
Phasemap Reverse.L Reverse.U *R2.fastq.gz >Rmap.txt
Fastrm Forward.* Reverse.*
> ls
Fish_L1_R1.fastq.gz    Fish_L3_R1.fastq.gz    Phase.L.ktab
Fish_L1_R2.fastq.gz    Fish_L3_R2.fastq.gz    Phase.U.ktab
Fish_L2_R1.fastq.gz    Fish_L4_R1.fastq.gz    Fmap.txt
Fish_L2_R2.fastq.gz    Fish_L4_R2.fastq.gz    Rmap.txt
```
