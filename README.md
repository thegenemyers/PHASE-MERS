# Phasmer: A FastK K-mer phase chaining utility  
<font size ="4">**_Author:  Haynes Heaton, Richard Durbin, Gene Myers_**<br>
**_First:   April 5, 2021_**<br>

This is a repository under initial development so nothing is guaranteed to wokr or may be
broken at any time :-)

&nbsp;

  - [Phasemer](#phasemer): produce a table or listing of het-kmers and optional homo-kers
  
  - [Phasemap](#phasemap): maps the phase-mers to locations along a read data set 

&nbsp;


<a name="phasemer"></a>
```
1. Phasemer [-L[s]] [-h<int>:<int>] [-m<read:%> [-d<int>:<int>]]
                        [-N<path_name>] <source>[.ktab]
```

In a scan of the FastK table \<source>, Phasemer identifies all bundles of 2&#8209;4 k&#8209;mers that differ only in their middle base, i.e. the &lfloor;k/2&rfloor;<sup>th</sup> base.
These k&#8209;mers we call **het&#8209;mers** and all others we call **hom&#8209;mers**
(pronounced \`homers\`).

If the &#8209;h option is given then only k&#8209;mers whose count is in the given
range (inclusive) are considered legitimate het&#8209;mers, otherwise all are accepted.
If the &#8209;d option is given then only k&#8209;mers whose count is in the given
range (inclusive) are considered legitimate hom&#8209;mers, otherwise all are accepted.

Phasemer either produces tables or prints to stdout depending on the setting of
its' options as described below.

###Listing Options:

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
`272/g<` indicates the given k-mer comes from the 272'nd het-bundle and is the
variant that has a 'g' as its middle base.

###Table Building Options:

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
```
2. Phasemap [-D<self>[.prof]] <upper>[.prof] <lower>[.prof]
```

Given the upper and lower profiles built from the 2 tables output by Phasemer, Phasemap
scans the reads in order and outputs the locations and id's of all the het- and hom-mers
found in each read in order along the read.  (*Not complete*).


&nbsp;


Overall one does the following:

```
FastK -t? -p Foo               //  Build a table Foo.ktab and profile Foo.prof of HiFi data set Foo
Phasemer -m? Foo -NMers        //  Build tables Mers.U.ktab and Mers.L.ktab of het- and hom-mers
FastK -p:Mers.U Foo -NMap.U    //  Build relative profiles Map.U.prof of Foo w.r.t Mers.U
FastK -p:Mers.L Foo -NMap.L    //  Build relative profiles Map.L.prof of Foo w.r.t Mers.L
Phasemap -DFoo Map.U Map.L     //  Scan all profiles to output site locations in each read


