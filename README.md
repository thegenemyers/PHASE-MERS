# Phaser: A K-mer phase chain algorithm  
<font size ="4">**_Author:  Haynes Heaton, Richard Durbin, Gene Myers_**<br>
**_First:   April 5, 2021_**<br>

  - [Phasemer](#phasemer): produce a table or listing of het-kmers and optional homo-kmers

&nbsp;


<a name="phasemer"></a>
```
Phasemer [-L[s]] [-h<int>:<int>] [-m<read:%> [-d<int>,<int>]] <source>[.ktab]
```

In a scan of the FastK table \<source>, Phasemer identifies all bundles of 2&#8209;4 k&#8209;mers that differ only in their middle base, i.e. the &lfloor;k/2&rfloor;<sup>th</sup> base.
These k&#8209;mers we call **het&#8209;mers** and all others we call **hom&#8209;mers**.

If the &#8209;h option is given then only k&#8209;mers whose count is in the given
range (inclusive) are considered legitimate het&#8209;mers, otherwise all are accepted.
If the &#8209;d option is given then only k&#8209;mers whose count is in the given
range (inclusive) are considered legitimate hom&#8209;mer, otherwise all are accepted.

Phasemer either produces tables or pints to stdout depending on the setting of
its' options as desscribed below.

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
in response to <code>Phasemer -L -h8:40 CBS.ktab</code> where CBS is a 50X HiFi data set of
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
<code>272/g</code> indicates the given k-mer comes from the 272'nd het-bundle and is the
variant that has a 'g' as its middle base.

###Table Building Options:

When a -L option is not set, then Phasemer produces 2 tables, for now 

