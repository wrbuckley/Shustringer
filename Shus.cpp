//
// Shustringer.cpp : Defines the entry point for the console application.
//
// Generalised Byte-Oriented (Binary Octet) Shustring Computation Tool
//
//    This software computes 
//
//      the shortest subsequence that for any locus of a given symbol 
//      sequence, is found exactly once within that symbol sequence.
//
// (C) 2013 - 2015 by William R. Buckley
//                    California Evolution Institute
//                    wrb@calevinst.org
//                    All Rights Reserved.
//
//  Change Log
// =========================================================================================
//    date     who     what 
// ========= ======= =======================================================================
//  20151027   wrb   Adapted the GenomeQuadTree software, removed UDS and other code paths, 
//                     to produce the SHUStringer sequence analysis tool.
//  20151110   wrb   Completed translation of specific-case data structures (trees of nodes)
//                     and their build routines, to a node generalised for 256 possible 
//                     symbols and the corresponding set of build routines.
//  20151203   wrb   Date of first satisfactory behavior over wide range of datasets.
//
//   endlog

#include      "stdafx.h"

#include       <stdio.h>
#include      <stdlib.h>
#include      <string.h>
#include       <ctype.h>
#include        <time.h>
#include        <math.h>

#include      "getopt.h"

#define MAXSEQUENCELENGTH  100000000                                                                        // One Hundred Million entries
#define HSTFINALREPORT1    "================================\n\n"                                           //
#define HSTTITLE1          "\n   Shustring Length Histogram\n================================\n"            //
#define HSTTITLE2          "       length          count\n      --------      ---------\n"                  //
                                                                                                            //
enum machinestates {START, AGAIN, NEXT};                                                                    //
                                                                                                            //
typedef struct bytenode {                                                                                   //
       struct bytenode *    lnk;                                                                            //
       signed int           lst[256];                                                                       //
     unsigned int           cnt[256];                                                                       //
}                                                                                                           //
bytnode;                                                                                                    //
                                                                                                            //
    unsigned char    sequen[MAXSEQUENCELENGTH];                                                             // <     c            >    symbol - actually, a character (byte) converted to int by sign extension - characters are unsigned octets
      signed int     seqnxt[MAXSEQUENCELENGTH];                                                             // <  next            >    index of     next tile in candidate UDS
    unsigned int     seqcap[MAXSEQUENCELENGTH];                                                             // <  disp, cap, len  >    far end of shustring; length of shustring = seqcap[i] - i
                                                                                                            //
    unsigned int     hstgrm[MAXSEQUENCELENGTH];                                                             //
    unsigned int     hstrjct;                                                                               //
    unsigned int     hsttot;                                                                                //
                                                                                                            //
    FILE           * seq;                                                                                   // Source file of genomic data
    FILE           * rpt;                                                                                   // CSV file listing locus and length of each unique sequence
                                                                                                            //
    bytnode        * freelist;                                                                              //
    bytnode        * bytree;                                                                                //
    bytnode        * bytreelist;                                                                            //
                                                                                                            //                           
    time_t           starttime      = 0;                                                                    // time of entry-point execution                                                              Performance Metrics            
    time_t           loadtime       = 0;                                                                    //
    time_t           treetime       = 0;                                                                    //
    time_t           generatetime   = 0;                                                                    //
    time_t           endtime        = 0;                                                                    //
    time_t           shustart       = 0;                                                                    //
    time_t           shustop        = 0;                                                                    //
    time_t           histstart      = 0;                                                                    //
    time_t           histstop       = 0;                                                                    //
                                                                                                            //
      signed int     threshold      = 0;                                                                    // 
                                                                                                            //
    unsigned int     treealloc      = 0;                                                                    // 
    unsigned int     maxalloc       = 0;                                                                    // 
    unsigned int     notfound       = 0;                                                                    // Report on missing shustrings
    unsigned int     kindflag       = 0;                                                                    // 
    unsigned int     shustringflag  = 0;                                                                    // Generate shustrings for symbol sequence
    unsigned int     outfileflag    = 0;                                                                    // 
    unsigned int     histogramflag  = 0;                                                                    // Generate histogram for symbol sequence or shustrings
    unsigned int     symbollistflag = 0;                                                                    //
    unsigned int     thresholdflag  = 0;                                                                    //
    unsigned int     operationCode  = 0;                                                                    //
    unsigned int     sequenceLength = 0;                                                                    //
                                                                                                            //
             char    seqfname[257];                                                                         //
             char    rptfname[257];                                                                         //
             char    outfname[257];                                                                         //
                                                                                                            //
int arg_to_int(const char* arg, int min, int max, int defalt, int test){                                    // arg     - string to be converted
    int rv, i = defalt;                                                                                     // min/max - the minimum/maximum allowed value, inclusive
	                                                                                                        // defalt  - the default value, in case of an error
                                                                                                            // test    - flag to perform bounds checking on conversion result
    if(!arg) return defalt;                                                                                 // no argument means we use the default value
    rv = sscanf(arg, "%d", &i);                                                                             // make sure we got an integer argument
    if(rv != 1) return defalt;                                                                              // on scan error, return the default value
    if(test != 0){                                                                                          // on test request -
        if(i < min || max < i){                                                                             //   make sure the integer argument is within the desired range
            return defalt;}}                                                                                //   - return default otherwise
	return i;}                                                                                              // on all conditions satisfied, return converted result

void about(void){
    printf(
    "\n  Measurement of relative distinction between subsequences\n"
    "  found within a sequence of symbols is the primary measure\n"
    "  provided by the SHUStringer computation tool.  SHUStringer\n"
    "  compares all regions of the sequence to all other regions of\n"
    "  the sequence and computes the length of subsequence beginning\n"
    "  at each locus such that the subsequence is distinct form the\n"
    "  subsequences so computed for all other loci.  Hence, the length\n"
    "  measure serves as a crude indicator of relative distinction\n"
    "  between regions of a sequence of symbols.\n\n");}

void copyright(void){
	printf(
    "\nSHUStringer v1.0.1\n"
    "Copyright (C) 2013, 2014 & 2015 by William R. Buckley\n"
    "California Evolution Institute\n"
    "All Rights Reserved\n\n"
    );}

void help(void){
	printf(
    "Usage: shus <option list> <filename>\n\n"
    "   Operations:\n"
    " ----------------------------------------------------------------------------\n"
    "   -j                      generate histogram of shustring lengths\n"
 //   "   -n                      report on shustrings that are not found\n"
    "   -o <filename>           set name of output file\n"
 //   "   -s                      generate symbol list\n"
     "\n"
   "   Parameters:\n"
    " ----------------------------------------------------------------------------\n"
 //   "   -k                      kind of histogram - 0 for shustrings (default)\n"
 //   "                                               1 for symbols\n"
    "   -r                      replace existing output file\n"
    "   -t <count>              output shustrings longer(+) or shorter(-)\n"
    "                             than the <count> threshold\n"
        "\n"
    "   Commentaries: NB - use of these options brings an immediate exit from shus\n"
    " ----------------------------------------------------------------------------\n"
    "   -a                      about this software\n"
    "   -h                      help menu (this message list)\n"
    "   -l                      license information\n"
    //"   -m                      print a shus user's manual\n"
    "\n"
    "   Warnings!\n"
    " ----------------------------------------------------------------------------\n\n"
    "   Maximum Sequence Length is %d symbols\n\n"
    "   <count> is a positive integer value\n\n\n",MAXSEQUENCELENGTH
    );}

void license(void){
	printf(
    "The rights statement of William R. Buckley:\n"
    "===========================================\n"
    "THIS SOFTWARE IS A COPYRIGHTED WORK AND MAY NOT BE REDISTRIBUTED\n"
    "WITHOUT THE EXPRESS WRITTEN PERMISSION OF THE AUTHOR.\n\n\n"
    "Portions of this software are\n"
    "Copyright (c)2002-2003 Mark K. Kim\n"
    "All rights reserved.\n"
    "\n"
    "The rights statement of Mark K. Kim:\n"
    "------------------------------------\n"
    " Redistribution and use in source and binary forms, with or without\n"
    " modification, are permitted provided that the following conditions\n"
    " are met:\n"
    "\n"
    "   * Redistributions of source code must retain the above copyright\n"
    "     notice, this list of conditions and the following disclaimer.\n"
    "\n"
    "   * Redistributions in binary form must reproduce the above copyright\n"
    "     notice, this list of conditions and the following disclaimer in\n"
    "     the documentation and/or other materials provided with the\n"
    "     distribution.\n"
    "\n"
    "   * Neither the original author of this software nor the names of its\n"
    "     contributors may be used to endorse or promote products derived\n"
    "     from this software without specific prior written permission.\n"
    "\n"
    " THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n"
    " *AS IS* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n"
    " LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS\n"
    " FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE\n"
    " COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,\n"
    " INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,\n"
    " BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS\n"
    " OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED\n"
    " AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,\n"
    " OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF\n"
    " THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH\n"
    " DAMAGE.\n"
    );}
                                                                                                            //
void initbytnode(bytnode * r){                                                                              //
    unsigned int   i;                                                                                       //
                                                                                                            //
    r->lnk = 0;                                                                                             //
    for(i=0;i<256;i++){                                                                                     //
        r->lst[i] = -1;                                                                                     //
        r->cnt[i] = 0;}}                                                                                    //
                                                                                                            //
bytnode * mallocbytnode(void){                                                                              //
    bytnode * r;                                                                                            //
                                                                                                            //
    if((*freelist).lnk){                                                                                    //
        r = freelist;                                                                                       //
        freelist = (*freelist).lnk;                                                                         //
        initbytnode(r);                                                                                     //
        return(r);                                                                                          //
    }else{                                                                                                  //
        r = (bytnode *) malloc(sizeof(bytnode));                                                            //
        treealloc += sizeof(bytnode);                                                                       //
        initbytnode(r);                                                                                     //
        return(r);}}                                                                                        //
                                                                                                            //
void freebytnode(bytnode * r){                                                                              // r - a free bytnode : put the bytnode in the LIFO free list - link the newly free 
    r->lnk = (*freelist).lnk;                                                                               //   node to the list of other free nodes, by pointing at the list
    (*freelist).lnk = r;}                                                                                   // adjust the head pointer of the free list - add a new head to the list
                                                                                                            //
void buildbytroot(bytnode * r){                                                                             // r - an initialised, empty bytnode
    unsigned int    i;                                                                                      // i - an index on the sequence
    unsigned char   c;                                                                                      // c - the i'th character of the sequence
                                                                                                            //
    initbytnode(r);                                                                                         // initialise the bytnode r                                                                   Redistribute the indicies of the elements of a list based upon the next character of the element.
    fseek(seq,0L,SEEK_SET);                                                                                 // seek to beginning of sequence
    for(i=0;i<MAXSEQUENCELENGTH;i++){                                                                       //
        seqnxt[i] = -1;}                                                                                    //
    i = 0;                                                                                                  //
    while(1){                                                                                               // while symbols remain to be read from the sequence
        c = fgetc(seq);                                                                                     // get the next symbol of the sequence
        if(feof(seq)){                                                                                      // NB - all files (in the real world) are of finite size - hence, the forever loop must
            sequenceLength = i;                                                                             //   eventually execute this block:   note the sequence length
            return;}                                                                                        //                                    exit the routine
        seqnxt[i] = r->lst[(unsigned int)c];                                                                // build a list of backward pointers, and so
        r->lst[(unsigned int)c] = i;                                                                        //   group the symbols by character
        r->cnt[(unsigned int)c]++;                                                                          // increment the character count
        seqcap[i] = i;                                                                                      // set initial length of cap to zero
        sequen[i] = c;                                                                                      // add character to sequence
        i++;}}                                                                                              // keep track of the sequence length (count of symbols in sequence)
                                                                                                            //
void buildbytnode(bytnode * r, unsigned int t){                                                             // r    - an as-yet-unused bytnode                                                            Sort a subsequence t into node r - the list of elements of the subsequence are linked via seqnxt[].
                                                                                                            // t    - is the index of a sortable subsequence                                              The elements of a subsequence generally are discontiguous with respect to the locus in the sequence
    unsigned int   hold;                                                                                    // hold - keeps track of the next index in the sortable subsequence!                            at which these elements are found.
    unsigned int   s;                                                                                       // s    - surrogate, to make code more readable (a sign-extended UINT32)
                                                                                                            //
    initbytnode(r);                                                                                         // initialise the bytnode r                                                                   Redistribute the indicies of the elements of a list based upon the next character of the element.
    while(t != -1){                                                                                         //
        if((seqcap[t] - t) < sequenceLength){                                                               // no cap can exceed the length of the entire sequence
            hold = seqnxt[t];                                                                               // hold the index to the next element of the subsequence
            seqcap[t]++;                                                                                    // increase the cap length of subsequence t
            if(seqcap[t] < sequenceLength){                                                                 // if the cap length of subsequence t is less than the length of the entire sequence          The cap is a length plus the index of the cap.  It may easily in value exceed the sequenceLength, and
                s = (unsigned int) sequen[seqcap[t]];                                                       //   then we need make no adjustment to that cap length in order to use it as an index          this fact must be taken into account when conducting the sort.
                seqnxt[t] = r->lst[s];                                                                      // add the subsequence to the list that corresponds to the symbol that is found at the
                r->lst[s] = t;                                                                              //   end of the subsequence - via the element found at the cap of that subsequence
                r->cnt[s]++;                                                                                // keep track of how many subsequences are in that s list
            }else{                                                                                          //
                s = (unsigned int) sequen[(seqcap[t] - sequenceLength)];                                    // if the cap length of the subsequence t is greater than or equal to the length of the       This is the case where the cap plus the index of the cap (the sum) has a value that is larger than is 
                seqnxt[t] = r->lst[s];                                                                      //   entire sequence, then we need to make an adjustment to the cap length in order to          the value of the sequenceLength.
                r->lst[s] = t;                                                                              //   use it as an index; add the subsequence to the list that corresponds to the symbol
                r->cnt[s]++;}                                                                               //   that is found at the end of the subsequence - indicated by the element found at the
            t = hold;                                                                                       //   cap of that subsequence; keep track of how many subsequences are in that s list
        }else{                                                                                              //
            t = seqnxt[t];}}}                                                                               //
                                                                                                            //
void buildbytree(bytnode * r,unsigned int level){                                                           // r     - a previously-populated bytnode                                                     This routine is constructed in state-machine format, and is recursively called.
                                                                                                            // level - the depth of the current node within the suffix tree
    unsigned int   i;                                                                                       // i     - general index, ranges over all 256 byte values                                                                                            
      signed int   nzi  = -1;       /* not zero index   */                                                  // nzi   - index to the list of the last symbol having at least one element in that list        NB - generate here, in this routine, reports respecting missing shustrings.   
    unsigned int   nzc  =  0;       /* not zero count   */                                                  // nzc   - number of lists that have at least one element in that list                                
    unsigned int   bbts =  START;   /* buildbytreestate */                                                  // bbts  - the state of the machine                                                                   
    bytnode      * s;                                                                                       // s     - pointer to new bytnode - facilitates recursive call
                                                                                                            //
    while(1){                                                                                               // while not done - being done is equivalent to reaching the NEXT state, when only one        
        switch(bbts){                                                                                       //   element is in the only list that has any elements, or by reaching an error state         
            case START: for(i=0;i<256;i++){                                                                 // for each symbol in the ASCII character set - each byte value, [0..255]                     Determine the number of characters for which the corresponding symbol_list is not of length zero.
                            if(r->cnt[i] > 1){                                                              //   count up the number of symbols for which the list_count is more than one
                                nzc++;                                                                      //
                                nzi = i;}}                                                                  //   and keep track of which symbol last reported a list_count of more than one
                        if(nzc == 1){                                                                       // if only one symbol has a not_zero list_count, then this node can be re-used,               Detect and service the special state of all symbols appearing in one list
                            bbts = AGAIN;                                                                   //   else, for the case of only one symbol in one list, we are done (nothing to               For AGAIN state, no need to allocate new node; simply resort the list on the current node.
                        }else{                                                                              //   otherwise,                                                                               For the case of more than one list containing substrings,                                                            
                            if(nzc > 1){                                                                    //   where more than one bucket has a list_count greater than one,                              and given a reasonable state for the node (some number of nodes greater than one have non-empty lists),
                                bbts = NEXT;                                                                //   a new node must be built                                                                 we need to sort each non_empty list upon a new node, at depth node_depth + 1.
                            }else{                                                                          //   and, where no bucket has a list_count greater than one,                                  
                                return;}}                                                                   //   we are done with this node                                                                
                        break;                                                                              //
                                                                                                            //
            case AGAIN: buildbytnode(r, r->lst[nzi]);                                                       // we re-sort the node without recursion when only one list of the node has elements          Re-sorting a node occurs when buildbytnode() is called a second or subsequent time for any one invocation 
                        bbts = START;                                                                       //                                                                                              of buildbytree() - called within that single invocation.
                        nzc = 0;                                                                            // reset the counter to zero, so as to facilitate generation of a correct, new total
                        break;                                                                              //
                                                                                                            //
            case  NEXT: s = mallocbytnode();                                                                // when there are at least two nodes with elements (even if they are both single 
                        for(i=0;i<256;i++){                                                                 //   elements) then we need to sort each list in a new node.  So, for each character 
                            if(r->cnt[i] > 1){                                                              //   amongst the elements, we see if the corresponding count is greater than one, and 
                                buildbytnode(s, r->lst[i]);                                                 //   if so, we build a new node (buildbytnode) and then complete the tree (buildbytree).
                                buildbytree(s, level + 1);}}                                                //
                        freebytnode(s);                                                                     // after the node is serviced, we may free it, as it will otherwise serve no useful 
                        return;                                                                             //   purpose.
                                                                                                            //
            default:    break;}}}                                                                           // report error condition - we should never get to the default case.
                                                                                                            //
int main( int argc, char *argv[], char **envp){                                                             //
    unsigned int   t,tt;                                                                                    //
    unsigned int   i;                                                                                       //
                                                                                                            //
    starttime = clock();                                                                                    //
    freelist = (bytnode *) malloc(sizeof(bytenode));                                                        //
    initbytnode(freelist);                                                                                  //
    copyright();                                                                                            //
    for(i=0;i<257;i++) *(seqfname+i) = *(rptfname+i) = *(outfname+i) = '\0';                                //
    while(1){                                                                                               //
        int c = getopt(argc, argv, "-ahjklmnsf:o:t:");                                                      //
        if(c == -1) break;                                                                                  //
		switch(c){                                                                                          //
            case 'a':                                                                                       //
                about();                                                                                    //
                exit(0);                                                                                    //
			case 'h':                                                                                       //
                help();                                                                                     //
                exit(0);                                                                                    //
            case 'j':                                                                                       //
                histogramflag = 1;                                                                          //
                break;                                                                                      //
            case 'k':                                                                                       //
                kindflag = arg_to_int(optarg, 0, 0, 0, 0);                                                  //
                if((kindflag < 0) || (kindflag > 1)){                                                       //
                    printf("Invalid histogram type code\n");                                                //
                    exit(0);}                                                                               //
                break;                                                                                      //
            case 'l':                                                                                       //
                license();                                                                                  //
                exit(0);                                                                                    //
            case 'm':                                                                                       //
                exit(0);                                                                                    //
                break;                                                                                      //
            case 'n':                                                                                       //
                notfound = 1;                                                                               //
                break;                                                                                      //
			case 'o':                                                                                       //
                strcpy(outfname,optarg);                                                                    //
                break;                                                                                      //
            case 's':                                                                                       // Shuffle a Uniformly Disordered Sequence
                symbollistflag = 1;                                                                         //
                break;                                                                                      //
            case 't':                                                                                       //
                threshold = arg_to_int(optarg, 0, 0, 0, 0);                                                 //
                thresholdflag = 1;                                                                          //
                break;                                                                                      //
			case   1:                                                                                       //
                strcpy(seqfname,optarg);                                                                    //
                break;}}                                                                                    //
    if(*seqfname){                                                                                          //
        if(0 != (seq = fopen(seqfname,"rb"))){                                                              //
            bytree = (bytnode *)malloc(sizeof(bytnode));                                                    //
            buildbytroot(bytree);                                                                           // Initialise the root of the "suffix tree" analog
            fclose(seq);                                                                                    //
            loadtime = clock();                                                                             //
            buildbytree(bytree,1);                                                                          // Build the "suffix tree" analog
            treetime = clock();                                                                             //
            if(histogramflag == 1){                                                                         // This histogram service is for post processing of generated UDS, Permutation and Shustrings.
                hsttot = 0;                                                                                 //
                for(t=0;t<MAXSEQUENCELENGTH;t++){                                                           // Ideally, this is governed by cardiality of the symbol set; i.e., ||S||
                    hstgrm[t] = 0;}                                                                         //
                hstrjct = 0;                                                                                //
                printf(HSTTITLE1);                                                                          //
                printf(HSTTITLE2);                                                                          //
                for(t=0;t<sequenceLength;t++){                                                              //
                    tt = (seqcap[t] - t) + 1;                                                               //
                    if(tt < MAXSEQUENCELENGTH){                                                             //
                        hstgrm[tt]++;                                                                       //
                    }else{                                                                                  //
                        hstrjct++;}}                                                                        //
                for(t=0;t<MAXSEQUENCELENGTH;t++){                                                           // Print generated histogram
                    if(hstgrm[t] != 0){                                                                     //
                        printf("  %12d   %12d\n",t,hstgrm[t]);                                              //
                        hsttot += hstgrm[t];}}                                                              //
                printf(HSTFINALREPORT1);                                                                    //
            }else{                                                                                          //
                if(*outfname){                                                                              //
                    if(0 != (rpt = fopen(outfname,"w"))){                                                   //
                        for(t=0;t<sequenceLength;t++){                                                      //
                            fprintf(rpt,"%12d,%12d,%c\n",t,(1 + (seqcap[t] - t)),sequen[t]);}               //
                        fclose(rpt);                                                                        //
                    }else{                                                                                  //
                        printf("Error opening output file:.%s.\n",rptfname);}                               //
                }else{                                                                                      //
                    for(t=0;t<sequenceLength;t++){                                                          //
                        i = (seqcap[t] - t) + 1;                                                            //
                        if(threshold){                                                                      //
                            if(threshold < (signed int)i){                                                  //
                                printf("%12d,%12d,%c\n",t,i,sequen[t]);}                                    //
                        }else{                                                                              //
                            printf("%12d,%12d,%c\n",t,i,sequen[t]);}}                                       //
                    printf("\n\n");}}                                                                       //
        }else{                                                                                              //
            printf("Unable to open input file:%s\n",seqfname);}                                             //
    }else{                                                                                                  //
        printf("Name of input file not given.\n");}                                                         //
    generatetime = clock();                                                                                 //
    endtime = clock();                                                                                      //
    printf(" Sequence Length:%12d\n",sequenceLength);                                                       //
    printf("     Start  time:%12d\n",starttime);                                                            //
    printf("     Load   time:%12d\n",loadtime);                                                             //
    printf("     Tree   time:%12d\n",treetime);                                                             //
    printf("     Report time:%12d\n",generatetime);                                                         //
    printf("     Stop   time:%12d\n\n",endtime);                                                            //
    printf("       treealloc:%12ld\n",treealloc);                                                           //
    printf("sizeof(treenode):%12d\n\n\n",sizeof(bytnode));}                                                 //
                        
 