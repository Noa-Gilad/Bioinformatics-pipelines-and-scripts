#!/usr/bin/env bash
###############################################################################
# window2fbe.sh – scan Skipper windows for FBF-like motifs and bracket them
# Author: Noa Gilad
#
# Description:
#   Takes Skipper window-level TSV and:
#       1) Converts to BED and extracts sequences from a reference genome
#       2) Scans each sequence for FBF-like variants (8-mer and 7-mer motifs)
#       3) Adds "[" and "]" around motif occurrences in the FASTA output
#
# Usage:
#   ./window2fbe.sh windows.tsv[.gz] genome.fa output.fa [MID_MM]
#
#   MID_MM (optional):
#       - Maximum mismatches allowed in the variable positions of the FBF motif
#       - Default: 2
#
# Requirements:
#   - bedtools (getfasta)
#   - awk (with basic function support)
#   - zcat (for gzipped input)
#
# Output:
#   - FASTA file where sequences contain bracketed FBF-like motif hits.
###############################################################################
# Usage: window2fbe.sh windows.tsv[.gz] genome.fa output.fa [MID_MM]
# MID_MM = mismatches allowed ONLY in the variable positions (default 2)

set -euo pipefail
if [ $# -lt 3 ]; then
  echo "Usage: $0 <windows.tsv[.gz]> <genome.fa> <output.fa> [MID_MM]" >&2
  exit 1
fi

IN="$1"; GENOME="$2"; OUT="$3"; MID_MM="${4:-2}"

# TSV columns (1-based)
CHR_COL=1; START_COL=2; END_COL=3; STRAND_COL=6; GENE_COL=14

TMP_BED=$(mktemp)
TMP_TAB=$(mktemp)
trap 'rm -f "$TMP_BED" "$TMP_TAB"' EXIT

case "$IN" in *.gz) ZCMD="zcat" ;; *) ZCMD="cat" ;; esac

# 1) TSV -> BED6 (0-based start)
$ZCMD "$IN" | awk -F"\t" -v c=$CHR_COL -v s=$START_COL -v e=$END_COL -v st=$STRAND_COL -v g=$GENE_COL '
  NR==1 { next }
  {
    printf "%s\t%d\t%d\t%s|%s:%s-%s(%s)\t0\t%s\n",
           $c, $s-1, $e, $g, $c, $s, $e, $st, $st
  }' > "$TMP_BED"

# 2) Extract sequences
bedtools getfasta -fi "$GENOME" -bed "$TMP_BED" -s -name -tab > "$TMP_TAB"

# 3) Scan & bracket motifs
awk -v mid_mm="$MID_MM" '
function U(s){return toupper(s)}
function lowComp(s){return (s ~ /(A{6,}|T{6,}|C{6,}|G{6,})/)}

# count mismatches vs allowed sets at given positions
function miscnt(seq,p,pos,sets,   i,n,pp,ss,b,mm){
    mm=0
    n=split(pos,pp,","); split(sets,ss,",")
    for(i=1;i<=n;i++){
        b=U(substr(seq,p+pp[i]-1,1))
        if(index(ss[i],b)==0) mm++
    }
    return mm
}

function is8(seq,p,   mm){
    # TGT D H H A T
    if(U(substr(seq,p,3))!="TGT") return 0
    if(U(substr(seq,p+6,1))!="A" || U(substr(seq,p+7,1))!="T") return 0
    mm=miscnt(seq,p,"4,5,6","AGT,ACT,ACT")
    if(mm>mid_mm) return 0
    if(lowComp(substr(seq,p,8))) return 0
    return 8
}
function is7A(seq,p,  mm){
    # TGT D H H A
    if(U(substr(seq,p,3))!="TGT") return 0
    if(U(substr(seq,p+6,1))!="A") return 0
    mm=miscnt(seq,p,"4,5,6","AGT,ACT,ACT")
    if(mm>mid_mm) return 0
    if(lowComp(substr(seq,p,7))) return 0
    return 7
}
function is7AT(seq,p, mm){
    # TGT D H AT  (vars at 4,5)
    if(U(substr(seq,p,3))!="TGT") return 0
    if(U(substr(seq,p+5,1))!="A" || U(substr(seq,p+6,1))!="T") return 0
    mm=miscnt(seq,p,"4,5","AGT,ACT")
    if(mm>mid_mm) return 0
    if(lowComp(substr(seq,p,7))) return 0
    return 7
}

# collect ALL hits (allow overlaps)
function collect(seq,   i,L,len){
    n=0; len=length(seq)
    for(i=1;i<=len;i++){
        L=is8(seq,i);   if(L){n++; st[n]=i; ln[n]=L}
        L=is7A(seq,i);  if(L){n++; st[n]=i; ln[n]=L}
        L=is7AT(seq,i); if(L){n++; st[n]=i; ln[n]=L}
    }
    return n
}

# sort by start asc, length desc
function sort_hits(n,   i,j,tmp){
    for(i=1;i<=n;i++)
      for(j=i+1;j<=n;j++)
        if( st[j] < st[i] || (st[j]==st[i] && ln[j]>ln[i]) ){
            tmp=st[i]; st[i]=st[j]; st[j]=tmp;
            tmp=ln[i]; ln[i]=ln[j]; ln[j]=tmp;
        }
}

# reset arrays (mawk-safe)
function reset_marks(){
    for(i in o) delete o[i]
    for(i in c) delete c[i]
}

# place bracket markers
function mark(seq,n,   i,s,e){
    reset_marks()
    for(i=1;i<=n;i++){
        s=st[i]; e=s+ln[i]-1
        o[s]= (s in o ? o[s]"[" : "[")
        c[e]= (e in c ? c[e]"]" : "]")
    }
}

# rebuild string with markers
function rebuild(seq,   out,i,ch){
    out=""
    for(i=1;i<=length(seq);i++){
        if(i in o) out = out o[i]
        ch = substr(seq,i,1)
        out = out ch
        if(i in c) out = out c[i]
    }
    return out
}

{
    hdr=$1; seq=$2
    sub(/::.*/,"",hdr)
    delete st; delete ln

    n=collect(seq)
    if(n>1) sort_hits(n)
    if(n>0){ mark(seq,n); seq_out=rebuild(seq) } else seq_out=seq

    print ">" hdr
    print seq_out
}' "$TMP_TAB" > "$OUT"

echo "Done -> $OUT"
