#!/bin/bash

# A simple tool that runs HGAP using 1.4 software.
# Must be run from the same directory as your input.xml
# Directory must be writeable.
# Must have your smrtpipe environment setup.

DEBUG=0
VERBOSE=0
SEYMOUR_HOME=${SEYMOUR_HOME:?"SMRT Pipe enviroment not detected."}
srcdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

usage() {
    echo -e "USAGE: $(basename $0) [params] <input.xml>\n"          \
            "-p    Optional path to a preassembler params file\n"   \
            "-r    Optional path to a resequencing params file\n"   \
            "-s    Optional path to a celera-assembler spec file\n"
    exit 1
}

debug() {
    if [ $DEBUG -eq 1 ]
    then
        echo $1 >&2
    fi
}

if [ $# -lt 1 -o "$1" == "--help" ]
then
    usage
fi

while getopts ":p:r:s:dv" opt
do
    case $opt in
      p)
        p_preasm="$OPTARG"
        ;;
      r)
        p_reseq="$OPTARG"
        ;;
      s)
        p_caspec="$OPTARG"
        ;;
      d)
        DEBUG=1
        ;;
      v)
        VERBOSE=1
        ;;
      \?)
        echo "Invalid option: -$OPTARG" >&2
        usage
        ;;
      :)
        echo "Option -$OPTARG requires an argument." >&2
        usage
        ;;
    esac
done

input=${@:$OPTIND:1}

p_preasm=${p_preasm="${srcdir}/params_preasm.xml"}
p_reseq=${p_reseq="${srcdir}/params_reseq.xml"}
caopts=

if [ \! -z $p_caspec ]
then
    caopts="-s $p_caspec"
fi

timeit() {
    desc=$1
    shift
    echo "Running $desc ..." >&2 
    debug "$@"
    
    cmd="$@"
    /usr/bin/time -f"${desc} time: %E real, %U user, %S sys" $cmd
    ec=$?
    if [ $ec -gt 0 ]
    then
        echo "$desc failed" >&2
        exit $ec
    fi 
}

debug "p_preasm = ${p_preasm}"
debug "p_reseq = ${p_reseq}"
debug "input = ${input}"
debug "caopts = ${caopts}"

if [ -z "$input" ]
then
    echo "INVALID cmd line: Must provide an input.xml."
    exit 1
fi

timeit "PreAssembler" "smrtpipe.py --distribute -D MAX_THREADS=60 -D HEARTBEAT_FREQ=-1 --params=${p_preasm} xml:${input}" > /dev/null
if [ \! -s data/corrected.fasta ] 
then
    echo "No corrected reads, consider lowering minLongReadLength in $p_preasm" >&2 
    exit 1
fi
timeit "CA prep" "fastqToCA -technology sanger -type sanger -reads data/corrected.fastq -libraryname reads" > reads.frg
timeit "CA" "runCA reads.frg -d assembly -p assembly $caopts" > ca.log
ln -s assembly/9-terminator/assembly.scf.fasta reference.fasta
timeit "Create Reference Repository" "referenceUploader -p. -f reference.fasta -c -n reference" > refupload.log
timeit "Resequencing" "smrtpipe.py --distribute -D MAX_THREADS=60 -D HEARTBEAT_FREQ=-1 --params=${p_reseq} xml:${input}" > /dev/null
