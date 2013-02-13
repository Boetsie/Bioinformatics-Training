#!/bin/bash

# A simple tool that runs HGAP using 1.4 software.
# Must be run from the same directory as your input.xml
# Directory must be writeable.
# Must have your smrtpipe environment setup.

QUEUE='secondary'
DEBUG=0
VERBOSE=0
SEYMOUR_HOME=${SEYMOUR_HOME:?"SMRT Pipe environment not detected."}
RunCA=${RunCA:-'runCA'}

srcdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
qsub=$(which qsub)

maybe_qsub=
if [ -x "$qsub" ]
then
    maybe_qsub="${qsub} -cwd -sync y -S /bin/bash -V -q ${QUEUE} -N CA -o ./CA.err -b y -j y -pe smp 15 "
fi

usage() {
    echo -e "USAGE: $(basename $0) [params] <input.xml>\n"          \
            "-p    Optional path to a preassembler params file\n"   \
            "-r    Optional path to a resequencing params file\n"   \
            "-s    Optional path to a celera-assembler spec file\n" \
            "-x    Override default options to smrtpipe\n"
    exit 1
}

debug() {
    if [ $DEBUG -eq 1 ]
    then
        echo $1 >&2
    fi
}

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

if [ $# -lt 1 -o "$1" == "--help" ]
then
    usage
fi

while getopts ":p:r:s:x:d" opt
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
      x)
        x_opts="$OPTARG"
        ;;
      d)
        DEBUG=1
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
x_opts=${x_opts="--distribute -D MAX_THREADS=60 -D HEARTBEAT_FREQ=-1"}
ca_opts=

if [ \! -z $p_caspec ]
then
    ca_opts="-s $p_caspec"
else
    if [ -e "${srcdir}/ca_default.spec" ]
    then 
        ca_opts="-s ${srcdir}/ca_default.spec"
    fi
fi

[ $DEBUG -eq 1 ] && x_opts="--debug ${x_opts}"

debug "p_preasm = ${p_preasm}"
debug "p_reseq = ${p_reseq}"
debug "input = ${input}"
debug "ca_opts = ${ca_opts}"
debug "x_opts = ${x_opts}"
debug "qsub = ${maybe_qsub}"
exit

if [ -z "$input" ]
then
    echo "INVALID cmd line: Must provide an input.xml."
    exit 1
fi

timeit "PreAssembler" "smrtpipe.py ${x_opts} --params=${p_preasm} xml:${input}" > /dev/null
if [ \! -s data/corrected.fasta ] 
then
    echo "No corrected reads, consider lowering minLongReadLength in $p_preasm" >&2 
    exit 1
fi
timeit "CA prep" "fastqToCA -technology sanger -type sanger -reads data/corrected.fastq -libraryname reads" > reads.frg
timeit "CA" "${maybe_qsub} ${RunCA} reads.frg -d assembly -p assembly ${ca_opts}"
ln -s assembly/9-terminator/assembly.scf.fasta reference.fasta
timeit "Create Reference Repository" "referenceUploader -p. -f reference.fasta -c -n reference"
timeit "Resequencing" "smrtpipe.py ${x_opts} --params=${p_reseq} xml:${input}" > /dev/null
