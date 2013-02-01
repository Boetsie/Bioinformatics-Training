HGAP 1.4 wrapper script
-----------------------
This is a convenience wrapper script to execute the HGAP assembly process using  
the 1.4 SMRT Pipe software release.  Its designed to execute the entire pipeline  
in one go with minimal input required from the user.

Files
-----
* hgap-1.4.sh - The main wrapper script
* params_preasm.xml - Default parameters for PreAssembler
* params_reseq.xml - Default parameters for Resequencing

Download these to a directory and include it in your PATH environment variable.   


Command Line
------------
Be sure to have your SMRT Pipe installation and environment ready prior to  
execution.  

    > hgap-1_4.sh --help
    USAGE: hgap-1_4.sh [params] <input.xml>
     -p    Optional path to a preassembler params file
     -r    Optional path to a resequencing params file
     -s    Optional path to a celera-assembler spec file

If the -p and/or -r options are provided, the default parameter files are   
overridden.  You must provide an input.xml which lists the base files to   
include in the analysis.

Output
------
The final output of HGAP will be located in the data/ directory:
* consensus.fasta.gz
* consensus.fastq.gz
