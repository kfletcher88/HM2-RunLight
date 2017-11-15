#!/bin/bash
#Wrapper to run HaploMerger2 batch jobs hm.batchA1-3, hm.batchB1-5, hm.batchD1-3

#Please refer to the original HaploMerger2 documentation for instructions on how to run additional scaffolding (hm.batch C1-2) and gap-filling (hm.batch E1) jobs.
#./manual.v3.4.pdf/manual.v3.4.pdf
#/github.com/mapleforest/HaploMerger2/blob/master/manual.v3.4.pdf

#MIT liscense
#Copyright 2017
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction,
#including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
#subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
#IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
#OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

usage="($basename "$0") [-f] [-t] [-a] [-p] [-h] [-n] [-N] -- A Shell wrapper to run HaploMerger2 batch jobs hm.batchA1-3, hm.batchB1-5, hm.batchD1-3
	Additional scaffolding and gap-filling steps may be run as part of the original pipeline, available at:
	https://github.com/mapleforest/HaploMerger2

Options:
	-h show this help message/

Required Arguments:
	-f An input assembly file, preferably soft-masked and gzipped (critical).
	-t Number of threads for multi-threaded stages to use (default = 1).
	-a Instructs the script to run the final steps, (hm.batchD1-3; remove tandem assemblies) on the alternative assembly, as well as the reference assembly.
	-p Prefix for outputs, deposited in directory HM2results (default = HM2out).
	-n Do not clean up the intermediate files produced by HaploMerger2. This can add up to Gigabytes of disk space [enforces N].
	-N Do not clean up intermediate assembly files, but remove everything else.
Outputs:
	Two assemblies.
	1. [InputAssemblyName]_mp_ref_rt.fa.gz
	2. Without specifying the -a flag:
		[InputAssemblyName]_mp_alt.fa.gz
	   Specifying the -a flag:
		[InputAssemblyName]_mp_alt_rt.fa.gz
"

while getopts ':hanNt:f:p:' option; do
        case "$option" in
                h)  echo "$usage"
                         exit
                         ;;
		f)  Fasta=$OPTARG
			 ;;
		t)  Threads=$OPTARG
			 ;;
		p)  Prefix=$OPTARG
			 ;;
		a)  Alt=1
			 ;;
		n)  NoClean1=1
			 ;;
		N)  NoClean2=1
			 ;;
		\?)printf "illegal option: -%s\n\n" "$OPTARG" >&2
                    echo "$usage"
                    exit 1
                         ;;
        esac
done

#Check critical arguments are provided.
if [[ -z $Fasta ]];then
       echo "
       ERROR: You must provide a FASTA file
"
       echo "$usage"
       exit 1
fi

#Set some defaults if not provided
if [[ -z $Threads ]]; then
Threads=1
echo "
Threads not specified, will proceed with default [1]
"
fi

if [[ -z $Prefix ]]; then
Prefix=HM2out
echo "No Prefix provided, outputs will be prefixed HMout in HM2results/
"
else
echo "Results will be deposited in HM2results with prefix of $Prefix
"
fi
#Output message on cleanliness

if [[ -n $NoClean1 ]]; then
echo "
All files output by Haplomerge2 will be retained
"
fi

if [[ -n $NoClean2 ]]; then
echo "
All intermediate assembly files output by HaploMerger2 will be retained
"
fi

#SetPath
exe_bin=$(dirname "$0")
#Check executable are there
HMexeTest=$(command -v $exe_bin/mapleforest/bin/initiation.pl)
if [[ $HMexeTest == "" ]]; then
	echo "The Haplomerger 2 bin is not accessible. Was the script moved? The bin should be located at ./mapleforest/bin/"
	exit 1
else
	echo "All looks good, let's proceed"
fi

#PATH=$exe_bin/mapleforest/bin/:$PATH
#echo "Set PATH to include $exe_bin/mapleforest/bin/ and $exe_bin/mapleforest/HMbatches"

#Set-up workspace
mkdir -p HM2temp

#Check if fasta file is gzipped.
if [[ $Fasta =~ .gz$ ]]; then
echo "
It looks like $Fasta is gzipped."
ScaffC=$(zcat $Fasta | grep -c '^>')
	if [[ $ScaffC == 0 ]]; then
	echo "No scaffolds detected, is $Fasta a gzipped fasta file? We cannot proceed"
	exit 1
	else
	echo "$Fasta contains $ScaffC scaffolds"
	cp $Fasta HM2temp/Input.fa.gz
	fi
else
echo "
It looks like $Fasta is not gzipped."
ScaffC=$(cat $Fasta | grep -c '^>')
        if [[ $ScaffC == 0 ]]; then
        echo "No scaffolds detected, is $Fasta a gzipped fasta file? We cannot proceed"
        exit 1
        else
        echo "$Fasta contains $ScaffC scaffolds"
        cp $Fasta HM2temp/Input.fa
	gzip HM2temp/Input.fa
	fi
fi
cp -R $exe_bin/mapleforest/bin HM2temp
cp -R $exe_bin/mapleforest/HMbatches/ HM2temp
cp $exe_bin/mapleforest/HMbatches/all_lastz.ctl $exe_bin/mapleforest/HMbatches/scoreMatrix.q HM2temp/
cp -R $exe_bin/mapleforest/chainNet_jksrc20100603_ubuntu64bit HM2temp
mkdir -p HM2results
cd HM2temp
PATH=./bin/:./HMbatches/:./chainNet_jksrc20100603_ubuntu64bit:$PATH


#Begin Haplomerger2
echo "running Haplomerger batch A modules"
echo "setting file open limit to 655350"
ulimit -n 655350
hm.batchA1.initiation_and_all_lastz Input $Threads
hm.batchA2.chainNet_and_netToMaf Input $Threads
hm.batchA3.misjoin_processing Input
hm.batchB1.initiation_and_all_lastz Input_mp $Threads
hm.batchB2.chainNet_and_netToMaf Input_mp $Threads
hm.batchB3.haplomerger Input_mp
hm.batchB4.refine_unpaired_sequences Input_mp $Threads
hm.batchB5.merge_paired_and_unpaired_sequences Input_mp
hm.batchD1.initiation_and_all_lastz Input_mp_ref $Threads
hm.batchD2.chainNet_and_netToMaf Input_mp_ref $Threads
hm.batchD3.remove_tandem_assemblies Input_mp_ref
mv Input_mp_ref_rt.fa.gz ../HM2results/$Prefix.mp_ref_rt.fa.gz
#Run batchD on alternative if flag a provided
if [[ $Alt == 1 ]]; then
hm.batchD1.initiation_and_all_lastz Input_mp_alt $Threads
hm.batchD2.chainNet_and_netToMaf Input_mp_alt $Threads
hm.batchD3.remove_tandem_assemblies Input_mp_alt
mv Input_mp_alt_rt.fa.gz ../HM2results/$Prefix.mp_alt_rt.fa.gz
else
mv Input_mp_alt.fa.gz ../HM2results/$Prefix.mp_alt.fa.gz
echo "

the final scripts were not run on the alternative assembly. The alternative assembly, prior to tandemn removal has been deposited in HM2results/
"
fi

cd ../
if [[ -z $NoClean1 && -z $NoClean2 ]]; then
rm HM2temp
exit
else
	if [[ -n $NoClean1 ]]; then
	exit
	fi
	if [[ -n $NoClean2 ]]; then
	echo "Intermediate assemblies will be deposited in HM2intermediates
"
	mkdir -p HM2intermediates
		if [[ -n $Alt ]]; then
		mv HM2temp/Input_mp_alt.fa.gz HM2intermediates/$Prefix.mp_alt.fa.gz
		fi
	mv HM2temp/Input_mp.fa.gz  HM2intermediates/$Prefix.mp.fa.gz
	mv HM2temp/Input_mp_ref.fa.gz HM2intermediates/$Prefix.mp_ref.fa.gz
	mv HM2temp/Input_mp_refx.fa.gz HM2intermediates/$Prefix.mp_refx.fa.gz
	mv HM2temp/Input_mpx.fa.gz HM2intermediates/$Prefix.mpx.fa.gz
	mv HM2temp/Inputx.fa.gz HM2intermediates/$Prefix.x.fa.gz
	rm HM2temp -R
	exit
	fi
fi
