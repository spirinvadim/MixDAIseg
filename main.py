import argparse
import subprocess

parser = argparse.ArgumentParser(description='DAIseg')

subparser = parser.add_subparsers(dest='mode', required=True)

parser_daiseg = subparser.add_parser("daiseg")
parser_daiseg.add_argument('--bed', type=str, help='Region bed file')
parser_daiseg.add_argument('--EM', type=str, help='Whether or not to use EM algorithm')
parser_daiseg.add_argument('--EM_steps', type=str, help='number of EMsteps')
parser_daiseg.add_argument('--EM_samples', type=int, help='Number of samples used in training')
parser_daiseg.add_argument('--HMM_par', type= str, help='File with parameters')
parser_daiseg.add_argument('--o', type= str, help = 'Name of output file' )
parser_daiseg.add_argument('--EM_est', type= str, help = 'Make estimation of the all parameters or only coalescent times' )
parser_daiseg.add_argument('--prepared_file', type=str, help='dlkfjgljk')
parser_daiseg.add_argument('--arch_cover', type=str)
parser_daiseg.add_argument('--obs_samples', type=str, help='File with samples names')
parser_daiseg.add_argument('--decoding', type=str, help='Viterbi or aposteriory decoding')
parser_daiseg.add_argument('--cut_off', type=float, help='Decoding cut off')

parser_mex_daiseg = subparser.add_parser("mex_daiseg")
parser_mex_daiseg.add_argument('--bed', type=str, help='Region bed file')
parser_mex_daiseg.add_argument('--EM', type=str, help='Whether or not to use EM algorithm')
parser_mex_daiseg.add_argument('--EM_steps', type=str, help='number of EMsteps')
parser_mex_daiseg.add_argument('--EM_samples', type=int, help='number of samples used in EM algorithm')
parser_mex_daiseg.add_argument('--HMM_par', type= str, help='File with parameters')
parser_mex_daiseg.add_argument('--out_prefix', type= str, help = 'Prefix to output file' )
parser_mex_daiseg.add_argument('--prepared_file', type=str, help='dlkfjgljk')
parser_mex_daiseg.add_argument('--arch_cover', type=str)
parser_mex_daiseg.add_argument('--obs_samples', type=str, help='File with samples names')
parser_mex_daiseg.add_argument('--obs_type', type=str, help='Type of observaions simple/independent')
parser_mex_daiseg.add_argument('--transition_matrix', type=str, help='Chose the type of transition matrix')

parser_leo_daiseg = subparser.add_parser("leo_daiseg")
parser_leo_daiseg.add_argument("-ind",help="[required] ingroup/outgrop list (json file) or comma-separated list e.g. ind1,ind2", type=str, required = True)
parser_leo_daiseg.add_argument("-vcfOut",help="[required] path to list of comma-separated vcf/bcf file(s) containing the individuals in the outgroups", type=str, required = True)
parser_leo_daiseg.add_argument("-vcfIn",help="[required] path to list of comma-separated vcf/bcf file(s) containing the individuals in the ingroup", type=str, required = True)
parser_leo_daiseg.add_argument("-demo",help="demographic model file", type=str, required = True)
parser_leo_daiseg.add_argument("-weights", metavar='',help="file with callability (defaults to all positions being called)")
parser_leo_daiseg.add_argument("-out", metavar='',help="outputfile (defaults to stdout)", default = '/dev/stdout')
parser_leo_daiseg.add_argument("-ancestral", metavar='',help="fasta file with ancestral information - comma-separated list or wildcards like vcf argument (default none)", default='')
parser_leo_daiseg.add_argument("-refgenome", metavar='',help="fasta file with reference genome - comma-separated list or wildcards like vcf argument (default none)", default='')
parser_leo_daiseg.add_argument("-haploid",help="Change from using diploid data to haploid data (default is diploid)", action='store_true', default = False)
parser_leo_daiseg.add_argument("-conditional",help="Use conditional probability (default is False)", action='store_true', default = False)

args = parser.parse_args()

data = ['python3']
if args.mode == 'daiseg':
	data.append('src/daiseg/dai.seg.py')
	if args.bed:
		data.append('--bed')
		data.append(args.bed)
	if args.EM:
		data.append('--EM')
		data.append(args.EM)
	if args.EM_steps:
		data.append('--EM_steps')
		data.append(args.EM_steps)
	if args.EM_samples:
		data.append('--EM_samples')
		data.append(args.EM_samples)
	if args.HMM_par:
		data.append('--HMM_par')
		data.append(args.HMM_par)
	if args.o:
		data.append('--o')
		data.append(args.o)
	if args.EM_est:
		data.append('--EM_est')
		data.append(args.EM_est)
	if args.prepared_file:
		data.append('--prepared_file')
		data.append(args.prepared_file)
	if args.arch_cover:
		data.append('--arch_cover')
		data.append(args.arch_cover)
	if args.obs_samples:
		data.append('--obs_samples')
		data.append(args.obs_samples)
	if args.decoding:
		data.append('--decoding')
		data.append(args.decoding)
	if args.cut_off:
		data.append('--cut_off')
		data.append(args.cut_off)
elif args.mode == 'mex_daiseg':
	data.append('src/mex_daiseg/daiseg.mex.2.py')
	if args.bed:
		data.append('--bed')
		data.append(args.bed)
	if args.EM:
		data.append('--EM')
		data.append(args.EM)
	if args.EM_steps:
		data.append('--EM_steps')
		data.append(args.EM_steps)
	if args.EM_samples:
		data.append('--EM_samples')
		data.append(args.EM_samples)
	if args.HMM_par:
		data.append('--HMM_par')
		data.append(args.HMM_par)
	if args.out_prefix:
		data.append('--out_prefix')
		data.append(args.out_prefix)
	if args.prepared_file:
		data.append('--prepared_file')
		data.append(args.prepared_file)
	if args.arch_cover:
		data.append('--arch_cover')
		data.append(args.arch_cover)
	if args.obs_samples:
		data.append('--obs_samples')
		data.append(args.obs_samples)
	if args.obs_type:
		data.append('--obs_type')
		data.append(args.obs_type)
	if args.transition_matrix:
		data.append('--transition_matrix')
		data.append(args.transition_matrix)
elif args.mode == 'leo_daiseg':
	data.append('src/leo_daiseg/src/main.py')
	if args.ind:
		data.append('-ind=' + args.ind)
	if args.vcfOut:
		data.append('-vcfOut=' + args.vcfOut)
	if args.vcfIn:
		data.append('-vcfIn=' + args.vcfIn)
	if args.demo:
		data.append('-demo=' + args.demo)
	if args.weights:
		data.append('-weights=' + args.weights)
	if args.out:
		data.append('-out=' + args.out)
	if args.ancestral:
		data.append('-ancestral=' + args.ancestral)
	if args.refgenome:
		data.append('-refgenome=' + args.refgenome)
	if args.haploid:
		data.append('-haploid=' + args.haploid)
	if args.conditional:
		data.append('-conditional=' + args.conditional)

subprocess.run(data)
