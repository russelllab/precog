#!/usr/bin/env python2.7
# coding: utf-8
import os, sys, time, datetime, argparse, numpy
path = os.getcwd()
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from sklearn.externals import joblib as jl
print path
sys.path.insert(0, path+'/data/')
import load_functions
timestamp1 = time.time()

##########################
## To Parse the arguments
##########################

parser = argparse.ArgumentParser(description='PRECOG (PREdicting COupling probabilities Of G-protein coupled receptors)', epilog='End of help. Contact: gurdeep.singh@bioquant.uni-heidelberg.de')
parser.add_argument('fasta_file', help='path of input file(FASTA format); see data/sample.fasta')
parser.add_argument('--hmm', help='path to the input file\'s hmmsearch o/p(against 7tm_1); if absent, this script will generate one for itself using default settings of HMM')
parser.add_argument('--o', help='path to the output file; if absent, the output will be printed on the screen')
args = parser.parse_args()

fasta_file = args.fasta_file
if fasta_file[0] != '/':
	fasta_file = path + '/' + args.fasta_file
hmm_file = args.hmm
if hmm_file != None:
	if hmm_file[0] != '/':
		hmm_file = path + '/' + hmm_file
out_file = args.o
if out_file != None:
	if out_file[0] != '/':
		out_file = path + '/' + out_file
os.system('clear')
print 'PRECOG v1.0\n##########'
print '#####\n'
##########################

################################
## Base class for all sequences
################################
class base:

	def __init__(self, name, mut):
		self.name = name
		self.seq = ''
		self.til = ''
		self.til_start = 0
		self.til_end = 0
		self.ctl = ''
		self.ctl_start = 0
		self.mutation = {}
		self.mutation['WT'] = {}
		self.mutation['WT']['position'] = {}
		self.add_mut(mut)

	def add_mut(self, mut):
		if mut != '':
			if self.mutation.has_key(mut) == False:
				self.mutation[mut] = {}
				self.mutation[mut]['position'] = {}

	def add_fasta(self, seq):
		self.seq += seq

	def show(self):
		print self.name
		print self.seq
		print self.mutation
		print self.til
		print len(self.til.replace('-', ''))
		print self.til_start
		print self.til_end
		print self.ctl
		print self.ctl_start

#################################
## Read the input FASTA sequences
#################################
def read_fasta(fasta_file):
	obj = {}
	for line in open(fasta_file, 'r'):
		if line[0] == '#':
			continue
		if line[0] != '\n' or line[0] != '':
			if line[0] == '>':
				mutation = ''
				name = (line.split('>')[1].replace('\n', '').split())[0]

				if '/' in name:
					mutation = name.split('/')[1]
					name = name.split('/')[0]
					if mutation[-1] == 'a' or mutation[-1] == 'p' or mutation[-1] == 'X' or mutation[0].isalpha() == False or mutation[-1].isalpha() == False or mutation[1:-1].isdigit() == False:
						go = 0
						name = ''
						continue
				else:
					mutation = ''

				if obj.has_key(name) == False:
					obj[name] = base(name, mutation)
					go = 1
				else:
					obj[name].add_mut(mutation)
					go = 0

			else:
				if go == 1:
					obj[name].add_fasta(line.replace(' ', '').replace('\n', ''))
	return obj

obj = read_fasta(fasta_file)
print 'Read FASTA'
#################################

###############################################
## Read the hmmsearch o/p of input against 7tm1
## or generate one if not provided by the user
###############################################

def extract_other_features(file):
	v = 0; w = 0; proteins = []
	for line in open(file, 'r'):
		if line[0:2] == '>>':
			name = (line.split('>>')[1].replace('\n', '').split())[0]
			if '/' in name:
				name = name.split('/')[0]

			if obj.has_key(name) == True:
				if name not in proteins:
					proteins.append(name)
					proteins = list(set(proteins))
					v = 1
				else:
					v = 0
			else:
				v = 0
			continue

		if v == 1:
			if '7tm_1' in line[:25]:
				x = line.split()
			elif name in line[:25]:
				y = line.split()
				obj[name].ctl = obj[name].seq[int(y[-1]):].upper()
				obj[name].ctl_start = int(y[-1]) + 1
				w = 1

		if w == 1:
			start = int(x[1])
			count = 0

			for k, (i, j) in enumerate(zip(x[2], y[2])):
				xseq = x[2][:k].replace('-', '').replace('.', '')
				if len(xseq) + start >= 172 and len(xseq) + start <= 204:
					obj[name].til += j.upper()
				if len(xseq) + start >= 172 and obj[name].til_start == 0:
					obj[name].til_start = int(y[1]) + len((y[2])[:k].replace('-', '').replace('.', ''))
				if len(xseq) + start >= 172 and len(xseq) + start <= 204 and obj[name].til_end >= 0:
					if len((y[2])[:k].replace('-', '').replace('.', '')) > 0:
						obj[name].til_end = int(y[1]) + len((y[2])[:k+1].replace('-', '').replace('.', '')) - 1
						'''
						if name == 'OR7G2':
							print '----'
							print y[1], len((y[2])[:k].replace('-', '').replace('.', ''))
							print obj[name].til_start, obj[name].til_end, k, count+start
							print (y[2])[:k+1].replace('-', '').replace('.', '')
						'''
				if i!='-' and i!='.':
					count += 1

			w = 0

if os.path.exists(path+'/temp') == False:
	os.system('mkdir '+path+'/temp')

## If the file is not provided by user, generate one itself
if hmm_file == None:
	hmm_file = load_functions.hmm_search(obj, path)

extract_other_features(hmm_file)
print 'Done reading HMMSEARCH against 7tm1'
###########################################################


###################################
## Generate HMM search o/p of input
## against our Gproteins models
###################################

def spinning_cursor():
    while True:
        for cursor in '|/-\\':
            yield cursor

## Create a new fasta file containing only those genes
## that are present in both the fasta and hmmsearch
def new_fasta():
	l = ''
	for name in obj:
		l += '>' + str(name) + '\n' + str(obj[name].seq) + '\n'
	if os.path.exists(path+'/temp') == False:
		os.system('mkdir temp')
	open(path+'/temp/new_fasta_file.txt', 'w').write(l)

new_fasta()
spinner = spinning_cursor()

print 'Generating HMMSEARCH o/p of the input against Gproteins models.....',

os.chdir(path+'/data/hmm_models/')

for files in os.listdir('.'):
	if files.endswith('.hmm'):
		os.system('hmmsearch '+files+' '+path+'/temp/new_fasta_file.txt'+' >' +path+ '/temp/'+files.split('.hmm')[0]+'.out')
		sys.stdout.write(spinner.next())
		sys.stdout.flush()
		time.sleep(0.01)
		sys.stdout.write('\b')

print 'Completed.'
###########################################################

#####################################
## Functions for feature construction
#####################################

def construct_features(name, pos, neg, hmm_pos, hmm_neg, features, row, *mut):
	df = pd.read_table(path+'/data/selected_features.txt', sep = '\t', index_col = 0)
	map_position = {}
	for x, y in df[['MSA_Pos', 'Domain_Pos']].as_matrix().tolist():
		map_position[str(x)] = str(y.replace('(', '|').replace(')', ''))
	#print map_position

	mutation = ''
	for arg in mut:
		mutation = arg
	#print mutation
	mutation_position = -1
	if mutation != 'WT':
		mutation_position = int(mutation[1:-1])
	#print mutation_position
	#print features
	dom_position = {}

	for f in features:
		remarks = []

		if 'pos' in f:
			if pos[name].has_key(str(f[:-3])) == True:
				if pos[name][str(f[:-3])]['position'] == mutation_position:
					AA = mutation[-1]
					dom_position[map_position[str(f[:-3])]] = f
				else:
					AA = pos[name][str(f[:-3])]['aa']

				if AA != '-':
					row.append(1)
				else:
					row.append(0)
			else:
				row.append(0)

			#if f in ['966pos', '967pos'] and mutation in ['T655M', 'WT']:
			#	print f, row[-1], mutation, AA

		elif 'neg' in f:
			if neg[name].has_key(str(f[:-3])) == True:
				if neg[name][str(f[:-3])]['position'] == mutation_position:
					AA = mutation[-1]
					dom_position[map_position[str(f[:-3])]] = f
				else:
					AA = neg[name][str(f[:-3])]['aa']

				if AA != '-':
					row.append(1)
				else:
					row.append(0)
			else:
				row.append(0)

		if 'bip' in f:
			if pos[name].has_key(str(f[:-3])) == True:
				if pos[name][str(f[:-3])]['position'] == mutation_position:
					AA = mutation[-1]
					dom_position[map_position[str(f[:-3])]] = f
				else:
					AA = pos[name][str(f[:-3])]['aa']

				if AA != '-':
					row.append(float(hmm_pos[str(f[:-3])][AA.upper()]))
				else:
					row.append(float(load_functions.find_max_bits(hmm_pos, f)))
			else:
				row.append(float(load_functions.find_max_bits(hmm_pos, f)))

		elif 'bin' in f:
			if neg[name].has_key(str(f[:-3])) == True:
				if neg[name][str(f[:-3])]['position'] == mutation_position:
					AA = mutation[-1]
					dom_position[map_position[str(f[:-3])]] = f
				else:
					AA = neg[name][str(f[:-3])]['aa']

				if AA != '-':
					row.append(float(hmm_neg[str(f[:-3])][AA.upper()]))
				else:
					row.append(float(load_functions.find_max_bits(hmm_neg, f)))
			else:
				row.append(float(load_functions.find_max_bits(hmm_neg, f)))

		elif 'TILL' in f:
			if mutation_position >= obj[name].til_start and mutation_position <= obj[name].til_end:
				sequence = ''
				if obj[name].til != '':
					sequence = list(obj[name].til.replace('-', ''))
					#print sequence, mutation, obj[name].til_start, obj[name].til_end, name, 'til'
					sequence[mutation_position - obj[name].til_start] = mutation[-1]
					sequence = "".join(sequence)
				dom_position['TILL'] = f
			else:
				sequence = obj[name].til
			row.append(len(sequence))

		elif 'CTLL' in f:
			if mutation_position >= obj[name].ctl_start:
				sequence = ''
				if obj[name].ctl != '':
					sequence = list(obj[name].ctl.replace('-', ''))
					#print sequence, mutation, obj[name].ctl_start, name, 'ctl'
					sequence[mutation_position - obj[name].ctl_start] = mutation[-1]
					sequence = "".join(sequence)
				dom_position['CTLL'] = f
			else:
				sequence = obj[name].ctl
			row.append(len(sequence))

		elif '_CTL' in f:
			if mutation_position >= obj[name].ctl_start:
				sequence = ''
				if obj[name].ctl != '':
					sequence = list(obj[name].ctl.replace('-', ''))
					sequence[mutation_position - obj[name].ctl_start] = mutation[-1]
					sequence = "".join(sequence)
					dom_position[f] = f
			else:
				sequence = obj[name].ctl
			row.append(sequence.count(f.split('_')[0]))

		elif '_TIL' in f:
			if mutation_position >= obj[name].til_start and mutation_position <= obj[name].til_end:
				sequence = ''
				if obj[name].til != '':
					sequence = list(obj[name].til.replace('-', ''))
					#print sequence, mutation, obj[name].til_start, obj[name].til_end, name, 'til'
					sequence[mutation_position - obj[name].til_start] = mutation[-1]
					sequence = "".join(sequence)
				dom_position[f] = f
			else:
				sequence = obj[name].til
			row.append(sequence.count(f.split('_')[0]))

		#if mutation != 'WT':
			#print f, remarks
			#sys.exit()
	#print mutation, row
	return row, dom_position

def read_aln(pos, neg, hmm_pos, hmm_neg, l, features, gprotein):
	data = []
	for name in obj:
		if pos.has_key(name) == True and neg.has_key(name) == True:
			for mut in obj[name].mutation:
				if mut == 'WT' or int(mut[1:-1]) <= len(obj[name].seq):
					row = []
					row.append(str(name))
					row.append(str(mut))
					row, dom_position = construct_features(name, pos, neg, hmm_pos, hmm_neg, features, row, mut)

					for pi in dom_position:
						if obj[name].mutation[mut]['position'].has_key(pi) == False:
							obj[name].mutation[mut]['position'][pi] = {}
							obj[name].mutation[mut]['position'][pi][dom_position[pi]] = []
							obj[name].mutation[mut]['position'][pi][dom_position[pi]].append(str(gprotein))
						else:
							if obj[name].mutation[mut]['position'][pi].has_key(str(dom_position[pi])) == False:
								obj[name].mutation[mut]['position'][pi][dom_position[pi]] = []
								obj[name].mutation[mut]['position'][pi][dom_position[pi]].append(str(gprotein))
							else:
								obj[name].mutation[mut]['position'][pi][dom_position[pi]].append(str(gprotein))
								obj[name].mutation[mut]['position'][pi][dom_position[pi]] = list(set(obj[name].mutation[mut]['position'][pi][dom_position[pi]]))
					data.append(row)
	#print data[3]
	#print data[1]
	return numpy.array(data)

def read_gprotein_hmm_out(file, hmm):
	v = 0; w = 0; gpcr = {}
	for line in open(file, 'r'):
		if line[0:2] == '>>':
			name = (line.split('>> ')[1].replace('\n', '').split())[0]
			if '/' in name:
				name = name.split('/')[0]
			if obj.has_key(name) == True:
				gpcr[name] = {}
				v = 1
			else:
				name = ''
				v = 0
			continue
		if v == 1:
			file_name = file.split('.out')[0].split('/')[-1]
			file_name = file_name.split('_')[0] + '_subali_' + file_name.split('_')[1]
			if file_name in line[:25]:
				x = line.split()
			elif name in line[:25]:
				y = line.split()
				if y[1] != '-' and y[-1] != '-':
					w = 1
		if w == 1:
			start = int(x[1])
			count = 0
			#print x
			#print y
			for k, (i, j) in enumerate(zip(x[2], y[2])):
				if i!='-' and i!='.':
					if hmm.has_key(str(count+start)) == True:
						#if 'GNA12' in file and 'pos' in file:
							#print count+start,
						if int(hmm[str(count+start)]) >= 393 and int(hmm[str(count+start)]) <= 1002:
							#if 'GNA12' in file and 'pos' in file:
							#	print hmm[str(count+start)]
							gpcr[name][str(hmm[str(count+start)])] = {}
							gpcr[name][str(hmm[str(count+start)])]['aa'] = j
							if j == '-':
								gpcr[name][str(hmm[str(count+start)])]['position'] = '-'
							else:
								gpcr[name][str(hmm[str(count+start)])]['position'] = int(y[1]) + len((y[2])[:k].replace('-', '').replace('.', ''))
					count += 1
				#print k+int(y[1]),
			#print
			w = 0

	return gpcr

def extract_features(file):
	gprotein = files.split('_')[0]
	features = []
	for line in open(file, 'r'):
		if gprotein not in line:
			features.append(line.split('\t')[0])
	return features

def extract_model(gprotein):
	#for files in os.listdir('/net/netfile2/ag-russell/bq_gsingh/gpcr/update_2/output_VI/'):
	for files in os.listdir(path+'/data/output/'):
		if gprotein in files and 'model' in files:
			model = jl.load(files)
			break
	return model

def k_fold(file):
	df = pd.read_table(file, lineterminator = '\n', sep = '\t')
	col = list(df.columns.values)
	df[col[1:-1]] = df[col[1:-1]].astype(float)
	min_max_scaler_all = MinMaxScaler()
	min_max_scaler_all.fit_transform(df[col[1:-1]])
	return min_max_scaler_all


#print 'Making predictions for:'
gpcr_list = obj.keys()
prediction = {}
gprotein_list = ['GNAI3', 'GNAI1', 'GNAZ', 'GNAO1', 'GNA12', 'GNA13', 'GNAQ', 'GNA14', 'GNA15', 'GNAS', 'GNAL']
os.chdir(path+'/data/output/')
for gprotein in gprotein_list:
	#print gprotein + '... ',
	for files in os.listdir('.'):
		if gprotein in files and 'fweight' in files:
			features = extract_features(files)
			#print features

			hmm_pos, hmm_pos_positions = load_functions.read_hmm(path+'/data/hmm_models/'+gprotein+'_pos.hmm')
			hmm_neg, hmm_neg_positions = load_functions.read_hmm(path+'/data/hmm_models/'+gprotein+'_neg.hmm')

			pos = read_gprotein_hmm_out(path+'/temp/'+gprotein+'_pos.out', hmm_pos_positions)
			neg = read_gprotein_hmm_out(path+'/temp/'+gprotein+'_neg.out', hmm_neg_positions)

			l= 'GPCR\t'
			for f in features:
				l+=f+'\t'
			l+= '\n'
			data = read_aln(pos, neg, hmm_pos, hmm_neg, l, features, gprotein)
			#if gprotein == 'GNAI3':
			#	print data

			feature_matrix = data[:, 2:]
			model = extract_model(gprotein)

			min_max = k_fold(path+'/data/feature_files/'+str(gprotein)+'_train.txt')

			feature_matrix = min_max.transform(numpy.array(feature_matrix))
			Y = model.predict(feature_matrix)
			Y_prob = model.predict_proba(feature_matrix)

			for (name, mut), y in zip(data[:, :2], Y_prob):
				obj[name].mutation[mut][gprotein] = round(y[1], 3)

			#print 'Completed.'
			break

print '\nWriting the output...\n'

timestamp2 = time.time()
now = datetime.datetime.now()
l = '# PRECOG (v1.0)\n'
l += '# http://precog.russelllab.org\n'
l += '# Contact: gurdeep[dot]singh[at]bioquant[dot]uni-heidelberg[dot]de\n'
l += '# Run time: ' + str(round(timestamp2 - timestamp1, 2)) + ' sec\n'
l += '# Date: ' + str(now.day) + '-' + str(now.month) + '-' + str(now.year) + '\n'
l += '# Time: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second) + '\n'
l += '# Please cite PRECOG:' + '\n'
l += "# Gurdeep Singh, Asuka Inoue, J Silvio Gutkind, Robert B Russell, Francesco Raimondi,\n# PRECOG: PREdicting COupling probabilities of G-protein coupled receptors, Nucleic Acids Research,\n# Volume 47, Issue W1, 02 July 2019, Pages W395–W401, https://doi.org/10.1093/nar/gkz392\n#\n"

new_name = '#GPCR/MUT'
while len(new_name) <= 20:
	new_name += ' '
l += new_name
gprotein_list = ['GNAI3', 'GNAI1', 'GNAZ', 'GNAO1', 'GNA12', 'GNA13', 'GNAQ', 'GNA14', 'GNA15', 'GNAS', 'GNAL']
for gprotein in gprotein_list:
	l+='\t'+gprotein
l+='\t7TM1_POS/BW/ALN_POS\tMutation_Info\n'

d = load_functions.load_iuphar(path+'/data/IUPHAR_couplings.tsv')

mut_info = load_functions.load_mut_info(path+'/data/mechismo_input_uniprot_muts_mods_v5.txt')

for name in obj:
	l += load_functions.check_iuphar(name, gprotein_list, d)
	l += load_functions.check_aska(name, gprotein_list, path)
	new_name = str(name) + '/' + str('WT')
	while len(new_name) <= 20:
		new_name += ' '
	l += new_name
	for gprotein in gprotein_list:
		if obj[name].mutation['WT'].has_key(gprotein):
			l += '\t' + str(obj[name].mutation['WT'][gprotein])
		else:
			l += '\t' + '-'
	l+='\tWT\n'
	for mut in obj[name].mutation:
		if mut != 'WT':
			new_name = str(name) + '/' + str(mut)
			mutation = new_name
			while len(new_name) <= 20:
				new_name += ' '
			l += new_name
			for gprotein in gprotein_list:
				if obj[name].mutation[mut].has_key(gprotein):
					l += '\t' + str(obj[name].mutation[mut][gprotein])
				else:
					l += '\t' + '-'
			l += '\t'
			new_name = ''
			for dp in obj[name].mutation[mut]['position']:
				new_name += str(dp.replace('|', '/')) + '/'
				for key in obj[name].mutation[mut]['position'][dp]:
					new_name += str(key.replace('|', '/'))
					for item in obj[name].mutation[mut]['position'][dp][key]:
						new_name += '/' + str(item)
					new_name += ','
			l += new_name[:-1] + '\t'
			if mut_info.has_key(mutation) == True:
				l += mut_info[mutation]
			else:
				l += '\t'
			l += '\n'

if out_file == None:
	print l
else:
	open(out_file, 'w').write(l)
	print 'Output saved at: '+str(out_file)
print 'Completed.\n'
print 'End of PRECOG\n#####'
print '##########\n'
sys.exit()
