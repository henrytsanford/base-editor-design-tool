"""The design engine: sgRNA designs for a transcript, as rows.

Everything here is importable and free of process-wide state. Reference data
arrives as a `source` argument (an Ensembl REST client or a local bundle) and
ClinVar as a `clinvar` argument, so a caller can hold several open at once --
which is what the web service does, one per worker process.

The three outputs are the three files the CLI writes, as lists of rows under
DESIGN_COLUMNS, ERROR_COLUMNS and ANNOTATION_COLUMNS.
"""
import re
from dataclasses import dataclass

from .transcript_source import empty_variants

# The header of sgrna_designs_<name>.txt
DESIGN_COLUMNS = [
	'sgRNA sequence', 'sgRNA context sequence', 'Gene Symbol', 'Ensembl Gene ID',
	'Ensembl transcript ID', 'Gene strand', 'Genome assembly',
	'Transcript reference allele', 'Transcript alternate allele',
	'Genome reference allele', 'Genome alternate allele', 'Chromosome',
	'sgrna genomic position', 'sgRNA Strand', 'PAM', 'Edit', '# edits',
	'#silent edits', 'Nucleotide edits', 'Amino acid edits', 'Mutation category',
	'Clinical significance', 'BsmBI flag', '4T flag']

# The header of error_report.txt
ERROR_COLUMNS = ['Gene Symbol', 'Ensembl transcript ID', 'sgRNA', 'sgRNA Strand',
				 'Error']

# The header of clinvar_annotations_<name>.txt. Rows are ragged by design: a row
# carries all 28 fields when a ClinVar SNP matched the edit and stops after
# 'Mutation category' when none did.
ANNOTATION_COLUMNS = [
	'sgRNA sequence', 'sgRNA strand', 'sgRNA context sequence', 'Chromosome',
	'Gene', 'Gene strand', 'Edit', 'Transcript reference allele',
	'Transcript alternate allele', 'Genome reference allele',
	'Genome alternate allele', 'Edit nucleotide', 'Edit nucleotide position(s)',
	'sgRNA amino acid change', 'Original codon', 'Edited codon',
	'Mutation category', 'SNP amino acid change', 'SNP name',
	'SNP clinical significance', 'SNP nucleotide position',
	'SNP reference allele', 'SNP alternate allele', 'SNP review status',
	'Same nucleotide position', 'Same nucleotide change',
	'Same amino acid position', 'Same amino acid change']

# The files a run writes; '%s' takes the run's output name. The engine owns the
# names because it owns the formats. The CLI and the golden tests read them here.
DESIGNS_FILE = 'sgrna_designs_%s.txt'
ERRORS_FILE = 'error_report.txt'
ANNOTATIONS_FILE = 'clinvar_annotations_%s.txt'

# 'all' means both deaminase directions, run in this order.
ALL_EDITS = ('C-T', 'A-G')


class UnknownBaseEditor(ValueError):
	"""--be-type named an editor that is not in the preset table.

	A ValueError, so a caller can treat it as the bad argument it is, but named
	so the CLI can report it without swallowing unrelated ValueErrors too.
	"""


@dataclass(frozen=True)
class DesignParams:
	"""Everything about a run that is not the transcript or the reference data."""
	pam: str = 'NGG'
	window: str = '4-8'
	sg_len: int = 20
	edit: str = 'all'
	intron_buffer: int = 30
	filter_gc: bool = False

	@classmethod
	def from_preset(cls, be_type, **overrides):
		"""Resolves a base-editor name. Raises UnknownBaseEditor if it is not a preset."""
		pam, window, sg_len, edit = get_pam_window_len(be_type)
		return cls(pam=pam, window=window, sg_len=sg_len, edit=edit, **overrides)

	@property
	def edits(self):
		"""The deaminase passes this run makes."""
		return ALL_EDITS if self.edit == 'all' else (self.edit,)


def revcom(s):
	basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N','K':'M','M':'K','R':'Y','Y':'R','S':'S','W':'W','B':'V','V':'B','H':'D','D':'H'}
	letters = list(s[::-1])
	letters = [basecomp[base] for base in letters]
	return ''.join(letters)


'''
Ensembl's REST API rejects versioned transcript IDs (ENST00000294952.13 returns
HTTP 400), so drop the suffix. Also trims stray whitespace, which spreadsheet
exports tend to leave behind.
'''
def strip_tr_version(tr):
	return re.sub(r'^(ENS[A-Z]*[GTP]\d+)\.\d+$', r'\1', str(tr).strip())


def get_aa_map():
	aa_map = {'Phe': 'F', 'Leu': 'L', 'Ile': 'I', 'Met': 'M', 'Val': 'V', 'Ser': 'S', 'Pro': 'P', 'Thr': 'T',
			  'Ala': 'A', 'Tyr': 'Y', 'Ter': 'Ter', 'His': 'H', 'Gln': 'Q', 'Asn': 'N', 'Lys': 'K', 'Asp': 'D',
			  'Glu': 'E', 'Cys': 'C', 'Trp': 'W', 'Arg': 'R', 'Gly': 'G'}
	return aa_map


def get_codon_map():
	codon_map = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I',
				 'ATA':'I', 'ATG':'M', 'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
				 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCT':'A', 'GCC':'A',
				 'GCA':'A', 'GCG':'A', 'TAT':'Y', 'TAC':'Y', 'TAA':'Ter', 'TAG':'Ter', 'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
				 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'TGT':'C', 'TGC':'C',
				 'TGA':'Ter', 'TGG':'W', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
				 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}
	return codon_map


def get_pam_window_len(be):
	be_types = {'BE1':'NGG_4-8_20_C-T', 'BE2':'NGG_4-8_20_C-T', 'BE3':'NGG_4-8_20_C-T', 'HF-BE3':'NGG_4-8_20_C-T', 'BE4':'NGG_4-8_20_C-T', 'BE4max':'NGG_4-8_20_C-T',
			   'BE4-Gam': 'NGG_4-8_20_C-T', 'YE1-BE3':'NGG_4-7_20_C-T', 'EE-BE3':'NGG_5-6_20_C-T', 'YE2-BE3':'NGG_5-6_20_C-T', 'YEE-BE3':'NGG_5-6_20_C-T',
			   'VQR-BE3': 'NGAN_4-11_20_C-T', 'VRER-BE3': 'NGCG_3-10_20_C-T', 'SaBE3': 'NNGRRT_3-12_21_C-T', 'SaBE4': 'NNGRRT_3-12_21_C-T',
			   'SaBE4-Gam': 'NNGRRT_3-12_21_C-T', 'Sa(KKH)-BE3': 'NNNRRT_3-12_21_C-T', 'Target-AID': 'NGG_2-4_20_C-T', 'Target-AID-NG': 'NG_2-4_20_C-T',
			   'xBE3': 'NG_4-8_20_C-T', 'eA3A-BE3': 'NG_4-8_20_C-T', 'A3A-BE3': 'NG_4-8_20_C-T', 'BE-PLUS':'NGG_4-14_20_C-T', 'ABE7.9':'NGG_5-8_20_A-G',
			   'ABE7.10': 'NGG_4-7_20_A-G', 'xABE':'NG_4-7_20_A-G', 'ABESa':'NNGRRT_6-12_21_A-G', 'VQR-ABE':'NGA_4-6_20_A-G', 'VRER-ABE':'NGCG_4-6_20_A-G',
			   'Sa(KKH)-ABE':'NNNRRT_6-12_21_A-G'}
	if be in be_types.keys():
		pam, window, sg_len, edit = be_types[be].split('_')
	else:
		raise UnknownBaseEditor(
			'Unknown base editor %r. Please enter ONE of the following: %s'
			% (be, ','.join(be_types)))
	return pam, window, int(sg_len), edit


def get_pam_pattern(pam):
	code = {'N':'ACTG', 'R':'AG', 'Y':'CT', 'S':'GC', 'W':'AT', 'K':'GT', 'M':'AC', 'B':'CGT', 'D':'AGT', 'H':'ACT', 'V':'ACG'}
	pattern = ''
	for p in pam:
		if p in code.keys():
			pattern = pattern +'['+ code[p] + ']'
		else:
			pattern+=p
	return pattern


def check_ressite_4t(sg):
	res_flag, t4_flag = '', ''
	if ('CGTCTC' in sg or 'GAGACG' in sg or sg.startswith('TCTC') or sg.startswith('AGACG') or sg.endswith('GAGAC')):
		res_flag = 'yes'
	if 'TTTT' in sg:
		t4_flag = 'yes'
	return res_flag, t4_flag


'''
Returns information about gene, assembly, chromosome of specified Ensembl transcript; Also returns absolute values for gene with respect to
genomic locations; Flags genes for absence of UTRs;
'''
def get_tr_info(source, tr, input_type):
	tr_info = source.lookup(tr)
	gene_name = tr_info['display_name'].rsplit('-', 1)[0]
	assembly = tr_info['assembly_name']
	gene_strand = tr_info['strand']
	chromosome = tr_info['seq_region_name']
	gene_id = tr_info['Parent']
	if gene_strand == 1:
		gene_start = tr_info['start']
		gene_end = tr_info['end']
		length = gene_end - gene_start
	else:
		gene_start = tr_info['end']
		gene_end = tr_info['start']
		length = gene_start - gene_end
	abs_pos_map, fs = get_absolute_pos(gene_start,gene_end,gene_strand, input_type)
	exons, cds_map = get_exons(source,tr,length,gene_start,gene_strand)
	if exons != '':
		utr, cds_start_exon, utr5_flag, utr3_flag = get_utrs(tr_info,exons,gene_strand)
	else:
		utr = ''
		cds_start_exon = ''
		utr5_flag = ''
		utr3_flag = ''
	return gene_name, assembly, gene_strand, chromosome, gene_id, exons, cds_map, abs_pos_map, fs, utr, cds_start_exon, utr5_flag, utr3_flag


'''
Returns exon boundaries for transcript;Also returns absolute values for gene with respect to
genomic locations
'''
def get_exons(source, tr, length, gene_start, gene_strand):
	exons_all = source.cds_mappings(tr, length)
	if exons_all is None:
		# No CDS (e.g. a non-coding transcript)
		return '',''
	if gene_strand == 1:
		gene_start_pos = exons_all[0]['start'] - gene_start
	else:
		gene_start_pos = gene_start - exons_all[0]['end']
	exons = []
	for i,e in enumerate(exons_all):
		exons.append(str(e['start'])+':'+str(e['end']))
	#abs_pos_map = get_absolute_pos(exons,gene_start_pos,gene_strand)
	cds_map = get_cds_map(exons, gene_strand)
	return exons, cds_map


'''
Generates a hash of utr locations
'''
def get_utr_map(utr, start, end, strand,j):
	k=start
	if strand == 1:
		while k<=end:
			utr[k] = j
			k+=1
			j+=1
	elif strand == -1:
		while k>=end:
			utr[k] = j
			k-=1
			j+=1
	return utr,j


'''
Returns hash of coordinates of 5' and 3' UTR
'''
def get_utrs(tr_info, exons, gene_strand):
	utr_exons = tr_info['Exon']
	ce_1 = exons[0]
	ce_last = exons[-1]
	if gene_strand == 1:
		ce_start,ce_end= ce_1.split(':')
		ce_last_start, ce_last_end = ce_last.split(':')
	elif gene_strand == -1:
		ce_end,ce_start= ce_1.split(':')
		ce_last_end,ce_last_start = ce_last.split(':')
	utr = {}
	j=0
	utr5_flag = 0
	utr3_flag = 0
	for i,e in enumerate(utr_exons):
		if gene_strand == 1:
			if (e['start']<int(ce_start)):
				if (e['end']<int(ce_start)):
					utr,j = get_utr_map(utr,e['start'],e['end'],gene_strand,j)
				elif (e['end']>int(ce_start)):
					cds_start_exon = i+1
					utr,j = get_utr_map(utr,e['start'],int(ce_start)-1,gene_strand,j)
			if (e['end'] > int(ce_last_end)):
				if (e['start']>int(ce_last_end)):
					utr,j = get_utr_map(utr,e['start'],e['end'],gene_strand,j)
				elif (e['start']<int(ce_last_end)):
					utr,j = get_utr_map(utr,int(ce_last_end)+1,e['end'],gene_strand,j)
			if (e['start'] == int(ce_start)):
				if i == 0:
					utr5_flag = 1
					cds_start_exon = i+1
				else:
					cds_start_exon = i+1
			if (e['end'] == int(ce_last_end)):
				if i == len(utr_exons)-1:
					utr3_flag = 1
		elif gene_strand == -1:
			if (e['end']>int(ce_start)):
				if (e['start']>int(ce_start)):
					utr,j = get_utr_map(utr,e['end'],e['start'],gene_strand,j)
				elif (e['start']<int(ce_start)):
					cds_start_exon = i + 1
					utr,j = get_utr_map(utr,e['end'],int(ce_start)+1,gene_strand,j)
			if (e['start']<int(ce_last_end)):
				if (e['end']>int(ce_last_end)):
					utr,j = get_utr_map(utr,int(ce_last_end)-1,e['start'],gene_strand,j)
				elif (e['end']<int(ce_last_end)):
					utr,j = get_utr_map(utr,e['end'],e['start'],gene_strand,j)
			if (e['end'] == int(ce_start)):
				if i == 0:
					utr5_flag = 1
					cds_start_exon = i+1
				else:
					cds_start_exon = i+1
			if (e['start'] == int(ce_last_end)):
				if i == len(utr_exons) - 1:
					utr3_flag = 1
	return utr, cds_start_exon, utr5_flag, utr3_flag


'''
Generates absolute positions for genomic locations of transcript starting from 40 nucleotides flanking sequence preceding
the 5'UTR up to 40 nucleotides flanking sequence succeeding the end of 3'UTR 
'''
def get_absolute_pos(gene_start, gene_end, gene_strand, input_type):
	abs_pos_map = {}
	fs = {}                 #Genomic location of flanking sequences
	#i=gene_start
	j=0
	if gene_strand == 1:
		if input_type == 'tid':
			i = gene_start-40   #including 40 nucs of flanking sequence annotation to abs_pos_map
			end = gene_end+40
		else:
			i = gene_start
			end = gene_end
		while i<=end:
			abs_pos_map[i] = j
			if i < gene_start or i > gene_end:
				fs[i] = j
			i+=1
			j+=1

	else:
		if input_type == 'tid':
			i = gene_start+40   #including 40 nucs of flanking sequence annotation to abs_pos_map
			end = gene_end-40
		else:
			i = gene_start
			end = gene_end
		while i>=end:
			abs_pos_map[i] = j
			if i > gene_start or i < gene_end:
				fs[i] = j
			i-=1
			j+=1

	return abs_pos_map, fs


'''
Generates absolute cds positions for genomic locations of transcript
'''
def get_cds_map(exons, gene_strand):
	cds_map={}
	j=0
	for e in exons:
		if gene_strand == 1:
			e = e.split(':')
			i = int(e[0])
			while(i <= int(e[1])):
				cds_map[i] = j
				i=i+1
				j=j+1
		else:
			e = e.split(':')
			i = int(e[1])
			while(i >= int(e[0])):
				cds_map[i] = j
				i-=1
				j+=1
	#j=j-1
	return cds_map


'''
Returns sequence of transcript along with 5' and 3' UTR
'''
def get_tr_sequence(source, tr):
	return source.genomic_sequence(tr)


'''
Returns protein sequence of transcript
'''
def get_pro_sequence(source, tr):
	return source.protein_sequence(tr)


'''
Returns CDS sequence of transcript
'''
def get_cds_sequence(source, tr):
	return source.cds_sequence(tr)


'''
Translates sgRNA sequence and annotates frame
'''
def get_sgrna_translated_seq(sgrna, cds_map, abs_pos, pos_for_index, fs, sgrna_start_pos, gene_strand, sgrna_strand, utr, e, label):
	map_key = pos_for_index[sgrna_start_pos]
	sgrna_trans = {}
	if sgrna_strand == 'sense':
		for i,n in enumerate(sgrna):
			if gene_strand == 1:
				if map_key not in cds_map.keys():
					if map_key in utr.keys():
						sgrna_trans[n+str(i+1)] = '0_U'
					elif map_key in fs.keys():
						sgrna_trans[n + str(i + 1)] = '0_FS'
					else:
						if map_key < int(e[0]):
							sgrna_trans[n+str(i+1)] = '0_'+label+':-'+str(int(e[0]) - map_key)
						elif map_key > int(e[1]):
							sgrna_trans[n+str(i+1)] = '0_'+label+':+'+str(map_key - int(e[1]))
					map_key +=1
				else:
					sgrna_trans[n+str(i+1)] = str((cds_map[map_key]//3)+1)+'_'+str(cds_map[map_key]%3)
					map_key +=1
			else:
				if map_key not in cds_map.keys():
					if map_key in utr.keys():
						sgrna_trans[n+str(i+1)] = '0_U'
					elif map_key in fs.keys():
						sgrna_trans[n + str(i + 1)] = '0_FS'
					else:
						if map_key < int(e[0]):
							sgrna_trans[n+str(i+1)] = '0_'+label+':+'+str(int(e[0]) - map_key)
						elif map_key > int(e[1]):
							sgrna_trans[n+str(i+1)] = '0_'+label+':-'+str(map_key - int(e[1]))
					map_key-=1
				else:
					sgrna_trans[n+str(i+1)] = str((cds_map[map_key]//3)+1)+'_'+str(cds_map[map_key]%3)
					map_key -=1
	elif sgrna_strand == 'antisense':
		for i,n in enumerate(sgrna):
			if gene_strand == 1:
				if map_key not in cds_map.keys():
					if map_key in utr.keys():
						sgrna_trans[n+str(len(sgrna)-i)] = '0_U'
					elif map_key in fs.keys():
						sgrna_trans[n + str(len(sgrna) - i)] = '0_FS'
					else:
						if map_key < int(e[0]):
							sgrna_trans[n+str(len(sgrna)-i)] = '0_'+label+':-'+str(int(e[0]) - map_key)
						elif map_key > int(e[1]):
							sgrna_trans[n+str(len(sgrna)-i)] = '0_'+label+':+'+str(map_key - int(e[1]))
					map_key+=1
				else:
					sgrna_trans[n+str(len(sgrna)-i)] = str((cds_map[map_key]//3)+1)+'_'+str(cds_map[map_key]%3)
					map_key+=1
			else:
				if map_key not in cds_map.keys():
					if map_key in utr.keys():
						sgrna_trans[n+str(len(sgrna)-i)] = '0_U'
					elif map_key in fs.keys():
						sgrna_trans[n + str(len(sgrna) - i)] = '0_FS'
					else:
						if map_key < int(e[0]):
							sgrna_trans[n+str(len(sgrna)-i)] = '0_'+label+':+'+str(int(e[0]) - map_key)
						elif map_key > int(e[1]):
							sgrna_trans[n+str(len(sgrna)-i)] = '0_'+label+':-'+str(map_key - int(e[1]))
					map_key-=1
				else:
					sgrna_trans[n+str(len(sgrna)-i)] = str((cds_map[map_key]//3)+1)+'_'+str(cds_map[map_key]%3)
					map_key -=1
	return sgrna_trans


'''
Parses the ClinVar SNP name
'''
def parse_snp_name(snp_name, aa_map):
	snp_aa = ''
	snp_aa_to = ''
	snp_aa_from = ''
	snp_aa_num = ''
	# If the SNP is in an exon, pull amino acid change out of ClinVar SNP name
	if snp_name.find('(p.') != -1:
		snp_aa = snp_name[snp_name.find('(p.')+3:-1]
		three_letter_aa = re.split(r'(\D+)', snp_aa)
		snp_aa_from = three_letter_aa[1]
		snp_aa_to = three_letter_aa[3]
		snp_aa_num = three_letter_aa[2]
		# Silent mutations are annotated as '='
		if snp_aa_to == '=':
			snp_aa_to = snp_aa_from
		snp_aa = snp_aa_from + snp_aa_num + snp_aa_to
	return snp_aa, snp_aa_from, snp_aa_num, snp_aa_to


'''
Returns:
	codon_pos_list, a list of 3 genomic positions for each nucleotide in a codon
	edit_gen_pos_list, a list of the genomic positions for all edited nucleotides in the codon
'''
def get_genomic_pos_list(edit_indices, gene_strand, sgrna_strand, sg_gen_pos):
	edit_gen_pos_list = []
	codon_pos_list = []
	for edit_pos,frame in edit_indices.items():
		# if sgRNA is in + strand
		if ((gene_strand == 1) and (sgrna_strand == 'sense')) or ((gene_strand == -1) and (sgrna_strand == 'antisense')):
			edit_gen_pos = (sg_gen_pos+int(edit_pos)-1)
		# elif sgRNA is in - strand
		elif ((gene_strand == 1) and (sgrna_strand == 'antisense')) or ((gene_strand == -1) and (sgrna_strand == 'sense')):
			edit_gen_pos = (sg_gen_pos-int(edit_pos)+1)
		edit_gen_pos_list.append(str(edit_gen_pos))
		# For edits in exon, add the other two positions in the codon to the codon_pos_list
		# Therefore we can annotate any SNPs that affect that codon
		if (frame != 'U') and ('Exon' not in frame) and (frame!='FS'):
			frame = int(frame)
			if gene_strand == 1:
				codon_pos_list = [edit_gen_pos-frame,edit_gen_pos-frame+1,edit_gen_pos-frame+2]
			elif gene_strand == -1:
				codon_pos_list = [edit_gen_pos+frame,edit_gen_pos+frame-1,edit_gen_pos+frame-2]
		elif (frame == 'U') or ('Exon' in frame) or (frame == 'FS'):
			# For utr and introns, only check for snps at the location of the edit
			codon_pos_list.extend([edit_gen_pos])
	return codon_pos_list, edit_gen_pos_list


'''
Returns dataframe containing information about the pathogenicity and position of SNPs created by edit
'''
def get_snps(edit_map, edit, sg_gen_pos, gene_strand, sgrna_strand, gene_variant_df, aa_map):
	snp_type_list = []
	snp_info = []
	edit_nuc, edit_to = edit.split('-')
	all_snps = gene_variant_df['ClinVar_SNP_Position'].tolist()
	# Iterate through all amino acid changes
	for k,v in edit_map.items():
		temp_snp_type_list = []
		edit_indices = {}
		edits = k.split('_')
		sgrna_edit_nuc = edits[0]
		for i in edits[1:]:
			edit_index, frame = i.split('-',1)
			edit_indices[edit_index] = frame
		edit_cat = v.split('_')[1]
		# UTR or flanking sequence mutations
		if 'UTR' in v or 'Flanking' in v:
			aa_edit,edit_type = v.split('_')
			old_codon, new_codon = '', ''
			aa_from = ''
			aa_num = ''
			aa_to = ''
		# Coding mutations
		elif '_' in v and 'Exon' not in v:
			aa_edit, edit_type, old_codon, new_codon = v.split('_')
			aa_from = aa_edit[0:3]
			aa_num = aa_edit[3:-3]
			aa_to = aa_edit[-3:]
		# Intronic mutations
		elif 'Exon' in v:
			aa_edit = 'intron'
			edit_type, old_codon, new_codon = '', '', ''
			aa_from = ''
			aa_num = ''
			aa_to = ''
		codon_pos_list, edit_gen_pos_list = get_genomic_pos_list(edit_indices, gene_strand, sgrna_strand, sg_gen_pos)
		# Check if there is any overlap between codon_pos_list and all_snps
		if any(i in codon_pos_list for i in all_snps):
			clinvar_snps_df = gene_variant_df[gene_variant_df.ClinVar_SNP_Position.isin(codon_pos_list)].loc[:,['Name','ClinicalSignificance','ClinVar_SNP_Position','ReferenceAllele','AlternateAllele','ReviewStatus']]
			for index,row in clinvar_snps_df.iterrows():
				snp_aa, snp_aa_from, snp_aa_num, snp_aa_to = parse_snp_name(row.Name, aa_map)
				# First check for nucleotide position
				if str(row.ClinVar_SNP_Position) not in edit_gen_pos_list:
					same_nucleotide_pos = False
					same_nucleotide_change = False
				elif str(row.ClinVar_SNP_Position) in edit_gen_pos_list:
					same_nucleotide_pos = True
					# If sgRNA is in the forward strand, C>T or A>G SNPs will be created
					if ((gene_strand == 1) and (sgrna_strand == 'sense') or ((gene_strand == -1) and sgrna_strand == 'antisense')):
						if row.AlternateAllele == edit_to:
							same_nucleotide_change = True
						else:
							same_nucleotide_change = False
					# If sgRNA is in the reverse strand, G>A or T>C SNPs will be created
					elif ((gene_strand == 1) and (sgrna_strand == 'antisense') or ((gene_strand == -1) and sgrna_strand == 'sense')):
						if row.AlternateAllele == revcom(edit_to):
							same_nucleotide_change = True
						else:
							same_nucleotide_change = False
				if snp_aa != '':
					if str(aa_num) != str(snp_aa_num):
						same_aa_pos = False
						if (aa_from == snp_aa_from) and (aa_to == snp_aa_to):
							same_aa_change = True
						else:
							same_aa_change = False
					elif str(aa_num) == str(snp_aa_num):
						same_aa_pos = True
						if (aa_from == snp_aa_from) and (aa_to == snp_aa_to):
							same_aa_change = True
						else:
							same_aa_change = False
				else:
					same_aa_pos = 'N/A'
					same_aa_change = 'N/A'
				# Append the clinical significance of exact match SNPs to snp_type_list
				# Do not require that the aa number matches
				if same_nucleotide_pos and same_nucleotide_change:
					if same_aa_pos and same_aa_change:
						temp_snp_type_list.append(row.ClinicalSignificance)
					elif (same_aa_pos == 'N/A') and (same_aa_change == 'N/A'):
						temp_snp_type_list.append(row.ClinicalSignificance)
				snp_info.append([str(sgrna_edit_nuc+'_'+'_'.join(edit_indices.keys())),
									 ';'.join(edit_gen_pos_list),
									 aa_edit,
									 old_codon,
									 new_codon,
									 edit_cat,
									 snp_aa] +
									 row.tolist() +
									 [same_nucleotide_pos,
									 same_nucleotide_change,
									 same_aa_pos,
									 same_aa_change])
		else:
			snp_info.append([str(sgrna_edit_nuc+'_'+'_'.join(edit_indices.keys())), ';'.join(edit_gen_pos_list),
				aa_edit,
				old_codon,
				new_codon,
				edit_cat])
		if not temp_snp_type_list:
			temp_snp_type_list.append('None')
		snp_type_list.extend(temp_snp_type_list)		
	return snp_type_list, snp_info


'''
Returns clinical significances as ;-separated string
'''
def get_clinical_sig(snp_type_list):
	if not snp_type_list:
		clinical_sig = ''
	else:
		clinical_sig = ';'.join(snp_type_list)
	return clinical_sig


'''Filters out edits that are in a GC motif'''
def filter_gc_motifs(sgrna_context, context_index, sgrna_strand):
	if sgrna_strand == 'sense':
		motif = sgrna_context[context_index-1]
		if motif == 'G':
			return False
	elif sgrna_strand == 'antisense':
		motif = sgrna_context[len(sgrna_context) - context_index-2]
		if motif == 'G':
			return False			
	return True


def filter_gc_motifs_for_aa(sgrna_strand,sgrna_context,codon_start,k):
	if sgrna_strand == 'sense':
		if sgrna_context[codon_start+k-1] == 'G':
			return False
	elif sgrna_strand == 'antisense':
		if sgrna_context[len(sgrna_context) - (codon_start+k)-2] == 'G':
			return False
	return True


'''
Returns edits for sgRNA in specified window
Also returns number of silent edits in window
'''
def get_edits(edit_map, context, window, edit, sgrna_trans, codon_map, j, sgrna_strand, pam, window_start, window_end, aa_map, filter_gc, sgrna_context, gene_strand):
	error = ''
	num_silent = 0
	num_stop = 0
	edit_nuc, edit_to = edit.split('-')
	if sgrna_strand == 'antisense':
		edit_nuc, edit_to = revcom(edit_nuc), revcom(edit_to)
		context_index_track = len(sgrna_trans)-j+len(pam)+3
	for i,n in enumerate(window):
		if n == edit_nuc:
			if sgrna_strand == 'sense':
				nuc_edit_pos = edit_nuc+str(j+1)
				aa_num,frame = sgrna_trans[nuc_edit_pos].split('_')
				context_index = j+4
			elif sgrna_strand == 'antisense':
				nuc_edit_pos = edit_nuc+str(j)
				aa_num,frame = sgrna_trans[nuc_edit_pos].split('_')
				nuc_edit_pos = revcom(edit_nuc)+str(j)
				context_index = context_index_track

			# Check motif
			if filter_gc:
				proceed = filter_gc_motifs(sgrna_context, context_index, sgrna_strand)
			else:
				proceed = True
			if proceed:
				if frame == '0':
					codon_start = context_index
					codon_end = context_index+3
					old_codon = context[codon_start:codon_end]
				elif frame == '1':
					codon_start = context_index-1
					codon_end = context_index+2
					old_codon = context[codon_start:codon_end]
				elif frame == '2':
					codon_start = context_index-2
					codon_end = context_index+1
					old_codon = context[codon_start:codon_end]
				else:
					old_codon = ''
				if old_codon != '':
					if old_codon in codon_map.keys():
						old_aa = codon_map[old_codon]
						old_aa_3 = list(aa_map.keys())[list(aa_map.values()).index(old_aa)]
						new_codon = []
						edit_indices = []
						for k,x in enumerate(old_codon):
							if x == edit_nuc:
								if filter_gc:
									motif_check = filter_gc_motifs_for_aa(sgrna_strand,sgrna_context,codon_start,k)
								else:
									motif_check = True
								if motif_check:								
									if aa_num+'_'+str(k) in sgrna_trans.values():
										nuc_index = next(key for key, value in sgrna_trans.items() if value == aa_num + '_' + str(k))[1:]
										if (int(nuc_index) >= window_start) and (int(nuc_index) <= window_end):
											new_codon.append(edit_to)
											nuc_index = nuc_index + '-' + str(k)
											edit_indices.append(nuc_index)
										else:
											new_codon.append(x)
									else:
										new_codon.append(x)
								else:
									new_codon.append(x)
							else:
								new_codon.append(x)
					else:
						error = 'Codon '+old_codon+'not standard codon'
						return '', '', error, '', '', '', '', ''
					new_codon = ''.join(new_codon)
					new_aa = codon_map[new_codon]
					new_aa_3 = list(aa_map.keys())[list(aa_map.values()).index(new_aa)]
					#aa_edit = old_aa+str(aa_num)+new_aa
					aa_edit = old_aa_3 + str(aa_num) + new_aa_3
					edit_indices.sort(key=lambda x: int(x.split('-')[0]))
					if sgrna_strand == 'antisense':
						nuc_edit_pos = revcom(edit_nuc)+'_'+'_'.join(edit_indices)
					else:
						nuc_edit_pos = edit_nuc+'_'+'_'.join(edit_indices)
					if nuc_edit_pos not in edit_map.keys():
						if old_aa == new_aa:
							edit_cat = 'Silent'
							num_silent+=1
						elif old_aa != 'Ter' and new_aa == 'Ter':
							edit_cat = 'Nonsense'
							num_stop+=1
						else:
							edit_cat = 'Missense'
						edit_map[nuc_edit_pos] = aa_edit+'_'+edit_cat+'_'+old_codon+'_'+new_codon
				else:
					nuc_edit_pos = nuc_edit_pos[0]+'_'+nuc_edit_pos[1:]+'-'+frame
					if nuc_edit_pos not in edit_map.keys():
						if 'Exon' in frame:
							if (frame.split(':')[1] == '-1') or (frame.split(':')[1] == '-2'):
								edit_cat = 'Splice-acceptor'
							elif (frame.split(':')[1] == '+1') or (frame.split(':')[1] == '+2'):
								edit_cat = 'Splice-donor'
							else:
								edit_cat = 'Intron'
							edit_map[nuc_edit_pos] = frame+'_'+edit_cat
						elif frame == 'U':
							edit_map[nuc_edit_pos] = 'utr_UTR'
						elif frame == 'FS':
							edit_map[nuc_edit_pos] = 'flankseq_Flanking'

		if sgrna_strand == 'sense':
			j+=1
		else:
			j-=1
			context_index_track+=1
		transcript_ref_allele = edit_nuc
		transcript_alt_allele = edit_to
		if gene_strand == 1:
			genome_ref_allele = edit_nuc
			genome_alt_allele = edit_to
		elif gene_strand == -1:
			genome_ref_allele = revcom(edit_nuc)
			genome_alt_allele = revcom(edit_to)
	return edit_map, num_silent, error, num_stop, transcript_ref_allele,transcript_alt_allele,genome_ref_allele,genome_alt_allele


'''
Returns edits for sgRNA, also returns number of silent edits
'''
def get_edit_info(context, sgrna, sgrna_strand, edit, window, pam, sgrna_trans, codon_map, sg_gen_pos, gene_strand, gene_variant_df, aa_map, filter_gc, sgrna_context):
	window_start, window_end = window.split('-')
	sgrna_window = sgrna[int(window_start)-1:int(window_end)]
	j_window = int(window_start)-1
	if sgrna_strand == 'antisense':
		sgrna_window = revcom(sgrna_window)
		j_window = int(window_end)

	edit_map = {}
	edit_map, window_silent, cds_error, num_stop, transcript_ref_allele, transcript_alt_allele, genome_ref_allele, genome_alt_allele = get_edits(edit_map,context,sgrna_window,edit,sgrna_trans,codon_map,j_window,sgrna_strand,pam,int(window_start),int(window_end), aa_map,filter_gc, sgrna_context, gene_strand)
	if cds_error != '':
		# A non-standard codon makes get_edits return '' for edit_map, so there
		# is nothing for get_snps to match against. The caller checks cds_error
		# and writes an error row.
		return edit_map, window_silent, cds_error, num_stop, '', [], transcript_ref_allele, transcript_alt_allele, genome_ref_allele, genome_alt_allele
	snp_type_list, snp_info = get_snps(edit_map, edit, sg_gen_pos, gene_strand, sgrna_strand, gene_variant_df, aa_map)
	clinical_sig = get_clinical_sig(snp_type_list)
	return edit_map, window_silent, cds_error, num_stop, clinical_sig, snp_info, transcript_ref_allele, transcript_alt_allele, genome_ref_allele, genome_alt_allele


'''
Returns edits in format suitable for writing to file
'''
def get_print_edits(edit_map):
	nuc_edits = ''
	aa_edits = ''
	cat = ''
	old_codon = ''
	new_codon = ''
	num_edits = 0
	for k in sorted(edit_map, key= lambda x: int(x.split('_')[1].split('-')[0])):
		v = edit_map[k]
		if '_' in v:
			vals = v.split('_')
			aa_edits = aa_edits + vals[0] + ';'
			cat = cat+vals[1]+';'
			# len(vals) > 2 for coding sequence, <= 2 for non-coding (intron, UTR, flanking)
			if len(vals) > 2:
				old_codon = old_codon + vals[2] + ';'
				new_codon = new_codon + vals[3] + ';'
		else:
			aa_edits = aa_edits + v + ';'
		ne_edits = k.split('_')
		nuc_edits = nuc_edits+ne_edits[0]
		for ne in ne_edits[1:]:
			nuc_edits = nuc_edits +'_' + ne.split('-')[0]
		nuc_edits = nuc_edits + ';'
		num_edits += 1
	return nuc_edits, aa_edits, old_codon, new_codon, cat, num_edits


def get_context_for_trans(ct_index, ct_index_check, abs_pos, pos_for_index, cds_map, fs, sgrna_context, cds_sequence, utr):
	context_for_trans = ''
	flag = 0
	i_count = 0
	cds_pos = ''
	while ct_index < ct_index_check:
		gen_pos = pos_for_index[ct_index]
		if gen_pos in cds_map.keys():
			if flag == 0: #Check to see if ct_index has encountered CDS
				cds_pos = cds_map[gen_pos]
				flag = 1
			cds_str = cds_sequence[cds_map[gen_pos]:cds_map[gen_pos] + (len(sgrna_context) - len(context_for_trans))]
			context_for_trans += cds_str
			ct_index += len(cds_str)
		elif gen_pos in utr.keys():
			context_for_trans += 'U'
			ct_index += 1
		elif gen_pos in fs.keys():
			context_for_trans += 'F'
			ct_index += 1
		else:
			context_for_trans += 'I'
			if flag == 0:
				i_count += 1
			ct_index += 1

	return context_for_trans, cds_pos, i_count


'''
Designs sgRNAs for specified PAM sequence and writes to output    
'''
def design_sgrnas(gene_name, assembly, chromosome, gene_id, designs, gene_seq, abs_pos, pos_for_index, fs, cds_map, utr, t, pam, exons, gene_strand, edit, window, cds_sequence, pam_len, sg_len, cds_start_exon, errors, annotations, gene_variant_df, aa_map, codon_map, input_type, intron_buffer, filter_gc):
	current_exon = cds_start_exon
	for i,e in enumerate(exons):
		e = e.split(':')
		label = 'Exon'+str(current_exon)
		if gene_strand == 1:
			start_pos = abs_pos[int(e[0])]
			pos_end = abs_pos[int(e[1])]
		else:
			start_pos = abs_pos[int(e[1])]
			pos_end = abs_pos[int(e[0])]

		if input_type == 'tid':
			pos = start_pos - intron_buffer
			pos_end = pos_end + intron_buffer
		else:
			pos = start_pos+4

		while pos < pos_end:
			sg_str_anti = 0
			sg_str_sense = 0
			target = gene_seq[pos:(pos+pam_len+sg_len)]
			context = gene_seq[(pos-4):(pos+pam_len+sg_len+4)]
			if len(context) == sg_len+pam_len+8:  #Check if full context is available for target sequence
				m = re.search('[^ATCG]', context)
				if m is None:
					start = target[0:pam_len]
					finish = target[sg_len:sg_len + pam_len]
					m_fwd = re.search(get_pam_pattern(pam), finish)
					m_rev = re.search(get_pam_pattern(revcom(pam)), start)
					if m_rev is not None:
						sg_str_anti = 1
					if m_fwd is not None:
						sg_str_sense = 1
			if sg_str_anti == 1:
				error = ''
				sgrna_for_trans = target[pam_len:pam_len+sg_len]
				sgrna = revcom(sgrna_for_trans)
				res_flag, t4_flag = check_ressite_4t(sgrna)
				sgrna_strand = 'antisense'
				sgrna_pam = revcom(target[0:pam_len])
				sgrna_context = revcom(gene_seq[(pos-3):(pos+pam_len+sg_len+4)])
				sgrna_start_pos = pos+pam_len
				sgrna_end_pos = pos+sg_len+pam_len-1

				ct_index = pos-3
				context_for_trans, cds_pos, i_count = get_context_for_trans(ct_index, pos+sg_len+pam_len+4, abs_pos, pos_for_index, cds_map, fs, sgrna_context, cds_sequence, utr)
				m = re.search('[^F]',context_for_trans)
				if m is None:
					context_for_trans = ''
					error = 'Entire context in flanking sequence'


				if 'I' in context_for_trans:
					map_key_context_start = pos_for_index[pos - 3]
					map_key_context_end = pos_for_index[pos + pam_len + sg_len + 3]
					if map_key_context_start in cds_map.keys():
						# True if sgRNA is in first coding exon or inner exon, false if in last coding exon
						if (cds_map[map_key_context_start] + pam_len + sg_len + 7) <= len(cds_sequence):
							# Get CDS from following exons
							context_for_trans = cds_sequence[cds_map[map_key_context_start]:cds_map[map_key_context_start] + pam_len + sg_len + 7]
						else:
							# Get all following CDS, then fill in U's
							context_for_trans = cds_sequence[cds_map[map_key_context_start]:len(cds_sequence)]
							context_for_trans = context_for_trans + (sg_len + pam_len + 7 - len(context_for_trans))*'U'
					elif map_key_context_end in cds_map.keys():
						# True if sgRNA is in inner exon or last coding exon, false if in first coding exon
						if (cds_map[map_key_context_end] - sg_len - pam_len - 6) >= 0:
							# Get CDS from previous exons
							context_for_trans = cds_sequence[(cds_map[map_key_context_end] - sg_len - pam_len - 6):cds_map[map_key_context_end]+1]
						else:
							# Get all of preceding CDS, then fill in U's
							context_for_trans = cds_sequence[0:cds_map[map_key_context_end]+1] #+1 because ceiling is not included
							context_for_trans = (sg_len + pam_len + 7 - len(context_for_trans))*'U' + context_for_trans
					else:
						if cds_pos != '':
							context_for_trans = cds_sequence[cds_pos-i_count:cds_pos-i_count+len(sgrna_context)]
						else:
							context_for_trans = sgrna_context

				if context_for_trans == '':
					if error == '':
						error = 'No context_for_trans found'
					print(error)
					errors.append([gene_name, t, sgrna, sgrna_strand, error])

				if context_for_trans != '':
					sg_gen_pos = pos_for_index[sgrna_end_pos]
					sgrna_trans = get_sgrna_translated_seq(sgrna_for_trans, cds_map, abs_pos, pos_for_index, fs, sgrna_start_pos, gene_strand, sgrna_strand, utr, e, label)
					edit_map, window_silent, cds_error, num_stop, clinical_sig, snp_info, transcript_ref_allele, transcript_alt_allele, genome_ref_allele, genome_alt_allele = get_edit_info(context_for_trans, sgrna, sgrna_strand, edit, window, pam, sgrna_trans, codon_map, sg_gen_pos, gene_strand, gene_variant_df, aa_map, filter_gc, sgrna_context)
					if cds_error != '':
						errors.append([gene_name, t, sgrna, sgrna_strand, cds_error])
						return 0
					nuc_edits, aa_edits, old_codon, new_codon, cat, num_edits = get_print_edits(edit_map)
					designs.append([sgrna, sgrna_context, gene_name, gene_id,t, gene_strand, assembly, transcript_ref_allele, transcript_alt_allele,
								genome_ref_allele, genome_alt_allele, chromosome, sg_gen_pos, sgrna_strand,
								sgrna_pam, edit, num_edits, window_silent, nuc_edits, aa_edits, cat, clinical_sig, res_flag, t4_flag])
					for snp_row in snp_info:
						annotations.append([sgrna, sgrna_strand, sgrna_context, chromosome, gene_name, gene_strand, edit, transcript_ref_allele,
										 transcript_alt_allele, genome_ref_allele, genome_alt_allele] + snp_row)

			if sg_str_sense == 1:
				error = ''
				sgrna = target[0:sg_len]
				res_flag, t4_flag = check_ressite_4t(sgrna)
				sgrna_strand = 'sense'
				sgrna_pam = target[sg_len:sg_len+pam_len]
				sgrna_context = gene_seq[pos-4:pos+sg_len+pam_len+3]
				sgrna_start_pos = pos
				sgrna_end_pos = pos+sg_len-1

				ct_index = pos-4
				context_for_trans, cds_pos, i_count = get_context_for_trans(ct_index, pos+sg_len+pam_len+3, abs_pos, pos_for_index, cds_map, fs, sgrna_context, cds_sequence, utr)
				m = re.search('[^F]',context_for_trans)
				if m is None:
					context_for_trans = ''
					error = 'Entire context in flanking sequence'


				if 'I' in context_for_trans:
					map_key_context_start = pos_for_index[pos - 4]
					map_key_context_end = pos_for_index[pos + pam_len + sg_len + 3]
					if map_key_context_start in cds_map.keys():
						if (cds_map[map_key_context_start]+sg_len+pam_len+7) <= len(cds_sequence):
							context_for_trans = cds_sequence[cds_map[map_key_context_start]:cds_map[map_key_context_start]+sg_len+pam_len+7]
						else:
							context_for_trans = cds_sequence[cds_map[map_key_context_start]:len(cds_sequence)]
							context_for_trans = context_for_trans + (sg_len + pam_len + 7 - len(context_for_trans))*'U'
					elif map_key_context_end in cds_map.keys():
						if (cds_map[map_key_context_end]-pam_len-sg_len-7) >= 0:
							# Get CDS from previous exons
							context_for_trans = cds_sequence[cds_map[map_key_context_end]-pam_len-sg_len-7:cds_map[map_key_context_end]]
						else:
							context_for_trans = cds_sequence[0:cds_map[map_key_context_end]]
							context_for_trans = (sg_len + pam_len + 7 - len(context_for_trans))*'U' + context_for_trans
					else:
						if cds_pos != '':
							context_for_trans = cds_sequence[cds_pos - i_count:cds_pos - i_count + len(sgrna_context)]
						else:
							context_for_trans = sgrna_context

				if context_for_trans == '':
					if error == '':
						error = 'No context_for_trans found'
					print(error)
					errors.append([gene_name, t, sgrna, sgrna_strand, error])

				if context_for_trans != '':
					sg_gen_pos = pos_for_index[sgrna_start_pos]
					sgrna_trans = get_sgrna_translated_seq(sgrna, cds_map, abs_pos, pos_for_index, fs, sgrna_start_pos, gene_strand, sgrna_strand, utr, e, label)
					edit_map, window_silent, cds_error, num_stop, clinical_sig, snp_info, transcript_ref_allele, transcript_alt_allele, genome_ref_allele, genome_alt_allele = get_edit_info(context_for_trans, sgrna, sgrna_strand, edit, window, pam, sgrna_trans, codon_map, sg_gen_pos, gene_strand, gene_variant_df, aa_map, filter_gc, sgrna_context)
					if cds_error != '':
						errors.append([gene_name, t, sgrna, sgrna_strand, cds_error])
						return 0
					nuc_edits, aa_edits, old_codon, new_codon, cat, num_edits = get_print_edits(edit_map)
					designs.append([sgrna, sgrna_context, gene_name, gene_id,t, gene_strand, assembly, transcript_ref_allele, transcript_alt_allele,
								genome_ref_allele, genome_alt_allele, chromosome, sg_gen_pos, sgrna_strand,
								sgrna_pam, edit, num_edits, window_silent, nuc_edits, aa_edits, cat, clinical_sig, res_flag, t4_flag])
					for snp_row in snp_info:
						annotations.append([sgrna, sgrna_strand, sgrna_context, chromosome, gene_name, gene_strand, edit, transcript_ref_allele,
										 transcript_alt_allele, genome_ref_allele, genome_alt_allele] + snp_row)
			pos = pos+1
		current_exon += 1
	return 1  


def get_seq_info(seq):
	exons = []
	cds = ''
	cur_val = 0
	ex = '0'
	for i,p in enumerate(seq):
		if p.islower() and seq[i-1].isupper():
			ex = ex +':'+ str(i-1)
			exons.append(ex)
			cds = cds + seq[cur_val:i]
		elif p.islower() and seq[i+1].isupper():
			cur_val = i+1
			ex = str(cur_val)
		else:
			continue
	ex = ex +':' +str(i)
	cds += seq[cur_val:]
	exons.append(ex)
	return exons, cds


def check_sequences(seq):
	m = re.search('[^ACTGactg]',seq)
	if m is not None:
		error = 'Sequences contain non-ACTG characters'
	else:
		error = ''
	return error

class PosForIndex(object):
	"""Index -> genomic position, the exact inverse of abs_pos.

	The design loop needs this several times per candidate guide, so it has to
	be cheap. get_absolute_pos walks the gene one base at a time, so index and
	position differ by a constant step and the inverse is arithmetic -- no
	second map, which on a 2 Mb gene like DMD would cost another 84 MB.

	Out-of-range indexes raise KeyError rather than returning a position outside
	the gene, so a lookup that should not happen fails loudly.
	"""
	__slots__ = ('origin', 'step', 'count')

	def __init__(self, abs_pos, gene_strand):
		# get_absolute_pos assigns index 0 first and counts up by one
		self.origin = next(iter(abs_pos))
		self.step = 1 if gene_strand == 1 else -1
		self.count = len(abs_pos)

	def __getitem__(self, index):
		if not 0 <= index < self.count:
			raise KeyError(index)
		return self.origin + self.step * index


def design_transcript(source, clinvar, transcript_id, params=None):
	"""Designs guides for one Ensembl transcript.

	`source` supplies the reference data, `clinvar` the variants (or None to skip
	annotation). Returns (designs, errors, annotations): three lists of rows,
	under DESIGN_COLUMNS, ERROR_COLUMNS and ANNOTATION_COLUMNS.

	A transcript the reference cannot fully answer yields an error row and no
	designs, rather than raising.
	"""
	params = params or DesignParams()
	designs, errors, annotations = [], [], []

	gene_name, assembly, gene_strand, chromosome, gene_id, exons, cds_map, \
		abs_pos_map, fs, utr, cds_start_exon, utr5_flag, utr3_flag = \
		get_tr_info(source, transcript_id, 'tid')

	def failed(error):
		errors.append([gene_name, transcript_id, 'N/A', 'N/A', error])
		return designs, errors, annotations

	if exons == '':
		return failed('Transcript not found')

	tr_seq = get_tr_sequence(source, transcript_id)
	seq_error = check_sequences(tr_seq)
	if seq_error != '':
		return failed(seq_error)
	if tr_seq == '':
		return failed('Transcript sequence not found')

	if get_pro_sequence(source, transcript_id) == '':
		return failed('Protein sequence not found')

	cds_sequence = get_cds_sequence(source, transcript_id)
	seq_error = check_sequences(cds_sequence)
	if seq_error != '':
		return failed(seq_error)
	if cds_sequence == '':
		return failed('Coding sequence not found')

	gene_variant_df = (clinvar.variants_for_gene(gene_name) if clinvar
					   else empty_variants())

	_run_edits(designs, errors, annotations, params,
			   gene_name=gene_name, assembly=assembly, chromosome=chromosome,
			   gene_id=gene_id, gene_seq=tr_seq, abs_pos=abs_pos_map, fs=fs,
			   cds_map=cds_map, utr=utr, t=transcript_id, exons=exons,
			   gene_strand=gene_strand, cds_sequence=cds_sequence,
			   cds_start_exon=cds_start_exon, gene_variant_df=gene_variant_df,
			   input_type='tid')
	return designs, errors, annotations


def design_sequence(name, sequence, params=None):
	"""Designs guides for a raw nucleotide sequence, the FASTA input path.

	Exons are taken from the case of the sequence, as get_seq_info describes.
	There is no gene to look up, so no ClinVar annotation and no `source`.
	"""
	params = params or DesignParams()
	designs, errors, annotations = [], [], []

	seq_error = check_sequences(sequence)
	if seq_error != '':
		errors.append([name, name, 'N/A', 'N/A', seq_error])
		return designs, errors, annotations

	exons, cds_sequence = get_seq_info(sequence)
	abs_pos_map, fs = get_absolute_pos(0, len(sequence), 1, 'nuc')
	cds_map = get_cds_map(exons, 1)

	_run_edits(designs, errors, annotations, params,
			   gene_name=name, assembly='', chromosome='', gene_id='',
			   gene_seq=sequence, abs_pos=abs_pos_map, fs=fs, cds_map=cds_map,
			   utr={}, t=name, exons=exons, gene_strand=1,
			   cds_sequence=cds_sequence, cds_start_exon=1,
			   gene_variant_df=empty_variants(), input_type='nuc')
	return designs, errors, annotations


def _run_edits(designs, errors, annotations, params, **shared):
	"""Runs design_sgrnas once per deaminase pass, into the same three lists.

	`shared` is everything about the transcript that both passes reuse; its keys
	are design_sgrnas parameter names, so it forwards straight through.
	"""
	pos_for_index = PosForIndex(shared['abs_pos'], shared['gene_strand'])
	codon_map = get_codon_map()
	aa_map = get_aa_map()
	for edit in params.edits:
		design_sgrnas(
			designs=designs, errors=errors, annotations=annotations,
			pos_for_index=pos_for_index, aa_map=aa_map, codon_map=codon_map,
			edit=edit, pam=params.pam, pam_len=len(params.pam),
			window=params.window, sg_len=params.sg_len,
			intron_buffer=params.intron_buffer, filter_gc=params.filter_gc,
			**shared)
