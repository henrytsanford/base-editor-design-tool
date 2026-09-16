import sys

# NB: keep this block parseable by Python 2 so it can report the version
# instead of dying with a SyntaxError. No f-strings here.
if sys.version_info < (3, 9):
	sys.exit(
		"This tool requires Python 3.9 or newer (tested on 3.13); you are "
		"running %s.\nSee the README for setup instructions."
		% sys.version.split()[0]
	)


import argparse
import contextlib
import csv
import os
from datetime import datetime

import pandas as pd
from Bio import SeqIO

from bedesign import (
	ANNOTATION_COLUMNS,
	ANNOTATIONS_FILE,
	DESIGN_COLUMNS,
	DESIGNS_FILE,
	ERROR_COLUMNS,
	ERRORS_FILE,
	DesignParams,
	UnknownBaseEditor,
	design_sequence,
	design_transcript,
	strip_tr_version,
)
# Transcript reference data, from the Ensembl REST API or a local bundle.
from bedesign.transcript_source import (
	BundleNotFound,
	ClinVarNotFound,
	ClinVarSource,
	EnsemblUnavailable,
	TranscriptNotFound,
	find_clinvar_db,
	get_source,
)


class InputFileError(Exception):
	"""--input-file could not be read. Reported by the CLI, not raised by the engine."""


DEFAULTS = DesignParams()


def get_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('--input-file',
		type=str,
		help='File with Ensembl transcript IDs or fasta file with nucleotide sequences')
	parser.add_argument('--clinvar-db',
		type=str,
		default=None,
		help='ClinVar database built by tools/build_clinvar.py '
			 '(default: the newest clinvar-<date>.db under --refdata)')
	parser.add_argument('--no-clinvar',
		action='store_true',
		help='Skip ClinVar annotation')
	# Accepted only so that a run using it gets a pointer instead of argparse's
	# 'unrecognized arguments'. Handled in __main__.
	parser.add_argument('--variant-file',
		type=str,
		default=None,
		help=argparse.SUPPRESS)
	parser.add_argument('--input-type',
		type=str,
		help='tid for transcript IDs and nuc for nucleotide sequence')
	parser.add_argument('--be-type',
		type=str,
		default='',
		help='Type of base editor')
	# Defaults come from DesignParams so the CLI and the importable engine cannot
	# drift apart; only --be-type and the output options are the CLI's own.
	parser.add_argument('--pam',
		type=str,
		default=DEFAULTS.pam,
		help='PAM sequence for guide design')
	parser.add_argument('--edit-window',
		type=str,
		default=DEFAULTS.window,
		help='Editing window')
	parser.add_argument('--sg-len',
		type=int,
		default=DEFAULTS.sg_len,
		help='Length of sgRNA')
	parser.add_argument('--edit',
		type=str,
		default=DEFAULTS.edit,
		help='Edit')
	parser.add_argument('--intron-buffer',
		type=int,
		default=DEFAULTS.intron_buffer,
		help='How far to tile into introns (bp)')
	parser.add_argument('--filter-gc',
		type=str,
		choices=['True','False'],
		default=str(DEFAULTS.filter_gc),
		help='Whether to filter out edits in a GC motif')
	parser.add_argument('--output-name',
		type=str,
		help='Output name')
	parser.add_argument('--source',
		type=str,
		choices=['rest','local'],
		default='rest',
		help='Where transcript reference data comes from: the Ensembl REST API (rest) '
			 'or a local reference bundle (local)')
	parser.add_argument('--refdata',
		type=str,
		default='refdata',
		help='Directory holding the local reference bundle (used with --source local)')
	return parser

'''
Turns a reference-data failure into a readable message instead of a traceback.
Used as a context manager so the run's output files are closed before exiting.
Only the CLI exits: the exceptions themselves stay catchable by other callers.
'''
@contextlib.contextmanager
def reference_error_guard(output_folder):
	try:
		yield
	except TranscriptNotFound as e:
		sys.exit("%s\nPartial output remains in %s." % (e, output_folder))
	except EnsemblUnavailable as e:
		sys.exit(
			"%s\n"
			"The Ensembl REST API is not responding; this is usually a temporary "
			"outage on their end. Re-run with --source local to design from a local "
			"reference bundle instead. Partial output remains in %s." % (e, output_folder)
		)

def write_readme(output_folder, input_file, params, clinvar, source):
	with open(output_folder + '/README.txt', 'w') as o:
		w = csv.writer(o, delimiter='\t')
		w.writerow(['Input file: ' + input_file])
		w.writerow(['PAM: ' + params.pam])
		w.writerow(['Edit window: ' + params.window])
		w.writerow(['Edit: ' + params.edit])
		w.writerow(['Intron Buffer: ' + str(params.intron_buffer)])
		w.writerow(['Filter out GC motifs: ' + str(params.filter_gc)])
		w.writerow(['Reference: ' + source.describe()])
		w.writerow(['Variants: ' + (clinvar.describe() if clinvar else 'none')])
		w.writerow(['Output folder: ' + output_folder])


def read_args(args):
	"""Turns the command line into what the engine and the output folder need.

	Raises UnknownBaseEditor for an unrecognised --be-type, and ClinVarNotFound
	when ClinVar is wanted but no database can be resolved.
	"""
	input_type = args.input_type
	if input_type == 'tid':
		try:
			input_df = pd.read_table(args.input_file)
		except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as e:
			# Named separately from the engine's errors so an unreadable input
			# file says which file, rather than surfacing as a bare pandas message.
			raise InputFileError(
				"Could not read --input-file '%s': %s\n"
				"It should be a tab-separated file with a header row, then one "
				"transcript ID and gene symbol per line." % (args.input_file, e))
		input_df.iloc[:, 0] = input_df.iloc[:, 0].map(strip_tr_version)
	else:
		input_df = pd.DataFrame(columns=['Sequence', 'ID'])
		with open(args.input_file) as fh:
			for i, fasta in enumerate(SeqIO.parse(fh, 'fasta')):
				input_df.loc[i, 'ID'] = fasta.id
				input_df.loc[i, 'Sequence'] = str(fasta.seq)
	overrides = dict(intron_buffer=args.intron_buffer,
					 filter_gc=args.filter_gc == 'True')
	if args.be_type != '':
		params = DesignParams.from_preset(args.be_type, **overrides)
	else:
		params = DesignParams(pam=args.pam, window=args.edit_window,
							  sg_len=args.sg_len, edit=args.edit, **overrides)
	# Nucleotide input has no gene to look up, so it never needs the database open.
	if args.no_clinvar or input_type != 'tid':
		clinvar = None
	else:
		clinvar = ClinVarSource(find_clinvar_db(args.clinvar_db, args.refdata))
	output_folder = args.output_name + '_' + str(datetime.now().strftime("%y-%m-%d-%H-%M-%S"))
	if not os.path.exists(output_folder):
		os.makedirs(output_folder)
	return input_df, params, clinvar, output_folder


def main():
	args = get_parser().parse_args()
	if args.variant_file is not None:
		sys.exit(
			"--variant-file is not used any more. ClinVar SNPs are read from a "
			"database instead of variant_summary.txt.\n"
			"Build one with: python tools/build_clinvar.py --variant-summary %s\n"
			"Then re-run without --variant-file, or pass --clinvar-db <path>."
			% args.variant_file)
	if args.input_type not in ('tid', 'nuc'):
		sys.exit('Please enter a valid input type; tid for a list of transcripts '
				 'or nuc for fasta sequences')
	try:
		source = get_source(args.source, args.refdata)
		input_df, params, clinvar, output_folder = read_args(args)
	except (BundleNotFound, ClinVarNotFound, UnknownBaseEditor, InputFileError) as e:
		sys.exit(str(e))
	write_readme(output_folder, args.input_file, params, clinvar, source)

	def opened(stack, filename, columns):
		path = os.path.join(output_folder, filename % args.output_name
							if '%s' in filename else filename)
		writer = csv.writer(stack.enter_context(open(path, 'w')), delimiter='\t')
		writer.writerow(columns)
		return writer

	with contextlib.ExitStack() as stack:
		stack.enter_context(reference_error_guard(output_folder))
		w = opened(stack, DESIGNS_FILE, DESIGN_COLUMNS)
		w_error = opened(stack, ERRORS_FILE, ERROR_COLUMNS)
		# Nucleotide input has no gene to look up, so it never annotates and the
		# file is not created at all.
		w_clin = (opened(stack, ANNOTATIONS_FILE, ANNOTATION_COLUMNS)
				  if args.input_type == 'tid' else None)
		for _, r in input_df.iterrows():
			print('Designing for ' + r.iloc[1])
			if args.input_type == 'tid':
				designs, errors, annotations = design_transcript(
					source, clinvar, r.iloc[0], params)
			else:
				designs, errors, annotations = design_sequence(
					r.iloc[1], r.iloc[0], params)
			w.writerows(designs)
			w_error.writerows(errors)
			if w_clin is not None:
				w_clin.writerows(annotations)


if __name__ == '__main__':
	main()
