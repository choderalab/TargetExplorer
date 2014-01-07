import os, textwrap
#import choderalab
from lxml.builder import E

# Look for the kinome root directory, which should be two below the pylib/choderalab directory
# XXX DEPRECATED, since we will eventually move to a setuptools-based installation, in which case the Python libraries will be separated from the database directories
#two_dirs_below_choderalab = os.path.join(choderalab.__path__[0], '..', '..')
#if os.path.exists(two_dirs_below_choderalab):
#    kinome_rootdir = os.path.abspath( two_dirs_below_choderalab )
#else:
#    kinome_rootdir = None

def check_correct_working_dir():
    if not os.path.exists('database') or not os.path.exists('external-data') or not os.path.exists('analysis'):
        raise Exception, 'Working directory must be the root directory of the TargetExplorer installation.'

def seq2pretty_html(seq, aa_css_class_list=None):
    '''Pass a sequence string.
    Returns a list of lxml html td elements, colored according to residue type.
    Use the aa_css_class_list option to pass a list of custom aa_css_classes. Must be the same length as the seq list.
    '''
    if aa_css_class_list != None and len(aa_css_class_list) != len(seq):
        raise Exception, 'aa_css_class_list must be list of same length as seq list.'

    spans = []

    for i, aa in enumerate(seq):
        styled_aa = E.span(aa)
        aatype = aa_types[aa]
        if aa_css_class_list == None:
            aacolor = aa_css_classes[aatype]
        else:
            aacolor = aa_css_class_list[i]
        styled_aa.set('class','%s' % aacolor)
        spans.append( styled_aa )
    return spans

aa_css_classes = {
'A':'c1', # aromatic
'C':'c2', # cysteine
'H':'c0', # hydrophobic
'+':'c5', # positive
'-':'c4', # negative
'P':'c7', # polar
'0':'gr', # gap
'x':'bl'}

aa_types = {
'A':'H', # hydrophobic
'C':'C', # cysteine
'D':'-', # negative
'E':'-',
'F':'A', # aromatic
'G':'P', # polar
'H':'+', # positive
'I':'H',
'K':'+',
'L':'H',
'M':'H',
'N':'P',
'P':'H',
'Q':'P',
'R':'+',
'S':'P',
'T':'P',
'V':'H',
'W':'A',
'Y':'A',
'X':'x',
'-':'0', # gap
'x':'x', # UNK (unknown) - present in 3LZB
'a':'x', # lower case represents conflicting residues
'c':'x',
'd':'x',
'e':'x',
'f':'x',
'g':'x',
'h':'x',
'i':'x',
'k':'x',
'l':'x',
'm':'x',
'n':'x',
'p':'x',
'q':'x',
'r':'x',
's':'x',
't':'x',
'v':'x',
'w':'x',
'y':'x'}


def which(program):
    '''
    Searches $PATH environment variable (from parent shell) for a program.
    '''
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def parse_fasta_file(filepath):
    '''Should also work with Vienna format
'''
    with open(filepath, 'r') as fasta_file:
        seq_strings = fasta_file.read().split('>')[1:] # First element is '', so ignore
        seq_strings_split = [ seq_string.split('\n') for seq_string in seq_strings ]
        seq_ids = [ seq_string_lines[0] for seq_string_lines in seq_strings_split ]
        sequences = [ ''.join(seq_string_lines[1:]) for seq_string_lines in seq_strings_split ]
    return seq_ids, sequences

def parse_fasta_string(fasta_string):
    '''Should also work with Vienna format
'''
    seq_strings = fasta_string.split('>')[1:] # First element is '', so ignore
    seq_strings_split = [ seq_string.split('\n') for seq_string in seq_strings ]
    seq_ids = [ seq_string_lines[0] for seq_string_lines in seq_strings_split ]
    sequences = [ ''.join(seq_string_lines[1:]) for seq_string_lines in seq_strings_split ]
    return seq_ids, sequences

def sequnwrap(sequence):
    '''
    Unwraps a wrapped sequence string
    '''
    unwrapped = sequence.strip()
    unwrapped = ''.join(unwrapped.split('\n'))
    return unwrapped

def seqwrap(sequence, add_star=False):
    '''
    Wraps a sequence string to a width of 60.
    If add_star is set to true, an asterisk will be added
    to the end of the sequence, for compatibility with
    Modeller.
    '''
    if add_star:
        sequence += '*'
    wrapped = ''
    for i in range(0,len(sequence),60):
        wrapped += sequence[i: i+60] + '\n'
    return wrapped

def twrap(text):
    '''
    Wraps text to a width of 60
    '''
    return '\n' + textwrap.fill(text, width=60) + '\n'

def match_kinDB_ID(stringtomatch):
    import re
    matchpattern = '.+_.+_.+_PK[0-9]+'
    return re.match(matchpattern, stringtomatch)

def parse_target_args(arglist, revert_to_all_targets=False):
    '''
    Pass sys.argv
    Possible arguments for -target flag:
        "SRC_HUMAN_P12931_PK0"
        '["SRC_HUMAN_P12931_PK0", "ABL1_HUMAN_P00519_PK0"]'
    Returns a list in either case, as follows:
        ["SRC_HUMAN_P12931_PK0"]
        ["SRC_HUMAN_P12931_PK0", "ABL1_HUMAN_P00519_PK0"]
    '''
    # If -targets flag not found
    if '-targets' not in arglist:
        if revert_to_all_targets:
            return 'all_targets'
        else:
            raise Exception, '-targets flag is required'
    try:
        targets = arglist[ arglist.index('-targets') + 1 ]
        if targets[0] == '[' and targets[-1] == ']':
            from ast import literal_eval
            return literal_eval(targets)
        elif match_kinDB_ID(targets):
            return [targets]
        else:
            raise Exception, '-targets argument set to unrecognizable string'
    except ValueError:
        raise Exception, 'no argument after -targets flag.'

def count_lines_in_file(filepath):
    import subprocess
    wc_output = subprocess.check_output(['wc', '-l', filepath])
    nlines = int(wc_output.split()[0])
    return nlines

