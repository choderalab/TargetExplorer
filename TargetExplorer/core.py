import os, textwrap
#import choderalab
from lxml.builder import E

# ========
# Global package variables
# ========

datestamp_format_string = '%Y-%m-%d %H:%M:%S UTC'

# =========
# =========


def xpath_match_regex_case_sensitive(context, attrib_values, xpath_argument):
    ''' To be used as an lxml XPath extension, for regex searches of attrib values.
    '''
    import re
    # If no attrib found
    if len(attrib_values) == 0:
        return False
    # If attrib found, then run match against regex
    else:
        regex = re.compile(xpath_argument)
        return bool( re.search(regex, attrib_values[0]) )


def xpath_match_regex_case_insensitive(context, attrib_values, xpath_argument):
    ''' To be used as an lxml XPath extension, for regex searches of attrib values.
    '''
    import re
    # If no attrib found
    if len(attrib_values) == 0:
        return False
    # If attrib found, then run match against regex
    else:
        regex = re.compile(xpath_argument, re.IGNORECASE)
        return bool( re.search(regex, attrib_values[0]) )


def int_else_none(literal):
    try:
        return int(literal)
    except ValueError:
        return None








#############
# XXX OLD XXX
#############

# Look for the kinome root directory, which should be two below the pylib/choderalab directory
# XXX kinome_rootdir is DEPRECATED, since we will eventually move to a setuptools-based installation, in which case the Python libraries will be separated from the database directories
#two_dirs_below_choderalab = os.path.join(choderalab.__path__[0], '..', '..')
#if os.path.exists(two_dirs_below_choderalab):
#    kinome_rootdir = os.path.abspath( two_dirs_below_choderalab )
#else:
#    kinome_rootdir = None

def parse_nested_dicts(nested_dict, parse_keys):
    '''
    Pass a nested dict structure and a list of keys to parse hierarchically.
    '''
    if not isinstance(nested_dict, dict):
        raise Exception, 'parse_nested_dicts expected to be passed a dict.'

    child = nested_dict

    for parse_key in parse_keys:
        #current_dict = get_from_dict(current_dict, parse_key)
        child = child.get(parse_key)
        if child == None:
            return None
        elif parse_key == parse_keys[-1]:
            return child
        elif not isinstance(child, dict):
            return None

    #def get_from_dict(input_dict, parse_key):
    #    return result

def check_correct_working_dir():
    if not os.path.exists('database') or not os.path.exists('external-data') or not os.path.exists('analysis'):
        raise Exception, 'Working directory must be the root directory of the TargetExplorer installation.'

def seq2pretty_html(seq, aa_css_class_list=None):
    '''Pass a sequence string.
    Returns a list of lxml html spans elements, colored according to residue type.
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

def write_css_stylesheet(filepath):
    '''
    Write a CSS stylesheet containing classes with custom colors for displaying alignments.
    '''
    css_text = '''.ali {
    font-family:Lucida Console,monospace;
    font-size:16px;
    font-weight:normal;
    line-height:97%;
}
.tblheader {
    font-family:Trebuchet MS, Verdana;
    letter-spacing: -0.06em;
    line-height: 96%; 
    font-size:14px;
    font-weight:normal;
}
.tc0 { color: #3e3f61; margin-right:-0.5em }
.tc1 { color: #585989 }
.tc2 { color: #7792ba }
.tc3 { color: #a1b4cc }
.tc4 { color: #b1c4dc }
.tc0bg { background: #3e3f61 }
.tc1bg { background: #585989 }
.tc2bg { background: #7792ba }
.tc3bg { background: #a1b4cc }
.tc4bg { background: #d7e4f7 }

circle.node {
    cursor: pointer;
    stroke: #3182bd;
    stroke-width: 1.5px;
}
line.link {
    fill: none;
    stroke: #9ecae1;
    stroke-width: 1.5px;
}

/* mview / multalin / tailor coloring */
.gr  { color:grey;    }
.bl  { color:black;   }
.m0  { color:blue;    }
.m1  { color:red;     }
.c0  { color:#33cc00; }
.c1  { color:#009900; }
.c2  { color:#ffff00; }
.c3  { color:#33cc00; }
.c4  { color:#cc0000; }
.c5  { color:#0033ff; }
.c6  { color:#6600cc; }
.c7  { color:#0099ff; }
.c8  { color:#666666; }
.c9  { color:#999999; }
.t0  { color:#5858a7; }
.t1  { color:#6b6b94; }
.t2  { color:#64649b; }
.t3  { color:#2121de; }
.t4  { color:#9d9d62; }
.t5  { color:#8c8c73; }
.t6  { color:#0000ff; }
.t7  { color:#4949b6; }
.t8  { color:#60609f; }
.t9  { color:#ecec13; }
.t10 { color:#b2b24d; }
.t11 { color:#4747b8; }
.t12 { color:#82827d; }
.t13 { color:#c2c23d; }
.t14 { color:#2323dc; }
.t15 { color:#4949b6; }
.t16 { color:#9d9d62; }
.t17 { color:#c0c03f; }
.t18 { color:#d3d32c; }
.t19 { color:#ffff00; }

/* charged - positive */
.pc1 { background-color:red; }
/* charged - negative */
.pc2 { background-color:blue; }
/* polar- uncharged */
.pc3 { background-color:red; }
/* special cases */
.pc4 { background-color:yellow; }
/* hydrophobic */
.pc5 { background-color:green; }
'''
    with open(filepath, 'w') as fileobj:
        fileobj.write(css_text)

