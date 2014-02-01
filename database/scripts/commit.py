# Commit data from database-stage.xml to database.xml
#
# Add datestamps and versionIDs
#
# Daniel L. Parton <partond@mskcc.org> - 15 Nov 2013
#

#==============================================================================
# Imports
#==============================================================================

import sys, datetime, os, re
import choderalab as clab
from lxml import etree

#==============================================================================
# Parameters
#==============================================================================

if '-write' in sys.argv:
    run_mode = 'write'
else:
    run_mode = 'nowrite'

print 'Running in mode: %s' % run_mode

DB_dir = 'database'

DBstage_filename = 'database-stage.xml'
DB_filename = 'database.xml'
DBstage_filepath = os.path.join(DB_dir, DBstage_filename)
DB_filepath = os.path.join(DB_dir, DB_filename)

now = datetime.datetime.utcnow()
now_datestamp = now.strftime(clab.DB.datestamp_format_string)

#==============================================================================
# Read main and staged versions of database
#==============================================================================

parser = etree.XMLParser(remove_blank_text=True)

# We need database-stage.xml, so fail if this is not present
if not os.path.exists(DBstage_filepath):
    raise Exception, '%s not found' % DBstage_filepath
DBstage_root = etree.parse(DBstage_filepath, parser).getroot()

#==============================================================================
# Check whether the staged database is ready to be committed
#==============================================================================

DBstage_datestamps_dict = { key : datetime.datetime.strptime(DBstage_root.attrib[key], clab.DB.datestamp_format_string) for key in DBstage_root.keys() if re.match('.*_last_.*', key) }

# make sure that all scripts have been run since the last time the DB was modified
for gather_script_ID in clab.DB.gather_script_IDs:
    if clab.DB.gather_script_last_run_key_dict[gather_script_ID] not in DBstage_root.keys():
        print 'Script %s has not been run yet.' % gather_script_ID
        print 'Exiting.'
        sys.exit()

DBstage_gather_uniprot_last_run = DBstage_datestamps_dict['gather_uniprot_last_run']

# compare each of the *_last_run attribs within DBstage against gather_uniprot_last_run, to make sure that gather-uniprot was run before all others
for gather_script_ID in clab.DB.gather_script_IDs:
    gather_script_last_run_key = clab.DB.gather_script_last_run_key_dict[gather_script_ID]
    if DBstage_datestamps_dict[gather_script_last_run_key] < DBstage_gather_uniprot_last_run:
        print '"%s" is dated "%s", while "gather_uniprot_last_run" is dated "%s". gather-uniprot.py must be the earliest script staged, so please run and stage the other script before committing.' % (gather_script_last_run_key, DBstage_datestamps_dict[gather_script_last_run_key], DBstage_gather_uniprot_last_run)
        print 'Exiting.'
        sys.exit()

#==============================================================================
# If the main database doesn't exist yet, then it will be created from database-stage.xml. Script then exits.
#==============================================================================
    
if not os.path.exists(DB_filepath):
    print '%(DB_filename)s not found' % vars()
    DB_root = DBstage_root
    DB_root.set('versionID', '0')
    DB_root.set('commit_last_run', now_datestamp)
    DB_root.set('commit_last_modif', now_datestamp)
    if run_mode == 'write':
        print 'Creating %(DB_filename)s from %(DBstage_filename)s...' % vars()
        clab.DB.writeDB(DB_root, DB_filepath)
    print 'Done.'
    sys.exit()

# If the database.xml does already exist, then parse it and get the versionID
else:
    DB_root = etree.parse(DB_filepath, parser).getroot()
    versionID = int(DB_root.get('versionID'))

#==============================================================================
# Now that we know an existing version of the DB is present, we compare all *_last_run attribs between DBstage and DB, to make sure all scripts have been run since the last commit
#==============================================================================

DB_datestamps_dict = { key : datetime.datetime.strptime(DB_root.attrib[key], clab.DB.datestamp_format_string) for key in DB_root.keys() if re.match('.*_last_.*', key) }

for gather_script_ID in clab.DB.gather_script_IDs:
    gather_script_last_run_key = clab.DB.gather_script_last_run_key_dict[gather_script_ID]
    if gather_script_last_run_key not in DB_datestamps_dict.keys():
        continue
    if DBstage_datestamps_dict[gather_script_last_run_key] < DB_datestamps_dict[gather_script_last_run_key]:
        print 'DBstage "%s" is dated "%s", while DB "%s" is dated "%s". All scripts must have been run at least once since the last commit.' % (gather_script_last_run_key, DBstage_datestamps_dict[gather_script_last_run_key], gather_script_last_run_key, DB_datestamps_dict[gather_script_last_run_key])
        print 'Exiting.'
        sys.exit()

#==============================================================================
# Now we compare the contents of DB and DBstage
#==============================================================================

# First a quick check of the number of entries

nDB_entries = len(DB_root.findall('entry'))
nDBstage_entries = len(DBstage_root.findall('entry'))
diff_nentries = nDBstage_entries - nDB_entries

if diff_nentries != 0:
    print '%(DBstage_filename)s has %(diff_nentries)+d entries compared to %(DB_filename)s' % vars()
    data_modified = True

# If the numbers of entries match, move on to a full comparison using GNU diff (ignoring the first line (datestamps))

else:
    DB_comparison_string = etree.tostring(DB_root, pretty_print=True)
    DBstage_comparison_string = etree.tostring(DBstage_root, pretty_print=True)
    DB_comparison_string = '\n'.join( DB_comparison_string.splitlines()[1:] ) # Ignore the first line, which contains datestamps (don't want to include these in the comparison)
    DBstage_comparison_string = '\n'.join( DBstage_comparison_string.splitlines()[1:] ) # Ignore the first line, which contains datestamps (don't want to include these in the comparison)

    diff_output = clab.DB.diff_DB_comparison_strings(DB_comparison_string, DBstage_comparison_string)

    if len(diff_output) > 0:
        print 'Comparison of data in database and database-stage indicates changes. Lines in diff comparison: %s' % len(diff_output)
        data_modified = True
    else:
        print 'Comparison of data in database and database-stage indicates no changes. Will update gather_uniprot_last_run attrib, but other data will not be modified.'
        data_modified = False

#==============================================================================
# Update versionids and datestamps
#==============================================================================

DBstage_root.set('commit_last_run', now_datestamp)

if data_modified:
    versionID += 1
    DBstage_root.set('commit_last_modif', now_datestamp)
else:
    DBstage_root.set('commit_last_modif', DB_root.get('commit_last_modif'))

DBstage_root.set('versionID', str(versionID))

#==============================================================================
# Write database to file
#==============================================================================

if run_mode == 'write':
    print 'Writing %(DB_filename)s...' % vars()
    clab.DB.writeDB(DBstage_root, DB_filepath)

print 'Done.'

