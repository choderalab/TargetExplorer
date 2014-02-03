import os, subprocess, shutil, tempfile, datetime
from lxml import etree
import TargetExplorer

datestamp_format_string = '%Y-%m-%d %H:%M:%S UTC'

DB_script_IDs = ['UniProt', 'PDB', 'NCBI_Gene', 'BindingDB', 'cBioPortal', 'prioritization']

DB_script_last_run_attribs = {
'UniProt' : 'gather_uniprot_last_run',
'PDB' : 'gather_pdb_last_run',
'NCBI_Gene' : 'gather_ncbi_gene_last_run',
'BindingDB' : 'gather_bindingdb_last_run',
'cBioPortal' : 'gather_cbioportal_last_run',
'prioritization' : 'prioritization_last_run'
}
DB_script_last_modif_attribs = {
'UniProt' : 'gather_uniprot_last_run',
'PDB' : 'gather_pdb_last_run',
'NCBI_Gene' : 'gather_ncbi_gene_last_run',
'BindingDB' : 'gather_bindingdb_last_run',
'cBioPortal' : 'gather_cbioportal_last_run',
'prioritization' : 'prioritization_last_run'
}

external_data_dir = os.path.join('external-data')
external_data_metadata_filepath = os.path.join(external_data_dir, 'metadata.xml')

local_data_filenames = {
'UniProt' : 'uniprot-search.xml' ,
'NCBI_Gene' : 'gene2pubmed.gz' ,
'BindingDB' : 'BindingDB_All.tab' ,
'cBioPortal' : 'cbioportal-mutations.xml'
}

local_data_filepaths = { DB_script_ID : os.path.join(external_data_dir, DB_script_ID, local_data_filenames[DB_script_ID])  for DB_script_ID in local_data_filenames.keys() }

retrieve_methods = {
'BindingDB' : TargetExplorer.BindingDB.retrieve_all_BindingDB_data
}

days_elapsed_for_suggestdl_dict = {
'UniProt' : 7 ,
'NCBI_Gene' : 7 ,
'BindingDB' : 30 ,
'cBioPortal' : 7
}


def update_external_data_metadata(DB_script_ID, datestamp, filename, filepath):
    '''
    To be called by method retrieve_external_data
    '''
    parser = etree.XMLParser(remove_blank_text=True)
    metadata_root = etree.parse(external_data_metadata_filepath, parser).getroot()
    DB_script_type_node = metadata_root.find(DB_script_ID)
    if DB_script_type_node == None:
        DB_script_type_node = etree.SubElement(metadata_root, DB_script_ID)
    xml_query = 'local_data_file[@filename="%s"]' % filename
    local_data_file_node = DB_script_type_node.find(xml_query)
    if local_data_file_node == None:
        local_data_file_node = etree.SubElement(DB_script_type_node, 'local_data_file')
    local_data_file_node.set('filename', filename)
    local_data_file_node.set('filepath', filepath)
    local_data_file_node.set('datestamp', datestamp)

    with open(external_data_metadata_filepath, 'w') as metadata_file:
        metadata_file.write(etree.tostring(metadata_root, pretty_print=True))


def retrieve_external_data(DB_script_ID, forcedl=False, uniprot_search_string=None):
    # TODO currently only works with gather-BindingDB
    parser = etree.XMLParser(remove_blank_text=True)
    now = datetime.datetime.utcnow()
    now_datestamp = now.strftime(datestamp_format_string)

    local_data_filename = local_data_filenames[DB_script_ID]
    local_data_filepath = local_data_filepaths[DB_script_ID]
    # TODO decorate this (with a decorator selected from a dict of decorators keyed by DB_script_ID) to add any additional required args (e.g. uniprot_search_string)
    retrieve_new_data_file = retrieve_methods[DB_script_ID]
    days_elapsed_for_suggestdl = days_elapsed_for_suggestdl_dict[DB_script_ID]

    # First check if local data file already exists
    if os.path.exists(local_data_filepath):
        print 'Local data file found at:', local_data_filepath
    # If not, download it
    else:
        print 'Local data file not found.'
        print 'Retrieving new data file from external server...'
        retrieve_new_data_file(local_data_filepath)
        update_external_data_metadata(DB_script_ID, now_datestamp, local_data_filename, local_data_filepath)

    # Check when the local data file was retrieved and download a new one if it is older than days_elapsed_for_suggestdl
    external_data_metadata_root = etree.parse(TargetExplorer.DB.external_data_metadata_filepath, parser).getroot()
    xml_query = '%s/local_data_file[@filename="%s"]' % (DB_script_ID, local_data_filename)
    local_datafile_datestamp = external_data_metadata_root.find(xml_query).get('datestamp')
    local_datafile_datestamp = datetime.datetime.strptime(local_datafile_datestamp, TargetExplorer.DB.datestamp_format_string)
    time_elapsed = now - local_datafile_datestamp
    if (time_elapsed.days > days_elapsed_for_suggestdl) or (forcedl == True):
        if time_elapsed.days > days_elapsed_for_suggestdl:
            print 'Local data file is %d days old.' % time_elapsed.days
        if forcedl:
            print 'Forcing retrieval of new data file.'
            download_new_local_data_file = True
        else:
            while True:
                user_response = raw_input('Suggest retrieving new data file from external server. Proceed? [y] ')
                if user_response in ['y', '']:
                    download_new_local_data_file = True
                    break
                elif user_response == 'n':
                    download_new_local_data_file = False
                    break
                else:
                    print 'User input not understood. Please try again.'

        if download_new_local_data_file:
            print 'Retrieving new data file from external server...'
            retrieve_new_data_file(local_data_filepath)
            update_external_data_metadata(DB_script_ID, now_datestamp, local_data_filename, local_data_filepath)

    print ''


def format_target_info(target_info_dict=None, include_pseudogene_data=False, return_format_string=False):
    if return_format_string:
        if include_pseudogene_data:
            format_string = '{:<25}  {:<12}    {:<10}  {:<9}  {:<6}    {:<8}  {:<8}  {:<15}  {:<10}  {:<10}'.format('kinDB_id', 'target_score', 'family', 'pk_length', 'npubs', '%muts', 'npk_pdbs', 'ndisease_assocs', 'nbioassays', 'pseudogene')
        else:
            format_string = '{:<25}  {:<12}    {:<10}  {:<9}  {:<6}    {:<8}  {:<8}  {:<15}  {:<10}'.format('kinDB_id', 'target_score', 'family', 'pk_length', 'npubs', '%muts', 'npk_pdbs', 'ndisease_assocs', 'nbioassays')
        return format_string

    else:
        if include_pseudogene_data:
            target_info_string = '{0[kinDB_id]:<25}  {0[target_score]: 12.1f}    {0[family]:<10}  {0[pk_domain_length]:<9}  {0[pubs]:<6}    {0[muts]:<8}  {0[npk_pdbs]:<8}  {0[disease]:<15}  {0[nbioassays]:<10}  {0[pseudogene]:<10}'.format(target_info_dict)
        else:
            target_info_string = '{0[kinDB_id]:<25}  {0[target_score]: 12.1f}    {0[family]:<10}  {0[pk_domain_length]:<9}  {0[pubs]:<6}    {0[muts]:<8}  {0[npk_pdbs]:<8}  {0[disease]:<15}  {0[nbioassays]:<10}'.format(target_info_dict)
        return target_info_string


def target_info_dict(pk_domain_node):
    '''
    Returns a list of target info dicts - one for each target domain within an entry node.
    '''
    entry_node = pk_domain_node.getparent().getparent()

    # info for this protein entry
    target_score_node = entry_node.find('target_score')
    pubs_score = int(float(target_score_node.get('publications')))
    cbioportal_mutations_score = target_score_node.get('cBioPortal_mutations')
    disease_assocs = len(entry_node.findall('uniprot/disease_association'))
    nbioassays = len(entry_node.findall('bioassays/bioassay'))

    # domain specific info
    kinDB_id = pk_domain_node.get('kinDB_id')
    pk_domain_id = int(pk_domain_node.get('id'))
    pk_domain_score_node = entry_node.find('target_score/pk_domain[@kinDB_id="%(kinDB_id)s"]' % vars())
    target_score = float( pk_domain_score_node.get('target_score') )
    family = pk_domain_node.getparent().get('family')
    pk_domain_length = pk_domain_node.get('length')
    npk_pdbs = pk_domain_score_node.get('npk_pdbs')
    pseudogene_score = float(pk_domain_score_node.get('pseudogene'))
    if pseudogene_score == 0.:
        pseudogene = ''
    else:
        pseudogene = 'pseudo'

    # Add info the target_info dict
    target_info_dict = {}
    target_info_dict['kinDB_id'] = kinDB_id
    target_info_dict['pk_domain_id'] = pk_domain_id
    target_info_dict['target_score'] = target_score
    target_info_dict['family'] = family
    target_info_dict['pk_domain_length'] = pk_domain_length
    target_info_dict['pubs'] = pubs_score
    target_info_dict['muts'] = cbioportal_mutations_score
    target_info_dict['disease'] = disease_assocs
    target_info_dict['nbioassays'] = nbioassays
    target_info_dict['npk_pdbs'] = npk_pdbs
    target_info_dict['pseudogene'] = pseudogene

    return target_info_dict


def diff_DB_comparison_strings(old_comparison_string, new_comparison_string):

    # Create a tempdir
    tempdirpath = tempfile.mkdtemp()

    try:
        # Write the new and old DBs to the tempdir, with the versionIDs and datestamps removed
        tempoldpath = os.path.join(tempdirpath, 'old.xml')
        tempnewpath = os.path.join(tempdirpath, 'new.xml')
        with open(tempoldpath, 'w') as tempoldfile:
            tempoldfile.write(old_comparison_string)
        with open(tempnewpath, 'w') as tempnewfile:
            tempnewfile.write(new_comparison_string)

        #stdout_pipe = subprocess.Popen(['diff', '--ignore-matching-lines=datestamp', tempoldDBpath, tempnewDBpath], stdout=subprocess.PIPE)
        stdout_pipe = subprocess.Popen(['diff', tempoldpath, tempnewpath], stdout=subprocess.PIPE)
        diff_output = stdout_pipe.communicate()[0]
        #print '\n'.join(new_comparison_string.split('\n')[0:20])
        #print '\n'.join(old_comparison_string.split('\n')[0:200])

    finally:
        shutil.rmtree(tempdirpath)

    return diff_output

def writeDB(DB_root, DB_filepath):
    with open(DB_filepath, 'w') as DB_file:
        from lxml import etree
        DB_file.write( etree.tostring(DB_root, pretty_print=True) )

