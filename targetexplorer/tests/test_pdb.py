import os
from targetexplorer.flaskapp import models
from targetexplorer.tests.utils import projecttest_context
from targetexplorer.pdb import GatherPDB, extract_pdb_data, extract_sifts_seq
from targetexplorer.utils import get_installed_resource_filepath
from nose.plugins.attrib import attr


# TODO will probably want to refactor extract_pdb_data at some point
# def test_extract_pdb_data():
#     pdb_dict = {
#         'pdb_row_id': 1,
#         'pdbid': '1OPL',
#         'ac': u'P00519',
#         'entry_name': u'ABL1_HUMAN',
#         'seq': u'MLEICLKLVGCKSKKGLSSSSSCYLEEALQRPVASDFEPQGLSEAARWNSKENLLAGPSENDPNLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPSNYITPVNSLEKHSWYHGPVSRNAAEYLLSSGINGSFLVRESESSPGQRSISLRYEGRVYHYRINTASDGKLYVSSESRFNTLAELVHHHSTVADGLITTLHYPAPKRNKPTVYGVSPNYDKWEMERTDITMKHKLGGGQYGEVYEGVWKKYSLTVAVKTLKEDTMEVEEFLKEAAVMKEIKHPNLVQLLGVCTREPPFYIITEFMTYGNLLDYLRECNRQEVNAVVLLYMATQISSAMEYLEKKNFIHRDLAARNCLVGENHLVKVADFGLSRLMTGDTYTAHAGAKFPIKWTAPESLAYNKFSIKSDVWAFGVLLWEIATYGMSPYPGIDLSQVYELLEKDYRMERPEGCPEKVYELMRACWQWNPSDRPSFAEIHQAFETMFQESSISDEVEKELGKQGVRGAVSTLLQAPELPTKTRTSRRAAEHRDTTDVPEMPHSKGQGESDPLDHEPAVSPLLPRKERGPPEGGLNEDERLLPKDKKTNLFSALIKKKKKTAPTPPKRSSSFREMDGQPERRGAGEEEGRDISNGALAFTPLDTADPAKSPKPSNGAGVPNGALRESGGSGFRSPHLWKKSSTLTSSRLATGEEEGGGSSSKRFLRSCSASCVPHGAKDTEWRSVTLPRDLQSTGRQFDSSTFGGHKSEKPALPRKRAGENRSDQVTRGTVTPPPRLVKKNEEAADEVFKDIMESSPGSSPPNLTPKPLRRQVTVAPASGLPHKEEAGKGSALGTPAAAEPVTPTSKAGSGAPGGTSKGPAEESRVRRHKHSSESPGRDKGKLSRLKPAPPPPPAASAGKAGGKPSQSPSQEAAGEAVLGAKTKATSLVDAVNSDAAKPSQPGEGLKKPVLPATPKPQSAKPSGTPISPAPVPSTLPSASSALAGDQPSSTAFIPLISTRVSLRKTRQPPERIASGAITKGVVLDSTEALCLAISRNSEQMASHSAVLEAGKNLYTFCVSYVDSIQQMRNKFAFREAINKLENNLRELQICPATAGSGPAATQDFSKLLSSVKEISDIVQR',
#         'chain_data': [{'chain_id': u'A', 'chain_row_id': 1}, {'chain_id': u'B', 'chain_row_id': 2}],
#         'structure_dirs': None
#     }
#     extract_pdb_data(pdb_dict)


@attr('unit')
def test_gather_pdb():
    with projecttest_context(set_up_project_stage='uniprot'):
        pdb_and_sifts_structure_files_dir = get_installed_resource_filepath(
            os.path.join('resources', 'structures')
        )
        GatherPDB(structure_dirs=pdb_and_sifts_structure_files_dir)
        first_pdb_row = models.PDB.query.first()
        assert isinstance(first_pdb_row, models.PDB)
        all_pdb_rows = models.PDB.query.all()
        assert '1OPL' in [pdb.pdbid for pdb in all_pdb_rows]


@attr('network')
def test_gather_pdb_using_network():
    with projecttest_context(set_up_project_stage='uniprot'):
        GatherPDB()
        first_pdb_row = models.PDB.query.first()
        assert isinstance(first_pdb_row, models.PDB)
        all_pdb_rows = models.PDB.query.all()
        assert '1OPL' in [pdb.pdbid for pdb in all_pdb_rows]


@attr('unit')
def test_extract_sifts_seq():
    sifts_filepath = get_installed_resource_filepath(os.path.join(
        'resources', '4L00.xml.gz'
    ))

    seq = 'MQYLNIKEDCNAMAFCAKMRSSKKTEVNLEAPEPGVEVIFYLSDREPLRLGSGEYTAEEL\
CIRAAQACRISPLCHNLFALYDENTKLWYAPNRTITVDDKMSLRLHYRMRFYFTNWHGTN\
DNEQSVWRHSPKKQKNGYEKKKIPDATPLLDASSLEYLFAQGQYDLVKCLAPIRDPKTEQ\
DGHDIENECLGMAVLAISHYAMMKKMQLPELPKDISYKRYIPETLNKSIRQRNLLTRMRI\
NNVFKDFLKEFNNKTICDSSVSTHDLKVKYLATLETLTKHYGAEIFETSMLLISSENEMN\
WFHSNDGGNVLYYEVMVTGNLGIQWRHKPNVVSVEKEKNKLKRKKLENKHKKDEEKNKIR\
EEWNNFSYFPEITHIVIKESVVSINKQDNKKMELKLSSHEEALSFVSLVDGYFRLTADAH\
HYLCTDVAPPLIVHNIQNGCHGPICTEYAINKLRQEGSEEGMYVLRWSCTDFDNILMTVT\
CFEKSEQVQGAQKQFKNFQIEVQKGRYSLHGSDRSFPSLGDLMSHLKKQILRTDNISFML\
KRCCQPKPREISNLLVATKKAQEWQPVYPMSQLSFDRILKKDLVQGEHLGRGTRTHIYSG\
TLMDYKDDEGTSEEKKIKVILKVLDPSHRDISLAFFEAASMMRQVSHKHIVYLYGVCVRD\
VENIMVEEFVEGGPLDLFMHRKSDVLTTPWKFKVAKQLASALSYLEDKDLVHGNVCTKNL\
LLAREGIDSECGPFIKLSDPGIPITVLSRQECIERIPWIAPECVEDSKNLSVAADKWSFG\
TTLWEICYNGEIPLKDKTLIEKERFYESRCRPVTPSCKELADLMTRCMNYDPNQRPFFRA\
IMRDINKLEEQNPDIVSEKKPATEVDPTHFEKRFLKRIRDLGEGHFGKVELCRYDPEGDN\
TGEQVAVKSLKPESGGNHIADLKKEIEILRNLYHENIVKYKGICTEDGGNGIKLIMEFLP\
SGSLKEYLPKNKNKINLKQQLKYAVQICKGMDYLGSRQYVHRDLAARNVLVESEHQVKIG\
DFGLTKAIETDKEYYTVKDDRDSPVFWYAPECLMQSKFYIASDVWSFGVTLHELLTYCDS\
DSSPMALFLKMIGPTHGQMTVTRLVNTLKEGKRLPCPPNCPDEVYQLMRKCWEFQPSNRT\
SFQNLIEGFEALLK'

    pdb_chain_obj = extract_sifts_seq(sifts_filepath, 'P23458', 'JAK1_HUMAN', '4L00', 'A', seq)

    assert pdb_chain_obj['experimental_seq_aln_conflicts'] == '--------\
-----------------------------------------------------------------------\
-----------------------------------------------------------------------\
-----------------------------------------------------------------------\
-----------------------------------------------------------------------\
-----------------------------------------------------------------------\
-----------------------------------------------------------------------\
-----------------------------------------------------------------------\
---------------------------------------------------gstsAQEWQPVYPMSQLSFD\
RILKKDLVQGEHLGRGTRTHIYSGTLMDYKDDEGTSEEKKIKVILKVLDPSHRDISLAFFEAASMMRQVSH\
KHIVYLYGVCVRDVENIMVEEFVEGGPLDLFMHRKSDVLTTPWKFKVAKQLASALSYLEDKDLVHGNVCTK\
NLLLAREGIDSECGPFIKLSDPGIPITVLSRQECIERIPWIAPECVEDSKNLSVAADKWSFGTTLWEICYN\
GEIPLKDKTLIEKERFYESRCRPVTPSCKELADLMTRCMNYDPNQRPFFRAIMRDINKLEEQNPDIVSEKK\
-----------------------------------------------------------------------\
-----------------------------------------------------------------------\
-----------------------------------------------------------------------\
-----------------------------------------------------------------------\
----------'

    assert pdb_chain_obj['observed_seq_aln'] == '-----------------------\
-----------------------------------------------------------------------\
-----------------------------------------------------------------------\
-----------------------------------------------------------------------\
-----------------------------------------------------------------------\
-----------------------------------------------------------------------\
-----------------------------------------------------------------------\
-----------------------------------------------------------------------\
-------------------------------------------WQPVYPMSQLSFDRILKKDLVQGEHLGR\
GTRTHIYSGTLMD----------KKIKVILKVLDPSHRDISLAFFEAASMMRQVSHKHIVYLYGVCVRDVE\
NIMVEEFVEGGPLDLFMHRKSDVLTTPWKFKVAKQLASALSYLEDKDLVHGNVCTKNLLLAREGIDSECGP\
FIKLSDPGIPITVLSRQECIERIPWIAPECVEDSKNLSVAADKWSFGTTLWEICYNGEIPLKDKTLIEKER\
FYESRCRPVTPSCKELADLMTRCMNYDPNQRPFFRAIMRDINKLEE-------------------------\
-----------------------------------------------------------------------\
-----------------------------------------------------------------------\
-----------------------------------------------------------------------\
------------------------------------------------------------------'
