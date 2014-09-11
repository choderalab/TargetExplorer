def test_extract_sifts_seq():
    from targetexplorer.PDB import extract_sifts_seq
    sifts_filepath = 'tests/resources/4L00.xml.gz'
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
