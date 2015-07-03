from targetexplorer.flaskapp import models
from targetexplorer.tests.utils import projecttest_context
from targetexplorer.PDB import GatherPDB


def test_gather_pdb():
    with projecttest_context(set_up_project_stage='uniprot'):
        GatherPDB()
        first_pdb_row = models.PDB.query.first()
        assert isinstance(first_pdb_row, models.PDB)
        all_pdb_rows = models.PDB.query.all()
        assert '1OPL' in [pdb.pdbid for pdb in all_pdb_rows]
