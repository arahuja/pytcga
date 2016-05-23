from nose.tools import eq_
import pytcga

def test_load_mutations():

    luad_mutations = pytcga.load_mutation_data('LUAD')
    eq_(len(luad_mutations), 255972)

    eq_(luad_mutations.TCGA_ID.nunique(), 542)