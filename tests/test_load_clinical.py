from nose.tools import eq_
import pytcga

def test_load_clinical():
    luad = pytcga.load_clinical_data('LUAD')

    eq_(len(luad), 522)
