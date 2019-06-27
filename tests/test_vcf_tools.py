import pytest

from upd.vcf_tools import (check_samples, get_vcf, Vcf)

def test_check_samples():
    ## GIVEN a list three samples
    sids = ['TEST_PROBAND', 'TEST_MOTHER', 'TEST_FATHER']
    proband = 'TEST_PROBAND'
    mother = 'TEST_MOTHER'
    father = 'TEST_FATHER'
    
    ## WHEN checking if the samples are in sids
    res = check_samples(sids, proband, mother, father)
    
    ## THEN assert that res is true
    assert res == True

def test_check_samples_missing_individual():
    ## GIVEN a list of two samples
    sids = ['TEST_PROBAND', 'TEST_MOTHER']
    proband = 'TEST_PROBAND'
    mother = 'TEST_MOTHER'
    father = 'TEST_FATHER'
    
    ## WHEN checking if the samples are in sids
    res = check_samples(sids, proband, mother, father)
    
    ## THEN assert that res is false
    assert res == False

def test_get_vcf(vcf_path):
    ## GIVEN the path to a VCF and the individuals that should be there
    proband = 'TEST_PROBAND'
    mother = 'TEST_MOTHER'
    father = 'TEST_FATHER'
    
    ## WHEN getting the vcf
    vcf_obj = get_vcf(vcf_path, proband, mother, father)
    
    ## THEN assert that vcf_obj has the correct type
    assert isinstance(vcf_obj, Vcf)

def test_get_nonexisting_vcf():
    ## GIVEN the path to a non existing VCF and the individuals that should be there
    proband = 'TEST_PROBAND'
    mother = 'TEST_MOTHER'
    father = 'TEST_FATHER'
    vcf_path = 'hej.vcf'
    
    ## WHEN getting the vcf
    with pytest.raises(OSError):
        ## THEN assert a non OSError is raised
        vcf_obj = get_vcf(vcf_path, proband, mother, father)

def test_get_wrong_formated_vcf(ped_path):
    ## GIVEN the path to a non existing VCF and the individuals that should be there
    proband = 'TEST_PROBAND'
    mother = 'TEST_MOTHER'
    father = 'TEST_FATHER'
    
    ## WHEN getting the vcf
    with pytest.raises(SyntaxError):
    ## THEN assert a non OSError is raised
        vcf_obj = get_vcf(ped_path, proband, mother, father)

    
def test_get_vcf_wrong_individual(vcf_path):
    ## GIVEN the path to a VCF and the individuals that should be there
    proband = 'TEST_PROBAND'
    mother = 'TEST_MOTHER'
    father = 'JOHN DOE'
    
    ## WHEN getting the vcf
    with pytest.raises(SyntaxError):
        ## THEN assert a SyntaxError is raised
        vcf_obj = get_vcf(vcf_path, proband, mother, father)
    