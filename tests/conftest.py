import pytest

@pytest.fixture()
def vcf_path():
    return 'tests/fixtures/test.vcf.gz'

@pytest.fixture()
def regions_output():
    line = ("15\t22958118\t101910531\tORIGIN=PATERNAL;TYPE=HETERODISOMY;LOW_SIZE=78952412;"
            "INF_SITES=138;SNPS=2032;HET_HOM=889/1006;OPP_SITES=0;START_LOW=20170126;END_HIGH="
            "102516586;HIGH_SIZE=82346460")
    return line


@pytest.fixture()
def ped_path():
    return 'tests/fixtures/test.fam'
