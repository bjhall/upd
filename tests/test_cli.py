from click.testing import CliRunner

from upd.cli import cli
from upd.__version__ import __version__

def test_version():
    runner = CliRunner()
    result = runner.invoke(cli, ['--version'])
    
    assert result.exit_code == 0
    assert result.output == __version__+'\n'

def test_upd(vcf_path):
    runner = CliRunner()
    result = runner.invoke(cli, [
        '--vcf', vcf_path, '--proband', 'TEST_PROBAND', '--mother', 'TEST_MOTHER', '--father', 
        'TEST_FATHER', '--vep', 'regions'
    ])
    
    assert result.exit_code == 0

