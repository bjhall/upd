import logging
import gzip
import re

from codecs import (open, getreader)
from pprint import pprint as pp

LOG = logging.getLogger(__name__)


def open_file(filename):
    """Open a file and return a iterable with lines"""
    if filename.endswith('.gz'):
        LOG.info(f"{filename} is zipped")
        handle = getreader('utf-8')(gzip.open(filename), errors='replace')
    else:
        handle = open(filename, mode='r', encoding='utf-8', errors='replace')

    return handle


class Variant(object):
    """Implements a Variant class for VCF variants
    
    gt_types: 0=HOM_REF, 1=HET, 3=HOM_ALT, 2=other
    """
    def __init__(self, variant_line):
        super(Variant, self).__init__()
        self.variant_line = variant_line
        self.CHROM = None
        self.POS = None
        self.INFO = {}
        self.ALT = None
        self.is_snp = True
        self.gt_quals = []
        self.gt_types = []
        self._initialize()
    
    def _initialize(self):
        splitted_line = self.variant_line.split('\t')
        self.CHROM = splitted_line[0]
        self.POS = int(splitted_line[1])
        self.ALT = splitted_line[4].split(',')
        if len(splitted_line[3]) != len(splitted_line[4]):
            self.is_snp = False
        self.INFO = self._build_info(splitted_line[7])
        
        if not len(splitted_line) > 8:
            return
        
        self.gt_quals, self.gt_types = self._build_gt(splitted_line)
    
    def _build_gt(self, var_info):
        """Build the genotype information
        
        Collapse the FORMAT information from the 8th column with each individual genotype 
        information. Then retrieve the genotype call and genotype quality
        
        Args:
             var_info (list): A splited variant line
        
        Returns:
            gt_quals, gt_types (list)
        
        """
        gt_map = {'0/0':0, '0/1':1,'1/1':3}
        gt_types = []
        gt_quals = []
        
        form = var_info[8].split(':')
        for ind_info in var_info[9:]:
            gt_info = dict(zip(form, ind_info.split(':')))
            gq = 0
            try:
                gq = int(gt_info.get('GQ',0))
            except Exception as err:
                pass
            gt_quals.append(gq)
            genotype = gt_info.get('GT','./.')
            if not genotype in gt_map:
                gt_types.append(2)
                continue
            gt_types.append(gt_map[genotype])
        
        return gt_quals, gt_types

    def _build_info(self, info):
        """Build a info dictionary from a info str and set self.INFO
        
        Args:
            info (str): Raw vcf info string
        
        """
        info_dict = {}
        if info == '.':
            return info_dict
        for value in info.split(';'):
            vals = value.split('=')
            if not len(vals) == 2:
                info_dict[vals[0]] = True
                continue
            info_dict[vals[0]] = vals[1]
        return info_dict
    
    def __str__(self):
        return self.variant_line
    
    def __repr__(self):
        return f"{self.CHROM}:{self.POS}:{self.gt_types}"

class HREC(object):
    """VCF header record"""
    def __init__(self, hid, number=None, htype=None, description=None):
        super(HREC, self).__init__()
        self.id = hid
        self.number = number
        self.type = htype
        self.description = description
    
    def info(self):
        """Return the header info in a dictionary"""
        return {
            'ID': self.id, 
            'Number': self.number, 
            'Type': self.type, 
            'Description': self.description
        }
        

class Vcf(object):
    """Implements a simple vcf parser that mimics parts of cyvcf2.VCF"""
    def __init__(self, variant_file):
        super(Vcf, self).__init__()
        self.variant_file = iter(variant_file)
        self.raw_header = []
        self.samples = []
        self._current_variant = None
        self._header_keys = set()
        self._initialize()
    
    def _initialize(self):
        self._parse_header()
    
    def _parse_header(self):
        """docstring for _parse_header"""
        line = '#'
        while line.startswith('#'):
            line = next(self.variant_file)
            line = line.rstrip()
            if line.startswith('#'):
                self.raw_header.append(line)
                if not line.startswith('##'):
                    splitted_line = line.split('\t')
                    if not len(splitted_line) > 9:
                        raise SyntaxError("No individuals in VCF")
                    self.samples = splitted_line[9:]
        self._current_variant = line

    def contains(self, key):
        """Check if the header contains key"""
        if not self._header_keys:
            for rec in self.header_iter():
                self._header_keys.add(rec.id)
        
        return key in self._header_keys
            

    def header_iter(self):
        """Iterates over the header lines
        
        Creates header records (HREC) for each of the INFO headers
        
        Yields:
            header_record (HREC)
        """
        info_pattern = re.compile(r'''\#\#INFO=<
                    ID=(?P<id>[^,]+),
                    Number=(?P<number>-?\d+|\.|[AGR]),
                    Type=(?P<type>Integer|Float|Flag|Character|String),
                    Description="(?P<desc>[^"]*)"
                    >''', re.VERBOSE
                )
        
        for header in self.raw_header:
            if header.startswith('##INFO'):
                match = info_pattern.match(header)
                if not match:
                    raise SyntaxError(f"One of the INFO lines is malformed:{header}")
                
                header_record = HREC(match.group('id'), match.group('number'), 
                                        match.group('type'), match.group('desc'))
                yield header_record
    
    def __next__(self):
        current_variant = Variant(self._current_variant)
        self._current_variant = next(self.variant_file).rstrip()
        return current_variant
    
    def __iter__(self):
        return self
    
    def __repr__(self):
        return f"{self.__class__.__name__} ({self.samples})"

def check_samples(sids, proband, mother, father):
    """Check if proband, mother and father exists in vcf
    
    Args:
        sids (list): List of sample ids in VCF
        proband (str): ID of proband in VCF
        mother (str): ID of mother in VCF
        father (str): ID of father in VCF
    
    Returns:
        bool: if all samples exists in VCF
    """
    if not all(elem in sids for elem in [proband, mother, father]):
        return False
    
    return True

def get_vcf(vcf_path, proband, mother, father):
    """Check and open a VCF
    
    Args:
        vcf_path (str)
        proband (str): ID of proband in VCF
        mother (str): ID of mother in VCF
        father (str): ID of father in VCF
    
    Returns:
        vcf_reader (Vcf)
        
    """
    vcf_handle = open_file(vcf_path)
    vcf_reader = Vcf(vcf_handle)
    
    if not check_samples(vcf_reader.samples, proband, mother, father):
        raise SyntaxError("At least one of the given sample IDs do not exist in the VCF header")

    return vcf_reader

def get_header_desc(reader, header_id):
    """Get description field of an header field ID
    
    Args:
        reader (cyvcf2.VCF)
        header_id (str)
    
    Returns:
        str: Information from a vcf header description
    """
    for rec in reader.header_iter():
        d = rec.info()
        if d.get('ID') == header_id:
            return d.get('Description')

    return None

def parse_CSQ_header(reader):
    """Parse the order of VEP fields from the vcf header
    
    Args:
        reader (Vcf)

    Returns:
        csq_format (list(str)): A list with the VEP header
    """
    csq_str = get_header_desc(reader, "CSQ")
    
    if not csq_str:
        raise ValueError("CSQ header field missing. The VCF need to be annotated with VEP")
        
    _, csq_format_str = csq_str.split('Format: ')
    csq_format_str = csq_format_str.rstrip('"')
    csq_format = csq_format_str.split('|')
    return csq_format

def get_pop_AF(variant, vep_fields, af_tag):
    """Extract population frequency from VEP annotations.
    
    Args:
        variant (Variant)
        vep_fields (list): Description of VEP annotation
        af_tag (str): Name of AF field to parse
    
    Returns:
        freq (float): The annotated frequency, returns 0 if no data
    """
    freq = 0
    if vep_fields:
        vep_data = variant.INFO['CSQ']
        first_vep_str = vep_data.split(',')[0]
        data = first_vep_str.split('|')

        for i in range(len(data)):
            if vep_fields[i] == af_tag:
                freq = data[i]
    else:
        freq = variant.INFO.get(af_tag)
            
    return float(freq or 0)
