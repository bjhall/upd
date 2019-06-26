from cyvcf2 import VCF

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
        vcf_reader (cyvcf2.VCF)
        
    """
    vcf_reader = VCF(vcf_path)
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
        reader (cyvcf2.VCF)

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
        variant (cyvcf2.Variant)
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
