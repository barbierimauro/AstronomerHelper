#!/bin/python3

class SesameQuery:
    """
    A class to query astronomical objects from the SESAME service (SIMBAD, VizieR, NED).
    
    Example usage:
    from sesame_query import SesameQuery

    query = SesameQuery()
    result = query.search("HD42499A")
    print(result)
    """

    def __init__(self):
        self.SESAME_OUTPUT = "-ox"
        self.SESAME_SERVER = "http://vizier.u-strasbg.fr/cgi-bin/nph-sesame"

    @staticmethod
    def is_successful(xml_string):
        target_count = xml_string.count("<Target ")
        not_found_count = xml_string.count("Nothing found")
        return target_count - not_found_count

    @staticmethod
    def parse_values(xml_string):
        root = ElementTree.fromstring(xml_string)
        ra = root.findtext(".//jradeg", "")
        dec = root.findtext(".//jdedeg", "")
        nome = re.sub(r'[{}_]', '', root.findtext(".//oname", "")).strip()

        ra_e = dec_e = oid = otype = pmra = pmra_e = pmde = pmde_e = sptype = nrefs = plx = plx_e = rv = rv_e = ""

        if '<Target option="S">' in xml_string:
            ra_e = root.findtext(".//errRAmas", "")
            dec_e = root.findtext(".//errDEmas", "")
            oid = root.findtext(".//oid", "")
            otype = root.findtext(".//otype", "")
            pmra = root.findtext(".//pmRA", "")
            pmra_e = root.findtext(".//epmRA", "")
            pmde = root.findtext(".//pmDE", "")
            pmde_e = root.findtext(".//epmDE", "")
            sptype = root.findtext(".//spType", "")
            nrefs = root.findtext(".//nrefs", "")
            plx = root.findtext(".//plx/v", "")
            plx_e = root.findtext(".//plx/e", "")
            rv = root.findtext(".//Vel/v", "")
            rv_e = root.findtext(".//Vel/e", "")

        return ra, dec, ra_e, dec_e, nome, oid, otype, pmra, pmra_e, pmde, pmde_e, sptype, nrefs, plx, plx_e, rv, rv_e

    def query_sesame(self, q):
        for resolver in ["S", "V", "N"]:
            url = f"{self.SESAME_SERVER}/{self.SESAME_OUTPUT}/{resolver}?{q}"
            response = requests.get(url).text
            success = self.is_successful(response)

            if success == 1:
                return self.parse_values(response), resolver

        print("Error: All resolvers failed.")
        sys.exit(1)

    def search(self, object_name):
        q = quote_plus(object_name.strip()).upper()
        parsed_values, resolver = self.query_sesame(q)

        labels = ["RESOLV", "RA", "DEC", "RA_E", "DEC_E", "NOME", "OID", "OTYPE", "PMRA", "PMRA_E", "PMDE", "PMDE_E", "SPTYPE", "NREFS", "PLX", "PLX_E", "RV", "RV_E"]
