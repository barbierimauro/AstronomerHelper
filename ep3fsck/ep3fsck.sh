#!/bin/bash
# ep3fsck : eso phase 3 fits checker
#
author="Mauro Barbieri"
year=2022
email=maurobarbieri.science@gmail.com

#
# HISTORY
# -------
# MB 2022-05-30T18:00:00_v1.0  : first release
# MB 2022-06-07T12:00:00_v2.0  : improved LaTeX output
# MB 2022-06-10T09:00:00_v3.0  : new version, implemented checks on all keywords
# MB 2022-06-14T12:00:00_v4.0  : verbatim output in color, check kw if condition apply to primary
#                                header or extensions
# MB 2022-06-15T18:00:00_v5.0  : check on path of ancillary programs, defined command line,
#                                download SIMBAD information, download DSS cutout, added help,
#                                added creation of auto documentation,
#                                added download of standards documents
# MB 2022-06-20T12:00:00_v6.0  : fixed command line, increased control of input parameters
# MB 2022-06-21T12:00:00_v7.0  : fixed path basic commands, substitued realpath with readlink,
#                                added check on bash version, check on latex errors
# MB 2022-06-22T16:34:35_v7.1  : check version of fitsverify, create function for check presence of 
#                                dfits and fitsverify, improved help, added autoedit, removed useless
#                                command line parameters, removed usage of gnome-open, improved
#                                handling of report generation removing the test on the latex compilation,
#                                improved internal comments, move all definitions at beginning of code,
#                                added snapshot of images
# MB 2022-06-23T13:40:53_v7.2  : added table of keyword types, added output of kw of selected dataprod,
#                                checking the presence of extensions
# MB 2022-07-05T16:43:37_v7.3  : corrected the behaviour on kw with multiple values, added option for adding keyword table,
#                                printing only errors, added option for getting ASCII output

# corrected the function for testing the data kind
# added numerical checks on RA,DEC,EXPTIME,MJD-OBS
# added check on unsupported kw
# added check on datatype of every kw
# added table of notes, need to implement a function that show the notes in the latex output
#
# Nausicaa request: print phase3 SDPS keywords for ADP files (for generating templates headers)


#
# remember to copy the version from history !!!
#
ver=2022-07-27T13:35:02_v7.4
vern=${ver:21:23}
verd=${ver:0:19}


# DESCRIPTION
# -----------
# Check FITS files for compliance with ESO PHASE3 keywords (SDPS v8 March 2022)                       
#
#



###############################################################################
#                                                                             #
# COMMON SETTINGS                                                             #
#                                                                             #
#                                                                             #
#                                                                             #
#                                                                             #
###############################################################################

padding="........"
hdr=hdr.tmp
simb=simbad.out
collection="Report"
user=$USER
oaout=false
oafile=out.txt
ofile=out.tex
opdf="${ofile/.tex/.pdf}"
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}"; )" &> /dev/null && pwd 2> /dev/null; )";
SDPS="SDPS v8 March 2022"
simbad=false
allheaders=false
latex_=false
fitsverify_=false
fileinfo_=false
image_=false
printkwlist=false
proginfo=false
provflag=false
bibflag=false
fileinfo=false



###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################


# definition of dataprod, and keyword list and flags
data_types=("image" "image_mef" "mef_image" "image_fluxmap" "spectrum" "cube_ifs" "visibility" "srctbl" "m_catalog" "catalog_tile" "catalog")
data_types_short=("image" "imef" "mefi" "iflux" "spec" "cube" "visi" "srctbl" "mcat" "ctile" "catal")
#
basic_kwlist=(SIMPLE BITPIX NAXIS NAXIS1 NAXIS2 NAXIS3 EXTEND XTENSION EXTNAME)
kwlist=(PRODCATG ASSOC ASSON ASSOM ORIGIN TELESC INSTR FILTER OBJECT RA DEC EQUINOX RADESYS TIMESYS EXPTIME TEXPTIME MJD-OBS MJD-END PROG OBID NCOMBINE OBSTECH FLUXCAL PROCSOFT REFERENC PROV BUNIT GAIN DETRON EFFRON WEIGHT CRVAL CRPIX CTYPE CUNIT CD CSYER CRDER PHOTZP PHOTZPER PHOTSYS SPECSYS EXT_OBJ CONTNORM TOT_FLUX FLUXERR WAVELMIN WAVELMAX LAMRMS LAMNLIN SPEC_BIN SPEC_ERR SPEC_SYE RA_ERR DEC_ERR NELEM VOCLASS VOPUB TITLE APERTURE TELAPSE TMID SPEC_VAL SPEC_BW BNOISE MAPMODE FEBE CONTENT INSMODE BASE_MIN BASE_MAX NUM_CHAN VIS2ERR T3PHIERR STOKES HDUCLASS HDUCLAS HDUDOC HDUVERS SCIDATA ERRDATA QUALDATA CONFDATA BKGDATA BKGERR BKGCONF ABMAGLIM PIXNOISE MAGLIM ABMAGSAT PSF_FWHM ELLIPTIC SNR SPEC_RES SKY_RES SKY_RERR STREHL ARCFILE CHECKSUM DATASUM ORIGFILE P3ORIG NDIT NJITTER NOFFSETS NUSTEP FPRA FPDE SKYSQDEG M_EPOCH APMATCHD TXLNK TXRGF TXCTY NOESODAT TFIELDS TTYPE TFORM TCOMM TUNIT TUTYP TUCD TDMIN TDMAX TNULL EXTNAME TZERO TSCAL EXTVER EXTLEVEL)
##kwlist=("PRODCATG" "ASSOC   " "ASSON   " "ASSOM   " "ORIGIN  " "TELESC  " "INSTR   " "FILTER  " "OBJECT  " "RA      " "DEC     " "EQUINOX " "RADESYS " "TIMESYS " "EXPTIME " "TEXPTIME" "MJD-OBS " "MJD-END " "PROG    " "OBID    " "NCOMBINE" "OBSTECH " "FLUXCAL " "PROCSOFT" "REFERENC" "PROV    " "BUNIT   " "GAIN    " "DETRON  " "EFFRON  " "WEIGHT  " "CRVAL   " "CRPIX   " "CTYPE   " "CUNIT   " "CD      " "CSYER   " "CRDER   " "PHOTZP  " "PHOTZPER" "PHOTSYS " "SPECSYS " "EXT_OBJ " "CONTNORM" "TOT_FLUX" "FLUXERR " "WAVELMIN" "WAVELMAX" "LAMRMS  " "LAMNLIN " "SPEC_BIN" "SPEC_ERR" "SPEC_SYE" "RA_ERR  " "DEC_ERR " "NELEM   " "VOCLASS " "VOPUB   " "TITLE   " "APERTURE" "TELAPSE " "TMID    " "SPEC_VAL" "SPEC_BW " "BNOISE  " "MAPMODE " "FEBE    " "CONTENT " "INSMODE " "BASE_MIN" "BASE_MAX" "NUM_CHAN" "VIS2ERR " "T3PHIERR" "STOKES  " "HDUCLASS" "HDUCLAS1" "HDUDOC  " "HDUVERS " "SCIDATA " "ERRDATA " "QUALDATA" "CONFDATA" "BKGDATA " "BKGERR  " "BKGCONF " "ABMAGLIM" "PIXNOISE" "MAGLIM  " "ABMAGSAT" "PSF_FWHM" "ELLIPTIC" "SNR     " "SPEC_RES" "SKY_RES " "SKY_RERR" "STREHL  " "ARCFILE " "CHECKSUM" "DATASUM " "ORIGFILE" "P3ORIG  " "NDIT    " "NJITTER " "NOFFSETS" "NUSTEP  " "FPRA    " "FPDE    " "SKYSQDEG" "M_EPOCH " "APMATCHD" "TXLNK   " "TXRGF   " "TXCTY   " "NOESODAT" "TFIELDS " "TTYPE   " "TFORM   " "TCOMM   " "TUNIT   " "TUTYP   " "TUCD    " "TDMIN   " "TDMAX   " "TNULL   " "EXTNAME " "TZERO   " "TSCAL   " "EXTVER  " "EXTLEVEL")
kflag_legend=("MM" "ME" "MA" "RC" "NO" "OP" "na" "un" "RE" "SU")
kflag_legend_num=(1 2 3 4 5 6 7 8 9 0)
kflag_legend_expl=("MANDATORY" "MANDATORY ESO" "MANDATORY WHEN APPLICABLE, SEE SDPS" "RECOMMENDED" "KW NOT ALLOWED" "OPTIONAL" "KW NOT REQUIRED" "KW NOT REQUIRED" "KW RESERVED FOR PHASE3 PROCESS" "KW POLARIZATION")

kwtype=(S S S S S S S S S F F F S S F F F F S I I S S S S S S F F F F F F S S F F F F F F S L L L F F F F F F F F F F I S S S F F F F F F S S S S F F I F F S S S S S S S S S S S S F F F F F F F F F F F S S S S S I I I I F F F L L S S S L I S S S S S S F F I S U U U U)
# S = string
# F = float
# I = integer
# L = logical
# U = undefined



#1 = MM = MANDATORY
#2 = ME = MANDATORY ESO
#3 = MA = MANDATORY WHEN APPLICABLE
#4 = RC = RECOMMENDED
#5 = NO = NOT ALLOWED
#6 = OP = OPTIONAL
#7 = na = NOT APPLICABLE
#8 = un = UNKNOWN
#9 = RE = RESERVED
#0 = SU = ASK SUPPORT

#kw flag definitions
#kwflag_image =(11 13 13 13 11 11 11 11 11 11 11 13 11 13 11 11 11 11 12 12 12 11 11 11 11 12 11 16 16 16 16 11 11 11 11 11 14 14 11 14 11 17 17 17 17 17 17 17 17 17 17 17 17 15 15 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 10 17 17 17 17 17 17 17 17 17 17 17 11 17 17 11 11 11 17 17 17 17 13 19 31 31 19 19 13 13 13 13 17 17 17 17 17 17 17 17 13 17 17 17 17 17 17 17 17 17 17 23 18 18 18 18)
#kwflag_imef  =(11 13 13 13 11 11 11 11 11 11 11 13 11 13 11 11 11 11 12 12 12 11 11 11 11 12 21 16 16 16 16 21 21 21 21 21 24 24 21 24 21 17 17 17 17 17 17 17 17 17 17 17 17 15 15 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 10 21 21 21 21 23 23 23 23 23 23 23 21 17 17 21 21 21 17 17 17 17 13 19 31 31 19 19 13 13 13 13 17 17 17 17 17 17 17 17 13 17 17 17 17 17 17 17 17 17 17 23 18 18 18 18)
#kwflag_mefi  =(11 13 13 13 11 11 11 11 11 11 11 13 11 13 11 11 11 11 12 12 12 11 11 11 11 12 21 16 16 16 16 21 21 21 21 21 24 24 21 24 21 17 17 17 17 17 17 17 17 17 17 17 17 15 15 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 10 17 17 17 17 17 17 17 17 17 17 17 21 17 17 21 21 21 17 17 17 17 13 19 31 31 19 19 13 13 13 13 17 17 17 17 17 17 17 17 13 17 17 17 17 17 17 17 17 17 17 23 18 18 18 18)
#kwflag_iflux =(11 13 13 13 11 11 11 11 11 11 11 13 11 13 11 11 11 11 12 12 12 11 11 11 11 12 11 17 17 17 17 11 11 11 11 11 14 14 17 17 17 17 17 17 17 11 11 11 17 17 17 17 17 15 15 17 17 17 17 17 17 17 17 17 11 11 11 17 17 17 17 17 17 17 10 17 17 17 17 17 17 17 17 17 17 17 15 17 17 11 11 11 17 17 11 13 17 19 31 31 19 19 17 17 17 17 17 17 17 17 17 17 17 17 13 17 17 17 17 17 17 17 17 17 17 23 18 18 18 18)
#kwflag_spec  =(11 13 13 13 11 11 11 17 31 31 31 13 11 13 11 11 11 11 12 12 12 11 11 11 11 12 15 17 17 17 17 17 17 17 17 15 17 17 17 17 17 11 13 11 11 11 11 11 16 16 11 14 14 14 14 21 21 21 21 21 21 21 21 21 17 17 17 17 17 17 17 17 17 17 10 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 11 11 17 17 13 19 31 31 19 19 17 17 17 17 17 17 17 17 17 17 17 17 13 21 21 21 26 21 21 21 23 23 17 23 18 18 18 18)
#kwflag_cube  =(11 13 13 13 11 11 11 17 11 11 11 13 11 13 11 11 11 11 12 12 12 11 11 11 11 12 21 17 17 17 17 21 21 21 21 21 24 24 27 17 17 11 17 17 17 17 11 11 17 17 15 17 17 15 15 27 27 27 27 27 27 27 27 27 17 17 17 17 17 17 17 17 17 17 10 21 21 21 21 23 23 23 17 17 17 17 11 11 17 17 17 16 17 11 11 13 13 19 31 31 19 19 17 17 17 17 17 17 17 17 17 17 17 17 13 27 27 27 27 27 27 27 27 27 17 27 18 18 18 18)
#kwflag_visi  =(11 13 13 13 11 11 11 17 11 11 11 13 11 13 11 11 11 11 12 12 12 11 17 11 11 12 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 11 17 17 17 17 11 11 17 17 17 11 11 17 17 17 17 17 17 17 17 17 17 17 17 17 17 11 11 11 11 11 11 11 10 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 11 17 17 17 19 31 31 19 19 17 17 17 17 17 17 17 17 17 17 17 17 13 21 21 21 26 23 17 26 26 26 26 23 18 18 18 18)
#kwflag_srctbl=(11 13 13 13 11 11 11 11 11 11 11 13 11 13 11 11 11 11 12 12 12 11 17 11 11 12 17 17 17 17 17 17 17 17 17 17 17 17 23 24 21 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 10 17 17 17 17 17 17 17 17 17 17 17 21 17 17 21 21 21 17 17 17 17 13 19 31 31 19 19 17 17 17 17 17 17 17 17 17 17 17 17 13 21 21 21 26 23 17 26 26 26 26 23 18 18 18 18)
#kwflag_mcat  =(11 17 17 17 11 11 11 13 11 17 17 17 17 13 17 17 11 11 12 12 17 13 17 11 11 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 11 17 17 17 17 17 11 11 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 10 17 17 17 17 17 17 17 17 17 17 17 13 17 13 17 17 17 17 17 17 17 17 19 31 31 19 19 17 17 17 17 13 13 11 13 13 13 13 13 13 21 21 21 21 21 17 21 26 26 23 23 15 15 15 15)
#kwglag_ctile =(11 13 13 13 11 11 11 13 11 11 11 13 11 13 17 17 11 11 12 12 17 13 17 11 11 12 17 17 17 17 17 17 17 17 17 17 17 17 17 17 11 17 17 17 17 17 11 11 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 10 17 17 17 17 17 17 17 17 17 17 17 13 17 13 17 17 17 17 17 17 17 17 19 31 31 19 19 17 17 17 17 11 11 11 13 13 13 17 17 13 21 21 21 21 21 17 21 26 26 23 23 15 15 15 15)
#kwflag_catal =(11 13 13 13 11 11 11 13 11 13 13 13 13 13 17 17 11 11 12 12 17 13 17 11 11 12 17 17 17 17 17 17 17 17 17 17 17 17 17 17 11 17 17 17 17 17 11 11 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 13 13 13 17 17 17 17 17 17 17 10 17 17 17 17 17 17 17 17 17 17 17 13 17 13 17 17 17 17 17 13 13 17 19 31 31 19 19 17 17 17 17 13 13 13 13 13 13 13 13 13 21 21 21 21 21 17 21 26 26 23 23 15 15 15 15)
#

# convert many rows in to 1 row
# useful for convert flag tables
# perl -00 -lpe 's/\n/ /g' spectrum2.tab
#
# convert many rows in to 1 row
# cat flag_kind.tab | tr '\n' ' 
# convert one line in many lines
#echo 'foo bar boo you too' | tr ' ' '\n'

#image
#imef
#mefi
#iflux
#spec
#cube
#visi
#srctbl
#mcat
#ctile
#catal


function shownotes() {
#This function does nothing except printing all the notes and the table of the correspondances
#between kw an notes
cat <<EONOTES
31 Mandatory for the cases where ancillary files are provided in a non-FITS format.
32 Mandatory for the cases where ancillary files are provided in association with the scientific data.
33 The flux map must always be associated with the RMS noise map, or the SNR map, or both.
34 The IFS 3D cube must always be associated with the white-light image. The white-light imagemust be delivered using the ASSONi/ASSOCi mechanism.
35 Applicable to photometric catalogues as well as APEX catalogues.
36 There can be cases where that keyword does not apply, for instance in the case of target catalogues of spectroscopic surveys, e.g. PESSTO where no appropriate RA/DEC value can be defined to characterise the catalogue as a whole.
37 If RADESYS='FK5', EQUINOX=2000.0 is mandatory. If RADESYS = 'ICRS', EQUINOX is tolerated and its value needs to be 2000.0.
38 Must be present if the system used is other than UTC.
39 Does not apply to catalogues for which no unique value of OBSTECH can be identified.
40 If a refereed publication is not available at the time of the data release going public, the value can be left to an empty string.
41 Mandatory depending on whether fluxes or magnitudes are provided in the source list.
42 EXT_OBJ is mandatory for externally submitted data products (e.g. spectroscopic public surveys). EXT_OBJ is not applicable to data processed in an unsupervised way, for which the keyword value cannot be assessed and thus the property is not known.
43 FLUXERR applies to SCIENCE.SPECTRUM with FLUXCAL='ABSOLUTE'. In case of SCIENCE.SPECTRUM with FLUXCAL = 'UNCALIBRATED', the FLUXERR keyword shall be set to -1. The special value -2 is reserved for the case when the flux error cannot be determined.
44 Applicable to spectroscopic and APEX catalogues. For photometric catalogues, the value is calculated by the Phase 3 system unless the combination (INSTRi, FILTERi) is not unique, in which case please contact Phase 3 operations support staff at https://support.eso.org/ to assess the correct values of the WAVELMIN/MAX keywords, to be added in the headers.
45 For APEX catalogues only.
46 SCIDATA is mandatory for the cases where ancillary extensions are provided in association with the scientific data. ERRDATA / QUALDATA / CONFDATA / BKGDATA / BKGERR / BKGCONF shall be used if the corresponding extension is provided.
47 For photometric catalogues. And if there is more than one filter, it is not applicable. Use MAGLIMi instead.
48 For photometric catalogues with more than one filter.
49 For VIRCAM and OmegaCAM only.
50 Applicable to the case when SKY_RES is expected to vary within the data collection due to the way it is estimated (see footnote 18).
51 For AO observations only.
52 NIR image data products qualify for the keyword if, and only if, all exposures and observations contributing to the given product share the same value for the respective parameter. If, for example, the product has been created from exposures taken with different detector integration time, the keyword DIT should not be defined in the FITS header.
53 Not mandatory in case of complex footprints.
54 Does not apply to spectroscopic catalogues for which no coverage pattern exists.
55 Applicable to multi-epoch catalogues formatted according to section 12.5.1 only.
56 For aperture-matched catalogues only.
57 In case data link is used (see section 12.3).
58 In case provenance per catalogue record is used (see section 5.2.3).
59 Applicable to products originating or containing data from a non-ESO facility.
60 Keyword may be absent for columns representing quantities having no units of measurement, otherwise it must be present.
61 For quantities having no units of measurement, the value shall be set to an empty string.
62 In case the UType is not defined in the IVOA document [9] like e.g. for CONTINUUM, the corresponding TUTYPi keyword shall be set to an empty string.
63 Mandatory for the TTYPE1 array only (start/stop spectral coordinates).
64 See section 5.18 and section 12.2.6 for applicability.
65 A specific value is requested in particular cases: 'PHASE3PROVENANCE', 'PHASE3CATALOG', 'PHASE3FILELIST'.
TABLE NOTES
image 
 2  3  4 11 13 24 91 96 102 103 104 105 114 125 
31 32 31 37 38 40 49 51  52  52  52  52  59  65
imagemef
 2  3  4 11 13 24 79 80 81 82 83 84 85 91 96 102 103 104 105 114 125 
31 32 31 37 38 40 46 46 46 46 46 46 46 49 51  52  52  52  52  59  65
mefimage == image
 2  3  4 11 13 24 91 96 102 103 104 105 114 125 
31 32 31 37 38 40 49 51  52  52  52  52  59  65
fluxmap 
 2  3  4 11 13 24 95 114 125 
31 33 31 37 38 40 50  59  65
spectra
 2  3  4 11 13 24 42 45 96 114 119 120 122 123 125
31 32 31 37 38 40 42 43 51  59  61  62  63  63  65
cube
 2  3  4 11 13 24 79 80 81 95 96 114 125
31 34 31 37 38 40 46 46 46 50 51  59  65
visibility
 2  3  4 11 13 24 114 119 125 
31 33 31 37 38 40  59  60  65
srctbl
2 3 4 11 13 21 24 38 96 114 119 125
31 32 31 37 38 39 40 41 51  59  60  65
mcatalog
 8 13 21 24 46 47 86 88 106 107 109 110 111 112 113 114 119 124 125
35 38 39 40 44 44 47 48  53  53  55  56  57  57  57  59  61  64  65
catalogtile
 2  3  4  8 13 21 24 46 47 86 88 109 110   111 114 119 124 125
31 32 31 35 38 39 40 44 44 47 48  55  56 57+58  59  61  64  65
catalog
 2  3  4  8 10  11   12 13 21 24 46 47 64 65 66 86 88 94 95 106 107 108 109 110 111 112 113 114 119 124 125
31 32 31 35 36 36+37 36 38 39 40 44 44 45 45 45 47 48 45 50  54  54  54  55  56 578  57  57  59  61  64  65
EONOTES
}

function devnotes() {
cat <<EODEVNOTES

1)
nel basic check, dopo aver contato il numero di estensioni
dire all'utente se nel caso di spettri e cubi, i dati sono contenuti in una estensione
oppure se sono nell'HDU primario


2)
se la kw progid e' definita si possono recuperare i seguenti dati

BIBLIO in csv
https://telbib.eso.org/export.php?programid=67.B-0108&total=9&rows=30000&fl=author,title,journal,volume,pages,year,bibcode,telescope,instrument,programid,datastatus,survey,sciver&export_format=csv

ABSTRACT in html
prendo il testo dopo <b>Abstract</b>
http://archive.eso.org/wdb/wdb/eso/abstract/query?progid=67.B-0108

FILES in html
i nomi dei file sono dopo "value="
http://archive.eso.org/wdb/wdb/eso/eso_archive_main/query?prog_id=67.B-0108&max_rows_returned=10000
EODEVNOTES
}

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################



function logo() {
echo "+---------------------------------+"
echo "| ep3fsck                         |"
echo "| ver:" $vern"                        |"
echo "|                                 |"
echo "| "$author "(c)" $year   "        |"
echo "| "$email                        "|"
echo "+---------------------------------+"
#echo "This is free software; see the source for copying conditions.  There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE."
echo
}

function is_power_of_two () {
    declare -i n=$1
    (( n > 0 && (n & (n - 1)) == 0 ))
}

function toBinary(){
    local n bit
    for (( n=$1 ; n>0 ; n >>= 1 )); do  bit="$(( n&1 ))$bit"; done
    printf "%s\n" "$bit"
    echo ${#bit}
}

function convertIntvalToBase () # (Val Base)
{
   val=$1
   base=$2
   result=""
   while [ $val -ne 0 ] ; do
        result=$(( $val % $base ))$result #residual is next digit
        val=$(( $val / $base ))
   done
   echo $result
}

function varkind() {
#INPUT  variable_name
#OUTPUT vk : variable kind
#I = integer, F = float, E = exponential, C = complex, H = hexadecimal
#L = logical, NAN, NULL, INF, S = string, VARNOTDEF, U = undefined type 
#
#USAGE:
#x=12.34e7          #define the variable
#varkind x          #call the function
#echo $vk           #visualize the data kind
#
#caveat: if the variable is read from some file/process it coul contains
#leading and trailing spaces, that could interfere in the recognition of
#the correct data kind. Those spaces could be eliminated with awk
#in the following way creating a new variable named y
#y=$(echo $x | awk '{$1=$1};1')
#
local -n myvar=$1
vk="0x0"
sbugga=true
[ sbugga = 'true' ] && echo
[ sbugga = 'true' ] && echo

lenmyvar=${#myvar}


if [ "$myvar" -eq "$myvar" ] 2> /dev/null
then
    rc0=$(echo $myvar | grep -o 0 | wc -l)
    rc1=$(echo $myvar | grep -o 1 | wc -l)
    n01=$(echo $rc0'+'$rc1 | bc)
    if [ $n01 -eq $lenmyvar ] && [ $lenmyvar -gt 1 ]; then
       vk=B
       [ sbugga = 'true' ] && echo BINARY $myvar $vk
    else   
       vk=I
       [ sbugga = 'true' ] && echo INTEGER $myvar $vk
       if is_power_of_two "$myvar"; then
          true
          #printf "%s\n" "is power of 2"
	   fi
    fi 
elif [[ "$myvar" =~ ^[+-]?[0-9]+\.?[0-9]*$ ]]; then
    rc0=$(echo $myvar | grep -o 0 | wc -l)
    rc1=$(echo $myvar | grep -o 1 | wc -l)
    n01=$(echo $rc0'+'$rc1 | bc)
    if [ $n01 -eq $lenmyvar ]; then
       vk=B
       [ sbugga = 'true' ] && echo BINARY $myvar $vk
    else   
		vk=F
		[ sbugga = 'true' ] && echo FLOAT $myvar $vk
	fi	
elif [[ "$myvar" =~ ^[+-]?[0-9]+\.?[0-9]*[eE][+-]?[0-9]*$ ]]; then
    vk=E
    [ sbugga = 'true' ] && echo EXPONENTIAL $myvar
elif [[ "$myvar" =~ ^[+-]?[0-9]+\.?[0-9]*[+-]?[0-9]+\.?[0-9][iI]$ ]]; then
    vk=C
    [ sbugga = 'true' ] && echo COMPLEX $myvar  $vk
elif [[ "$myvar" =~ ^[+-]?[0-9]+\.?[0-9]*[+-]?[0-9]+[iI]$ ]]; then
    vk=C
    [ sbugga = 'true' ] && echo COMPLEX $myvar  $vk
elif [[ "$myvar" =~ ^[+-]?[0-9]*[+-]?[0-9][iI]$ ]]; then
    vk=C
    [ sbugga = 'true' ] && echo COMPLEX $myvar  $vk
fi    
[ sbugga = 'true' ] && echo 1 $vk

myvar1=$(echo $myvar | tr '[a-zA-Z]' '[A-ZA-Z]') # convert to uppercase for finding all possible values
if [ -n "$myvar1" ] && [ "$vk" = "0x0" ]; then
   [ sbugga = 'true' ] && echo 2 "$vk"
   if   [ "$myvar1" = "TRUE" ]  || [ "$myvar1" = "T"  ]; then
    vk=L
    [ sbugga = 'true' ] && echo LOGICAL $myvar  $vk
   elif [ "$myvar1" = "FALSE" ] || [ "$myvar1" = "F" ]; then
    vk=L
    [ sbugga = 'true' ] && echo LOGICAL $myvar  $vk
   elif [ "$myvar1" = "NAN" ]; then
    vk=NAN
    [ sbugga = 'true' ] && echo NAN $myvar  $vk
   elif [ "$myvar1" = "NULL" ]; then
    vk=NULL
    [ sbugga = 'true' ] && echo NULL $myvar  $vk
   elif [ "$myvar1" = "INF" ] || [ "$myvar1" = "+INF" ] || [ "$myvar1" = "-INF" ]; then
    vk=INF
    [ sbugga = 'true' ] && echo INF $myvar  $vk
   else
		if [[ "$myvar" =~ ^[0-9A-Fa-f]{1,}$ ]]; then
			vk=H
			[ sbugga = 'true' ] && echo HEX $myvar  $vk
		elif [[ "$myvar" =~ ^0x[0-9A-Fa-f]{1,}$ ]]; then
			vk=H
			[ sbugga = 'true' ] && echo HEX $myvar  $vk
		else
			vk=S
			[ sbugga = 'true' ] && echo STRG $myvar  $vk
		fi	
  fi
fi


if [ -z "$myvar" ]; then
    vk=ND
    [ sbugga = 'true' ] && echo NOTDEFINED $myvar  $vk
fi

if [ "$vk" = "0x0" ]; then
    vk=U
    [ sbugga = 'true' ] && echo UNRECOGNIZEDTYPE $myvar  $vk
fi

[ sbugga = 'true' ] && echo FINAL vk $vk
[ sbugga = 'true' ] && echo
[ sbugga = 'true' ] && echo
}




function bkpcode() {
echo "make a backup copy of the code"
cp -i $SCRIPT_DIR/ep3fsck.sh $SCRIPT_DIR/ep3fsck.$ver.sh
exit 1
}


function createmd5() {
echo "make a md5sum of the code"
md5sum  $SCRIPT_DIR/ep3fsck.sh >ep3fsck.$ver.md5sum
}

function checkmd5() {
md5sum --check ep3fsck.$ver.md5sum
}


function updatever() {
vernum=$(echo $vern + 0.1| bc)
newver=`date +%Y-%m-%dT%H:%M:%S`"_v"$vernum
echo
echo "ep3fsck new version timestamp:"
echo
echo "MB "$newver
echo
}

function editsrc() {
edname_set=0
for edname in efte geany kate gedit emacs nano vi
do
   which $edname  &>/dev/null 2>&1
   ret_code=$?
   if (( $ret_code == 0 )) && (( $edname_set == 0 )); then
      edname_set=1
      echo $edname $edname_set
      $edname $SCRIPT_DIR/ep3fsck.sh&
      exit 1
   fi
done
}



function bashcheck() {
currentver=${BASH_VERSION:0:5}
requiredver="4.2.1"
if [ "$(printf '%s\n' "$requiredver" "$currentver" | sort -V | head -n1)" = "$requiredver" ]; then 
       true
else
       echo "Less than ${requiredver}"
       echo "current          BASH version:" $currentver
       echo "minimum required BASH version:" $requiredver
       echo
       echo "ep3fsck cannot run"
       echo
       echo "exiting ..."
       echo
       exit 1
fi
}



function fitsswcheck() {
err=0
which dfits  &>/dev/null 2>&1
ret_code=$?
if [ $ret_code -gt 0 ]; then
   echo
   echo "***  ERROR "
   echo
   echo "dfits      not found in your PATH"
   echo 
   echo "please install it or modify your PATH variable"
   echo
   echo "PATH variable"
   echo $PATH
   echo
   echo "dfits sourcecode"
   echo ftp://ftp.eso.org/pub/qfits/qfits-6.2.0.tar.gz
   echo http://www.eso.org/~jpritcha/jFEROS-DRS/FEROS-DRS/src/dfits.c
   echo
   let err=$err+1
else
   mydfits=$(which dfits)
fi


which fitsverify &>/dev/null 2>&1
ret_code=$?
if [ $ret_code -gt 0 ]; then
   echo
   echo "***  ERROR "
   echo
   echo "fitsverify not found in your PATH"
   echo 
   echo "please install it or modify your PATH variable"
   echo
   echo "PATH variable"
   echo $PATH
   echo 
   echo "fitsverify sourcecode"
   echo https://heasarc.gsfc.nasa.gov/docs/software/ftools/fitsverify/
   echo
   let err=$err+1
else
   fitsver=$(which fitsverify)
   fitsver_num=$("$fitsver" _null_file_.fits 2>&1| grep "fitsverify" | awk '{printf"%s\n",$2}')
   #$fitsver _null_file_.fits 2>&1| grep "fitsverify" | awk '{printf"%s %s\n",$3,$4}' | sed -e 's/(//g' -e 's/)//g' -e 's/V//g'
fi

# deb astropy-utils fits2bitmap

which convert &>/dev/null 2>&1
ret_code=$?
if [ $ret_code -gt 0 ]; then
   convt=0
else
   convt=1
fi

which wget &>/dev/null 2>&1
ret_code=$?
if [ $ret_code -gt 0 ]; then
   wget=0
else
   wget=1
fi


if [ $err -gt 0 ]; then
   exit 1
fi
}


function autodoc () {
a2ps -o $SCRIPT_DIR/ep3fsck.ps $SCRIPT_DIR/ep3fsck.sh
ps2pdf  $SCRIPT_DIR/ep3fsck.ps $SCRIPT_DIR/ep3fsck.pdf
pdfview_set=0
for pdfview in evince okular xpdf gv
do
  which $pdfview  &>/dev/null 2>&1
  ret_code=$?
  if (( $ret_code == 0 )) && (( $pdfview_set == 0 )); then
     pdfview_set=1
     #echo $edname $edname_set
     $pdfview $SCRIPT_DIR/ep3fsck.pdf &>/dev/null 2>&1&
     exit 1
  fi
done
}



function docdownload() {
mkdir ep3_doc
echo "*************************************************************"
echo "*************************************************************"
echo "*************************************************************"
echo "*************************************************************"
echo "*************************************************************"
echo
echo
echo "Downloading documentation ..."
echo
echo
echo
if [ $wget -eq 1 ]; then
   wget https://www.eso.org/sci/observing/phase3/p3sdpstd.pdf
   wget https://archive.eso.org/cms/tools-documentation/dicb/ESO-044156_6DataInterfaceControlDocument.pdf
   wget https://fits.gsfc.nasa.gov/standard40/fits_standard40aa-le.pdf
   echo "done!"
   echo "Documentation was succesfully downloaded in the directory ep3_doc"
   echo
else
   echo "wget not found in your PATH"
   echo 
   echo "please install it or modify your PATH variable"
   echo
   echo "PATH variable"
   echo $PATH
   echo
fi
}



function basic_ckw (){
		i=-1
		l1end=$(grep -n ^END$ $hdr | sed -e 's/:END//g')
		#echo l1end $l1end
		for i2 in $l1end
		do
			let i=$i+1
			if [ $i -eq 0 ]; then 
			   i1=2
			   i2p=$i2
			else
			   let i1=$i2p+2
			fi
			let i3=$i2-$i1+1
			i2p=$i2
			#echo ext $i s $i1 e $i2 l $i3
			echo >> $ofile
			echo EXTENSION $i >> $ofile
			#echo ${basic_kwlist[@]}
			for bname in "${basic_kwlist[@]}";
			do
				kvalue_n=$(head -$i2 $hdr | tail -n $i3 | cut -d= -f1 | grep ^$bname | wc -l)
				kvalue=$(head -$i2 $hdr | tail -n $i3 | grep ^$bname\ | cut -d= -f2 | cut -d/ -f1)
				if [ ! -z "$kvalue" ]; then 
					echo $bname '=' $kvalue >>$ofile
				fi
				# if the keyword is present more than once, print all its values
				#seq 1 $kvalue_n | while read n
				#do
				#	kvalue=$(head -$i2 $hdr | tail -n $i3 | grep ^$bname | head -$n | tail -1)
				#	echo $kvalue >>$ofile
				#done
			done
			for bname in LONGSTR CONTINUE CDELT1 CDELT2 CDELT3 PC1_1 PC1_2 PC2_1 PC2_2 PC001001 PROV0 OBID0
			do
					kvalue_n=$(cut -d= -f1 $hdr | grep ^$bname | wc -l)
					if [ $kvalue_n -gt 0 ]; then
						echo ERROR KEYWORD $bname NOT ALLOWED !!!! 
						echo "\textcolor{red}{ERROR "${bname}" not allowed !!!}"  >>$ofile
					fi
			done
			if [ $i -eq 0 ]; then
				#check DATE
				kvalue_DATE=$(head -$i2 $hdr | tail -n $i3 | grep ^DATE\ | cut -d= -f2 | cut -d/ -f1)
				if [ -z "$kvalue_DATE" ]; then
						echo ERROR KEYWORD DATE NOT PRESENT !!!! 
						echo "\textcolor{red}{ERROR DATE not present !!!}"  >>$ofile
				else
					true
					#echo DATE present
				fi
				#check progid
				progid=$(grep ^XPROG_ID $hdr | cut -d= -f2 | cut -d/ -f1)
				lprogid=${#progid}
				if [ $lprogid -eq 0 ]; then
					progid=$(grep 'HIERARCH ESO OBS PROG ID' $hdr | cut -d= -f2 | cut -d/ -f1)
				fi
				if [ $lprogid -gt 0 ]; then
					progid1=$(echo $progid | sed -e "s/'//g" -e "s/ //g")
					progid2=${progid1%(*}
					echo PROG_ID $progid2
					if [ progid = "'MULTI'" ]; then
						true
						echo multiple PROGID
						# look for PROGID
					else
						true
						echo single PROGID
					fi 
					if [ $provflag = true ]; then
						wget -O flist.csv 'http://archive.eso.org/tap_obs/sync?REQUEST=doQuery&LANG=ADQL&MAXREC=200&FORMAT=csv&QUERY=SELECT%20dp_id,dp_cat,dp_tech,dp_type,telescope,instrument,filter_path,grat_path,gris_path,slit_path,ins_mode,date_obs,mjd_obs,exposure,object,target,ra,dec,ra_pnt,dec_pnt,tel_alt,tel_az,telescope,obs_mode%0aFROM%20dbo.raw%0awhere%20prog_id%20like%20%27'$progid2'%25%27' 2>/dev/null
						nfiles=$(wc -l flist.csv| awk '{print $1-1}')
						#cat flist.csv | awk -F, '{printf"%s\\\\ \n",$0}' | sed -e 's/,/\&/g' # convert the csv to latex table
						
						#echo nfiles $nfiles
					fi
					if [ $bibflag = true ]; then
						wget -O bib.csv 'https://telbib.eso.org/export.php?programid='$progid2'&total=9&rows=30000&fl=author,title,journal,volume,pages,year,bibcode,telescope,instrument,programid,datastatus,survey,sciver&export_format=csv' 2>/dev/null
						nbib=$(wc -l bib.csv| awk '{print $1-1}')
						awk '{ FPAT = "([^, ]+)|(\"[^\"]+\")" }{ print $7 }' bib.csv
						#echo nbib $nbib
					fi
					if [ $proginfo = true ]; then
						wget -O abstract.tmp 'http://archive.eso.org/wdb/wdb/eso/abstract/query?progid='$progid2 2>/dev/null
						r1=$(grep -n '<table width="90%"><tr><td align="justify">' abstract.tmp| cut -d: -f1)
						r2=$(grep -n ^'</td></tr></table>' abstract.tmp| cut -d: -f1)
						let r3=$r2-$r1+1
						title=$(head -$r2 abstract.tmp | tail -n $r3 | grep '<h1>' | sed -e 's/<h1>//g' -e 's/<\/h1>//g' -e 's/^[[:space:]]*//')
						pi=$(head -$r2 abstract.tmp | tail -n $r3 | grep PI | sed -e 's/<p>//g' -e 's/<b>//g' -e 's/<\/b>//g' -e 's/^[[:space:]]*//')
						abstract=$(head -$r2 abstract.tmp | tail -n $r3 | grep Abstract | sed -e 's/<p>//g' -e 's/<b>//g' -e 's/<\/b>//g' -e 's/^[[:space:]]*//')
					fi
				fi
			fi
		done
       telescop=$(grep ^TELESCOP $hdr | cut -d= -f2 | cut -d/ -f1)
       if [ telescop = "'MULTI'" ]; then
			true
			# look for TELESCi
	   else
			false
       fi 
       filter=$(grep ^FILTER $hdr | cut -d= -f2 | cut -d/ -f1)
       if [ filter = "'MULTI'" ]; then
			true
			# look for FILTERi
	   else
			false
       fi 
       instrume=$(grep ^INSTRUME $hdr | cut -d= -f2 | cut -d/ -f1)
       if [ instrume = "'MULTI'" ]; then
			true
			# look for INSTRUMi
	   else
			false
       fi 
       
       
}




function ckw (){

         if ($oaout); then
            echo
            echo "HEADER " $nh
            echo "HEADER " $nh >> $oafile
         fi

       # create a list of kw presence
       kwpresence_list=(`seq 1 130`)
       i=-1
       for item in "${kwpresence_list[@]}"
       do
         let i=$i+1
         kwpresence_list[$i]=0
       done

       echo "\begin{Verbatim}[commandchars=\\\\\{\}]"  >> $ofile
       i=-1
       j=0
       k=0
       v=0

       # make a loop over all the keywords, start with a for on the kwflags
       # and since it has the same number of elements of kwlist, for each iteration it increase the index i
       # to get the right label from kwlist and increase k to get the right number in table 8
       
       for vv in "${kwflag[@]}";
       do
         # start a loop over all the kw
         # i :: count the kw index
         # k :: is like i, but increased by 1 to respect the indexes in table 8 of sdps
         # v :: is the flag in the data type table

         hn=${vv:0:1}  # header to which apply the flag
         v=${vv:1:2}   # flag

         let i=$i+1
         let j=$v-1
         if [ $i -eq 10 ]; then
            true
         else
            let k=$k+1
         fi
    

         kname="${kwlist[$i]}"
         kflag="${kflag_legend[$j]}"
         kflag_e="${kflag_legend_expl[$j]}"
         # list of kw indexes that could have duplicates
         kmulval=(2 3 4 19 31 32 33 34 35 36 37 66 76 88 106 107 111 116 117 118 119 120 121 122 123 124 126 127)                                     
         kounter=0
         # check if the keyword could be multiple present
         for kk in "${kmulval[@]}";
         do
            if (( "$k" == "$kk" )); then
               kounter=1
            fi
         done
         # define the search pattern for single or multiple keywords
         if   [ $kounter -eq 1 ]; then
                search_kname="^""$kname""[1-9]"
         elif [ $kounter -eq 0 ]; then
                search_kname="^""$kname\b"
         fi

         #search pattern for cases with complex keywords
         if (( $k == 6 )); then
            search_kname="-e ^TELESCOP -e ^$kname[1-9]"
         fi
         if (( $k == 7 )); then
            search_kname="-e ^INSTRUME -e ^$kname[1-9]"
         fi
         if (( $k == 8 )); then
            search_kname="-e ^FILTER -e ^$kname[1-9]"
         fi
         if (( $k == 18 )); then
            search_kname="-e ^PROG_ID -e ^$kname[1-9]"
         fi
         if (( $k == 25 )); then
            search_kname="-e ^PROVXTN -e ^$kname[1-9]"
         fi

		

		#echo k $k kname $kname

		if [ $k -eq 10 ] && [ $kname = 'RA' ]; then
		   #echo kname $kname
		   RA=$(grep ^"$kname\b" $hdr | cut -d= -f2 | cut -d/ -f1|sed -e's/ //g')
		   lRA=${#RA}
			varkind RA
			vk10=$vk
			#echo vl $vk vk10 $vk10 RA $RA
			if [ $vk10 = "F" ]; then
				s1=$(echo $RA'>'0 | bc)
				s2=$(echo $RA'<'360 | bc)
				#echo s1 $s1 s2 $s2
			    if [ $s1 -lt 1 ]; then
			       echo "ERROR: RA<0 !!!"
			       echo RA = $RA
			    elif [ $s2 -lt 1 ]; then
			       echo "ERROR: RA>360 !!!"
			       echo RA = $RA
			    #else
			    #   echo RA ok
			    fi
			fi
		fi
		if [ $k -eq 10 ] && [ $kname = 'DEC' ]; then
		   #echo kname $kname
		   DEC=$(grep ^"$kname\b" $hdr | cut -d= -f2 | cut -d/ -f1|sed -e's/ //g')
		   lDEC=${#DEC}
			varkind DEC
			vk10=$vk
			#echo vl $vk vk10 $vk10 DEC $DEC
			if [ $vk10 = "F" ]; then
				s1=$(echo $DEC'>'-90 | bc)
				s2=$(echo $DEC'<'90  | bc)
				#echo s1 $s1 s2 $s2
			    if [ $s1 -lt 1 ]; then
			       echo "ERROR: DEC<-90 !!!"
			       echo DEC = $DEC
			    elif [ $s2 -lt 1 ]; then
			       echo "ERROR: DEC>90 !!!"
			       echo DEC = $DEC
			    #else
			    #   echo DEC ok
			    fi
			fi
		fi
         
		if (( $k == 14 )); then
		   EXPTIME=$(grep ^$kname $hdr | cut -d= -f2 | cut -d/ -f1|sed -e's/ //g')
		   lEXPTIME=${#EXPTIME}
			varkind EXPTIME
			vk14=$vk
			if [ $vk14 = "F" ] || [ $vk14 = "E" ] || [ $vk14 = "I" ]; then
				s1=$(echo $EXPTIME'>='0 | bc -l)
				#echo s1 $s1
			    if [ $s1 -lt 1 ]; then
			       echo "ERROR: EXPTIME<0 !!!"
			       echo EXPTIME = $EXPTIME
			    #else
			    #   echo EXPTIME ok
			    fi
			fi
		fi
		if (( $k == 15 )); then
		   TEXPTIME=$(grep ^$kname $hdr | cut -d= -f2 | cut -d/ -f1|sed -e's/ //g')
		   lTEXPTIME=${#TEXPTIME}
			varkind TEXPTIME
			vk15=$vk
			if [ $vk15 = "F" ] || [ $vk15 = "E" ] || [ $vk15 = "I" ]; then
				s1=$(echo $TEXPTIME'>'1E-15 | bc -l)
				s2=$(echo $TEXPTIME'<'1E15  | bc -l)
				#echo s1 $s1 s2 $s2
			    if [ $s1 -lt 1 ]; then
			       echo "ERROR: TEXPTIME<1E-15 !!!"
			       echo TEXPTIME = $TEXPTIME
			    elif [ $s2 -lt 1 ]; then
			       echo "ERROR: TEXPTIME>1E+15 !!!"
			       echo TEXPTIME = $TEXPTIME
			    #else
			    #   echo TEXPTIME ok
			    fi
				s1=$(echo $EXPTIME'<'$TEXPTIME+0.01| bc -l)
				#echo s1 $s1 
			    if [ $s1 -lt 1 ]; then
			       echo "ERROR: EXPTIME<TEXPTIME+0.01 !!!"
			       echo EXPTIME  = $EXPTIME
			       echo TEXPTIME = $TEXPTIME
			       echo EXPTIME - TEXPTIME: `echo $EXPTIME-$TEXPTIME|bc -l`      
			    #else
			    #   echo TEXPTIME ok
			    fi
			fi
		fi
		if (( $k == 16 )); then
		   MJDOBS=$(grep $kname $hdr | cut -d= -f2 | cut -d/ -f1|sed -e's/ //g')
		   lMJDOBS=${#MJDOBS}
			varkind MJDOBS
			vk16=$vk
			if [ $vk16 = "F" ]; then
				s1=$(echo $MJDOBS - 40587 | bc | grep -c '-')
			    if [ $s1 -ge 1 ]; then
			       echo "ERROR: MJD-OBS<40587 !!!"
			       echo MJD-OBS = $MJDOBS
			    #else
			    #   echo MJD-OBS ok
			    fi
			fi
		fi
		if (( $k == 17 )); then
		   MJDEND=$(grep $kname $hdr | cut -d= -f2 | cut -d/ -f1|sed -e's/ //g')
		   lMJDEND=${#MJDEND}
			varkind MJDEND
			vk17=$vk
		fi
		if [ $k -eq 17 ] && [ $lMJDOBS -gt 0 ] && [ $lMJDEND -gt 0 ]; then
			if [ $vk16 = "F" ] && [ $vk17 = "F" ]; then
			    s1=$(echo $MJDOBS - $MJDEND | bc | grep -c '-')
			    if [ $s1 -lt 1 ]; then
			       echo "ERROR: MJD-END<MJD-OBS !!!"
			       echo MJD-OBS = $MJDOBS
			       echo MJD-END = $MJDEND
			    #else
			    #   echo MJDs ok
			    fi
			fi
				s1=$(echo $TEXPTIME/86400.0'<='$MJDEND-$MJDOBS| bc -l)
				#echo s1 $s1 
			    if [ $s1 -lt 1 ]; then
			       echo "ERROR: TEXPTIME<=MJDEND-MJDOBS !!!"
			       echo TEXPTIME = $TEXPTIME
			       echo MJD-OBS  = $MJDOBS
			       echo MJD-END  = $MJDEND
			       echo MJD-OBS -MJD-END [s]: `echo '('$MJDEND-$MJDOBS')'*86400|bc -l`      
			    #else
			    #   echo TEXPTIME'<='MJDEND-MJDOBS ok
			    fi
		fi
         
         # for each kw count how many occurrences of the kw are present
         kvalue_n=$(cut -d= -f1 "$hdr" | grep -c -E $search_kname )


#1 = MM = MANDATORY
#2 = ME = MANDATORY ESO
#3 = MA = MANDATORY WHEN APPLICABLE
#4 = RC = RECOMMENDED
#5 = NO = NOT ALLOWED
#6 = OP = OPTIONAL
#7 = na = NOT APPLICABLE
#8 = un = UNKNOWN
#9 = RE = RESERVED
#0 = SU = ASK SUPPORT
         

         # if the keyword is present more than once, print all its values
#         if [ $kvalue_n -ge 0 ]; then
            if [ $kvalue_n -eq 0 ]; then
               maxval=1
            else
               maxval=$kvalue_n
            fi
            if [ $kvalue_n -ge 2 ]; then
               multiflag="\textcolor{red}{>1 times}"
            else
               multiflag=""
            fi
            seq 1 $maxval | while read n
            do 
              kvalue=$(grep $search_kname $hdr | head -"$n" | tail -1 | cut -d= -f2| cut -d/ -f1|awk '{$1=$1};1')
              #kvalue=$(grep $search_kname $hdr | head -"$n" | tail -1 | cut -d= -f2)
              #kvalue=$(cat "$hdr" | grep "$search_kname" | head -"$n" | tail -1 | cut -d= -f2)
              #kvalue=$(cat "$hdr" | grep $search_kname | head -"$n" | tail -1)
              #echo $k $maxval $search_kname $kvalue
              kvalue_len=${#kvalue}
              if [ $p3hdr = "true" ]; then
                 if [ $kvalue_len -ge 1 ]; then
					kvalue_full=$(grep $search_kname $hdr | head -"$n" | tail -1)
				 else
					kvalue_full=$kname" = unset value [written by ep3fsck]"
                 fi
                 if (( $v == 1)) || (( $v == 2)) || (( $v == 3)) || (( $v == 4)) || (( $v == 6)) || (( $v == 9)); then
                    if (( $nh == 1 )); then
					   echo $kvalue_full >> $p3hdrfile
					fi
                    if (( $nh == 2 )) && (( $hn >=2 )); then
					   echo $kvalue_full >> $p3hdrfile
					fi
                 fi
              fi
              
              #echo k $k kname $kname kvalue $kvalue
              varkind kvalue
              #echo '##################'
              if [ $vk = "C" ] || [ $vk = "H" ] || [ $vk = "U" ] || [ $vk = "NAN" ] || [ $vk = "NULL" ] || [ $vk = "INF" ]; then
                 echo $kname "=" $kvalue "(BAD KIND:" $vk")"
              fi
              if [ $vk = "NOTDEF" ]; then
                 kwnotdef=true
              fi
              if [ $vk = "F" ] || [ $vk = "E" ] || [ $vk = "I" ] ; then
                 kwnumb=true
              fi
              if [ $vk = "L" ]; then
                 kwlog=true
              fi
              if [ $vk = "S" ]; then
                 kwstr=true
                 #echo STRING $kname $kvalue
              fi
              ##echo ext $nh num $k where $hn flag $v len $kvalue_len
              
              if [ $kvalue_len -lt 1 ]; then
                    if (( $v == 7)) || (( $v == 8)) || (( $v == 9)) || (( $v == 0)); then
                      kvalue=$kflag_e
                      kvalue_desc="kw absent"
                      printyn=false
                    else
                      kvalue="\textcolor{red}{"$kflag_e"}"
                      kvalue_desc="\textcolor{red}{kw absent}"
                      printyn=true
                      #if [ $p3hdr = "true" ]; then
						#kvalue_full=$kname" = unset value [written by ep3fsck]"
						#echo $kvalue_full >> $p3hdrfile
					  #fi
                      if (( $v == 5 )); then
                         printyn=false
                      fi
                    fi
              else
                    if   (( $v == 7)) || (( $v == 8)); then
                       kvalue="\textcolor{orange}{WARNING kw not relevant}"
                       kvalue_desc="\textcolor{orange}{kw present}"$multiflag
                       printyn=false
                    elif (( $v == 5)); then
                       kvalue="\textcolor{red}{ERROR kw not allowed}"
                       kvalue_desc="\textcolor{red}{kw present}"$multiflag
                       printyn=true
                    elif (( $v == 9)); then
                       kvalue="\textcolor{magenta}{WARNING kw reserved}"
                       kvalue_desc="\textcolor{magenta}{kw present}"$multiflag
                       printyn=true
                    elif (( $v == 0)); then
                       kvalue="\textcolor{brown}{WARNING ask support}"
                       kvalue_desc="\textcolor{brown}{kw present}"$multiflag
                       printyn=true
                    else
                       kvalue_desc="\textcolor{green}{kw present}"$multiflag
                       printyn=false
                    fi
              fi
              
              if (( "$nh" != "$hn" )) &&  (( "$hn" < 3 )); then
                 kvalue="\textcolor{teal}{kw not relevant}"
                 if [ $kvalue_len -lt 1 ]; then
                    kvalue_desc="\textcolor{teal}{kw absent}"
                 fi
                 printyn=false
              fi
              
              strg1=$(printf "%3.3i | %s%s | %s | %s\n" "$k" "$kname" "${padding:${#kname}}" "$kvalue" "$kvalue_desc" | sed -e 's/_/\\_/g')
              if ($printyn); then
                echo $strg1 >>$ofile
				if ($oaout); then
					echo $strg1 | sed -e 's/\\textcolor{red}//g' -e 's/{//g' -e 's/}//g' -e 's/\\_/_/g'  
					echo $strg1 | sed -e 's/\\textcolor{red}//g' -e 's/{//g' -e 's/}//g' -e 's/\\_/_/g'  >>$oafile
				fi
              fi
            
            done
#         fi
       done
       echo "\end{Verbatim}"      >>$ofile
}


function readme() {
echo $"
HELP
ep3fsck.sh -h

single FITS
ep3fsck.sh -T out.txt -f fitsfile -t dataprod

all FITS in the current directory
ep3fsck.sh -T out.txt -f all -t dataprod

dataprod:
         image  =   IMAGE (single) associated file format
         imef   =   IMAGE (single) MEF format
         mefi   =   MEFIMAGE (mosaic) mef format
         iflux  =   IMAGE.FLUXMAP
         spec   =   SPECTRUM
         cube   =   CUBE.IFS
         visi   =   VISIBIlITY
         srctbl =   SRCTBL
         mcat   =   MCATALOG
         ctile  =   CATALOGTILE
         catal  =   CATALOG
"
}


function help() {
echo $"
Usage: ep3fsck OPTIONS
 -f fnames    : file names, possible values                                  
    all       = analyze all fits in the directory
    name.fits = analyze only this FITS file
 -h           : print this help and exit
 -T           : plain ASCII output file                        [out.txt]                 
 -t dataprod  : data type, possible values                                   
    image     = IMAGE    (single) file format
    imef      = IMAGE    (single) MEF  format
    mefi      = MEFIMAGE (mosaic) MEF  format
    iflux     = IMAGE.FLUXMAP
    spec      = SPECTRUM
    cube      = CUBE.IFS
    visi      = VISIBIlITY
    srctbl    = SRCTBL
    mcat      = MCATALOG
    ctile     = CATALOGTILE
    catal     = CATALOG
 -p name.hdr  : intendend for ESO ADP files,
                print only the phase3 SDPS headers
"
exit 1
}

function devhelp() {
echo $"
Usage: ep3fsck OPTIONS
(in squared brackets the default values)
 -a           : enable print of all headers   (true/false)       [false]
 -c           : collection name in the LaTeX report             [Report]
 -f fnames    : file names, possible values                                  
    all       = analyze all fits in the directory
    name.fits = analyze only this FITS file 
 -h           : print this help and exit
 -k           : print in the LaTeX report the Table 8 of SDPS
                                              (true/false)       [false]
 -l           : LaTeX output                                   [out.tex]
 -m           : download ESO phase 3 and FITS documentation                  
 -T           : plain ASCII output file                        [out.txt]                 
 -r           : print some usage examples                                    
 -s           : enable Simbad search on the object  (true/false) [false]
 -t dataprod  : data type, possible values                                   
    image     = IMAGE    (single) file format
    imef      = IMAGE    (single) MEF  format
    mefi      = MEFIMAGE (mosaic) MEF  format
    iflux     = IMAGE.FLUXMAP
    spec      = SPECTRUM
    cube      = CUBE.IFS
    visi      = VISIBIlITY
    srctbl    = SRCTBL
    mcat      = MCATALOG
    ctile     = CATALOGTILE
    catal     = CATALOG
 -p name.hdr  : intendend for ESO ADP files,
                print only the phase3 SDPS headers
 -u           : user name in the LaTeX report                   [\$USER]
DEVEL OPTIONS
 -B           : make version backup
 -C           : check  md5sum
 -D           : show source code                                     
 -E           : edit the source code
 -M           : create md5sum
 -V           : provide new version timestamp
How to update software version and make the new checksum:
1) Get the new version and timestamp   ./ep3fsck.sh -V
2) Manually update the source code 
3) create the new checksum             ./ep3fsck.sh -M
4) check if the checksum correspond    ./ep3fsck.sh -C
5) make a backup copy of the software  ./ep3fsck.sh -B
"
exit 1
}



###############################################################################
#                                                                             #
#                                                                             #
# M A I N   C O D E                                                           #
#                                                                             #
#                                                                             #
###############################################################################

logo

numarg=$#
if [ $numarg -eq 0 ]; then
  help
fi




while getopts ":p:a:c:f:k:T:o:t:u:h:H:m:r:s:B:C:D:E:M:V:" opt; do
  OPTARG_len=${#OPTARG}
  case $opt in
    a) allheaders="$OPTARG"
    ;;
    c) collection="$OPTARG"
    ;;
    f) files="$OPTARG"
    if [ $files == all ]; then
       files="*.fits"
    fi
    numfiles=$(ls $files | wc -l)
    if [ $numfiles -lt 1 ]; then
       echo "*** ERROR: no such files present!"
       echo "exiting ..."
       exit 1
    fi
    ;;
    o) 
       ofile="$OPTARG"
       opdf="${ofile/.tex/.pdf}"
       latex_=true
       rm -rf $ofile $opdf
    ;;
    p)
       p3hdr="true"
       p3hdrfile="$OPTARG"
       rm -rf $p3hdrfile
    ;;   
    T) 
       oaout="true"
       oafile="$OPTARG"
       rm -rf $oafile 
    ;;
    s) simbad='true'
    ;;
    k) printkwlist='true'
    ;;
    t) tipo="$OPTARG"
    ;;
    u) user="$OPTARG"
    ;;
    \?) echo "*** WARNING: Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac


  case $OPTARG in
    h) 
    help
    exit 1
    ;;
    m)
    docdownload
    exit 1
    ;;
    r)
    readme
    exit 1
    ;;
    H)
    devhelp
    exit 1
    ;;
    B)
    bkpcode
    exit 1
    ;;
    C)
    checkmd5
    exit 1
    ;;
    D)
    autodoc
    exit 1
    ;;
    E)
    editsrc
    exit 1
    ;;
    F)
    latex_=true
    simbad=true
    fitsverify_=true
    fileinfo_=true
    image_=true
    ;;
    M)
    createmd5
    exit 1
    ;;
    V)
    updatever
    exit 1
    ;;
  esac

  case $OPTARG in
    -*)
    echo >&2 'ep3fsck:' $"unrecognized option" "\`$1'"
    echo >&2 $"Try \`ep3fsck -h' for more information."
    exit 1
    ;;
  esac
done



#check the bash version
bashcheck

#check for presence of dfits and fitsverify
fitsswcheck


# if no filenames are provided, ask the user
if [ ${#files} -eq 0 ]; then
   echo
   echo "*** WARNING: no filenames provided!"
   echo "Do you want to use all FITS files in the current directory ? [y/n]"
   read useall
   case $useall in
   y|Y) 
        files="*.fits"
        numfiles=$(ls $files | wc -l)
        if [ $numfiles -lt 1 ]; then
           echo "*** ERROR: no such files present!"
           echo "exiting ..."
           exit 1
        fi
        ;;
   n|N) echo "Do you want to provide a filename ? [y/n]"
        read provfilename
        case $provfilename in
           y|Y)
           echo "insert filename"
           read files
           numfiles=$(ls $files | wc -l)
           if [ $numfiles -lt 1 ]; then
              echo "*** ERROR: no such files present!"
              echo "exiting ..."
              exit 1
           fi
           ;;
           n|N)
           echo "exiting ..."
           exit 1
           ;;
        esac
        ;;
   esac
fi


# if no dataprod is provided, ask the user
if [ ${#tipo} -eq 0 ]; then
   echo
   echo "*** WARNING: no file dataprod provided!"
   echo echo "Do you want to provide a file dataprod ? [y/n]"
   read usetipo
   case $usetipo in
   y|Y) echo $"
supported dataprod:
         image  =   IMAGE (single) associated file format
         imef   =   IMAGE (single) MEF format
         mefi   =   MEFIMAGE (mosaic) MEF format
         iflux  =   IMAGE.FLUXMAP
         spec   =   SPECTRUM
         cube   =   CUBE.IFS
         visi   =   VISIBIlITY
         srctbl =   SRCTBL
         mcat   =   MCATALOG
         ctile  =   CATALOGTILE
         catal  =   CATALOG
"
        read tipo
        ;;
   n|N) 
   echo "exiting ..."
   exit 1
   ;;
   esac
fi



# check the dataprod
co=0
for dt in "${data_types_short[@]}";
do
  if [ "$tipo" = "$dt" ]; then
      let co=$co+1
  fi
done
if   [ $co -gt 1 ]; then
     echo 'more than one dataprod correspondance'
     echo 'exiting ...'
     exit
elif [ $co -lt 1 ]; then
     echo 'no dataprod correspondance'
     echo 'exiting ...'
     exit
fi


#TABLE NOTES
#image 
# 2  3  4 11 13 24 91 96 102 103 104 105 114 125 
#31 32 31 37 38 40 49 51  52  52  52  52  59  65
#imagemef
# 2  3  4 11 13 24 79 80 81 82 83 84 85 91 96 102 103 104 105 114 125 
#31 32 31 37 38 40 46 46 46 46 46 46 46 49 51  52  52  52  52  59  65
#mefimage == image
# 2  3  4 11 13 24 91 96 102 103 104 105 114 125 
#31 32 31 37 38 40 49 51  52  52  52  52  59  65
#fluxmap 
# 2  3  4 11 13 24 95 114 125 
#31 33 31 37 38 40 50  59  65
#spectra
# 2  3  4 11 13 24 42 45 96 114 119 120 122 123 125
#31 32 31 37 38 40 42 43 51  59  61  62  63  63  65
#cube
# 2  3  4 11 13 24 79 80 81 95 96 114 125
#31 34 31 37 38 40 46 46 46 50 51  59  65
#visibility
# 2  3  4 11 13 24 114 119 125 
#31 33 31 37 38 40  59  60  65
#srctbl
#2 3 4 11 13 21 24 38 91 96 114 119 125
#31 32 31 37 38 39 40 49 51  59  60  65
#mcatalog
# 8 13 21 24 46 47 86 88 106 107 109 110 111 112 113 114 119 124 125
#35 38 39 40 44 44 47 48  53  53  55  56  57  57  57  59  61  64  65
#catalogtile
# 2  3  4  8 13 21 24 46 47 86 88 109 110   111 114 119 124 125
#31 32 31 35 38 39 40 44 44 47 48  55  56 57+58  59  61  64  65
#catalog
# 2  3  4  8 10  11   12 13 21 24 46 47 64 65 66 86 88 94 95 106 107 108 109 110 111   112 113 114 119 124 125
#31 32 31 35 36 36+37 36 38 39 40 44 44 45 45 45 47 48 45 50  54  54  54  55  56 57+58  57  57  59  61  64  65



# define the keywords flag for each data prod
if   [ $tipo == "image" ]; then
     kwflag=(11 13 13 13 11 11 11 11 11 11 11 13 11 13 11 11 11 11 12 12 12 11 11 11 11 12 11 16 16 16 16 11 11 11 11 11 14 14 11 14 11 17 17 17 17 17 17 17 17 17 17 17 17 15 15 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 10 17 17 17 17 17 17 17 17 17 17 17 11 17 17 11 11 11 17 17 17 17 13 19 31 31 19 19 13 13 13 13 17 17 17 17 17 17 17 17 13 17 17 17 17 17 17 17 17 17 17 23 18 18 18 18)
elif [ $tipo == "imef" ];  then
     kwflag=(11 13 13 13 11 11 11 11 11 11 11 13 11 13 11 11 11 11 12 12 12 11 11 11 11 12 21 16 16 16 16 21 21 21 21 21 24 24 21 24 21 17 17 17 17 17 17 17 17 17 17 17 17 15 15 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 10 21 21 21 21 23 23 23 23 23 23 23 21 17 17 21 21 21 17 17 17 17 13 19 31 31 19 19 13 13 13 13 17 17 17 17 17 17 17 17 13 17 17 17 17 17 17 17 17 17 17 23 18 18 18 18)
elif [ $tipo == "mefi" ];  then
     kwflag=(11 13 13 13 11 11 11 11 11 11 11 13 11 13 11 11 11 11 12 12 12 11 11 11 11 12 21 16 16 16 16 21 21 21 21 21 24 24 21 24 21 17 17 17 17 17 17 17 17 17 17 17 17 15 15 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 10 17 17 17 17 17 17 17 17 17 17 17 21 17 17 21 21 21 17 17 17 17 13 19 31 31 19 19 13 13 13 13 17 17 17 17 17 17 17 17 13 17 17 17 17 17 17 17 17 17 17 23 18 18 18 18)
elif [ $tipo == "iflux" ]; then
     kwflag=(11 13 13 13 11 11 11 11 11 11 11 13 11 13 11 11 11 11 12 12 12 11 11 11 11 12 11 17 17 17 17 11 11 11 11 11 14 14 17 17 17 17 17 17 17 11 11 11 17 17 17 17 17 15 15 17 17 17 17 17 17 17 17 17 11 11 11 17 17 17 17 17 17 17 10 17 17 17 17 17 17 17 17 17 17 17 15 17 17 11 11 11 17 17 11 13 17 19 31 31 19 19 17 17 17 17 17 17 17 17 17 17 17 17 13 17 17 17 17 17 17 17 17 17 17 23 18 18 18 18)
elif [ $tipo == "spec" ];  then
     kwflag=(11 13 13 13 11 11 11 17 31 31 31 13 11 13 11 11 11 11 12 12 12 11 11 11 11 12 15 17 17 17 17 17 17 17 17 15 17 17 17 17 17 11 13 11 11 11 11 11 16 16 11 14 14 14 14 21 21 21 21 21 21 21 21 21 17 17 17 17 17 17 17 17 17 17 10 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 11 11 17 17 13 19 31 31 19 19 17 17 17 17 17 17 17 17 17 17 17 17 13 21 21 21 26 21 21 21 23 23 17 23 18 18 18 18)
elif [ $tipo == "cube" ];  then
     kwflag=(11 13 13 13 11 11 11 17 11 11 11 13 11 13 11 11 11 11 12 12 12 11 11 11 11 12 21 17 17 17 17 21 21 21 21 21 24 24 27 17 17 11 17 17 17 17 11 11 17 17 15 17 17 15 15 27 27 27 27 27 27 27 27 27 17 17 17 17 17 17 17 17 17 17 10 21 21 21 21 23 23 23 17 17 17 17 11 11 17 17 17 16 17 11 11 13 13 19 31 31 19 19 17 17 17 17 17 17 17 17 17 17 17 17 13 27 27 27 27 27 27 27 27 27 17 27 18 18 18 18)
elif [ $tipo == "visi" ];  then
     kwflag=(11 13 13 13 11 11 11 17 11 11 11 13 11 13 11 11 11 11 12 12 12 11 17 11 11 12 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 11 17 17 17 17 11 11 17 17 17 11 11 17 17 17 17 17 17 17 17 17 17 17 17 17 17 11 11 11 11 11 11 11 10 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 11 17 17 17 19 31 31 19 19 17 17 17 17 17 17 17 17 17 17 17 17 13 21 21 21 26 23 17 26 26 26 26 23 18 18 18 18)
elif [ $tipo == "srctbl" ];then
     kwflag=(11 13 13 13 11 11 11 11 11 11 11 13 11 13 11 11 11 11 12 12 12 11 17 11 11 12 17 17 17 17 17 17 17 17 17 17 17 17 23 24 21 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 10 17 17 17 17 17 17 17 17 17 17 17 21 17 17 21 21 21 17 17 17 17 13 19 31 31 19 19 17 17 17 17 17 17 17 17 17 17 17 17 13 21 21 21 26 23 17 26 26 26 26 23 18 18 18 18)
elif [ $tipo == "mcat" ];  then
     kwflag=(11 17 17 17 11 11 11 13 11 17 17 17 17 13 17 17 11 11 12 12 17 13 17 11 11 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 11 17 17 17 17 17 11 11 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 10 17 17 17 17 17 17 17 17 17 17 17 13 17 13 17 17 17 17 17 17 17 17 19 31 31 19 19 17 17 17 17 13 13 11 13 13 13 13 13 13 21 21 21 21 21 17 21 26 26 23 23 15 15 15 15)
elif [ $tipo == "ctile" ]; then
     kwglag=(11 13 13 13 11 11 11 13 11 11 11 13 11 13 17 17 11 11 12 12 17 13 17 11 11 12 17 17 17 17 17 17 17 17 17 17 17 17 17 17 11 17 17 17 17 17 11 11 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 10 17 17 17 17 17 17 17 17 17 17 17 13 17 13 17 17 17 17 17 17 17 17 19 31 31 19 19 17 17 17 17 11 11 11 13 13 13 17 17 13 21 21 21 21 21 17 21 26 26 23 23 15 15 15 15)
elif [ $tipo == "catal" ]; then
     kwflag=(11 13 13 13 11 11 11 13 11 13 13 13 13 13 17 17 11 11 12 12 17 13 17 11 11 12 17 17 17 17 17 17 17 17 17 17 17 17 17 17 11 17 17 17 17 17 11 11 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 13 13 13 17 17 17 17 17 17 17 10 17 17 17 17 17 17 17 17 17 17 17 13 17 13 17 17 17 17 17 13 13 17 19 31 31 19 19 17 17 17 17 13 13 13 13 13 13 13 13 13 21 21 21 21 21 17 21 26 26 23 23 15 15 15 15)
else
  echo "unsupported data file"
  echo "exiting"
  exit
fi






# create the latex base document
cat >$ofile<<EOFA
\documentclass[a4paper,10pt]{article}
\usepackage{graphicx}
\usepackage{fancyvrb}
\usepackage{fvextra}
\usepackage{xcolor}
\usepackage{hyperref}
\hypersetup{
    colorlinks=true,
    linkcolor=cyan,
    filecolor=magenta,
    urlcolor=blue,
    }
\urlstyle{same}
EOFA

echo "\title{"$collection"}">>  $ofile
echo "\author{"$user"}">>  $ofile



cat >>$ofile<<EOFB
\date{\today}
\begin{document}
\maketitle
\tableofcontents
\clearpage
\section*{Relevant Documents}
ESO Phase 3 standards documentation (SDPS)\\\\
\url{https://www.eso.org/sci/observing/phase3/p3sdpstd.pdf}\\\\
~\\\\
ESO Data Interface Control\\\\
\url{https://archive.eso.org/cms/tools-documentation/eso-data-interface-control.html}\\\\
~\\\\
FITS standard:\\\\
\url{https://fits.gsfc.nasa.gov/fits_standard.html}\\\\
\clearpage
EOFB

if [ $printkwlist = 'true' ]; then

echo "\section*{Keywords list}" >>  $ofile


echo "\begin{Verbatim}[breaklines=true]"  >> $ofile
printf "| %s | %s | %s | %s | %s | \n" Idx KW TYPE HDR FLAG  >>  $ofile
seq 0 129 | while read i
do
   #echo $i
   let j=$i+1
   printf "| %3.3i | %-8s | %1s | %1i | %1i |\n" $j "${kwlist[$i]}" "${kwtype[$i]}" "${kwflag[$i]:0:1}" "${kwflag[$i]:1:2}" >>  $ofile
done
echo "\end{Verbatim}"  >> $ofile
echo "\clearpage" >>  $ofile

echo "\begin{Verbatim}[breaklines=true]"  >> $ofile
cat >>$ofile<<EOD1
KW FLAG values
1  = MANDATORY
2  = MANDATORY ESO
3  = MANDATORY WHEN APPLICABLE
4  = RECOMMENDED
5  = NOT ALLOWED
6  = OPTIONAL
7  = NOT APPLICABLE
8  = UNDEFINED
9  = RESERVED
0  = CONTACT SUPPORT
EOD1
echo "\end{Verbatim}"  >> $ofile

echo "\begin{Verbatim}[breaklines=true]"  >> $ofile
cat >>$ofile<<EOD2
HDR values
1  = PRIMARY HEADER
2  = EXTENSIONS HEADER
3  = PRIMARY AND EXTENSIONS HEADER
EOD2
echo "\end{Verbatim}"  >> $ofile

echo "\begin{Verbatim}[breaklines=true]"  >> $ofile
cat >>$ofile<<EOD3
TYPE values
L  = LOGICAL
F  = FLOAT
I  = INTEGER
S  = STRING
EOD3
echo "\end{Verbatim}"  >> $ofile

fi

echo "\cleardoublepage" >>  $ofile






# start the check

nfiles=$(ls -1 $files |  wc -l)

nff=0
ls -1 $files | while read f
do
   let nff=$nff+1
   echo "checking" "$f" "("$nff" of "$nfiles")"
   g=$(echo "$f" | sed -e 's/_/\\_/g')  
   echo "\section{"$g"}" >>  $ofile
   $mydfits -x 0 "$f" > $hdr

         if ($oaout); then
            echo "$f" >> $oafile
         fi

   
   # get SIMBAD information
   if [ $simbad == true ]; then
      echo "\subsection{SIMBAD object information}" >>  $ofile
      objname=$(grep ^OBJECT $hdr | head -1 | cut -d= -f2 | cut -d/ -f1 | sed -e "s/'//g" -e 's/"//g' -e 's/ //g')
      if [ ${#objname} -gt 0 ]; then
         ogif=$objname".gif"
         ojpg=$objname".jpg"
      else
         ogif=dump.gif
         ojpg=dump.jpg
      fi
      RA=$(grep ^RA $hdr | grep -v RADECSYS | grep -v RADESYS | grep -v RA_ |  head -1 | cut -d= -f2 | cut -d/ -f1 | sed -e "s/'//g" -e 's/"//g' -e 's/ //g')
      DEC=$(grep ^DEC $hdr | head -1 | cut -d= -f2 | cut -d/ -f1 | sed -e "s/'//g" -e 's/"//g' -e 's/ //g')
      if (( ${#RA} > 0 )) && (( ${#DEC} > 0 )); then
      if [ $wget -eq 1 ]; then
         #simbad query by URL
         ##wget -O $simb 'http://simbad.cds.unistra.fr/simbad/sim-coo?Coord='$RA'+'$DEC'&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=30&Radius.unit=arcsec&submit=submit+query&%20id&output.format=ASCII&output.max=1' 2>/dev/null
         wget -O $simb 'http://simbad.cds.unistra.fr/simbad/sim-coo?Coord='$RA'+'$DEC'&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=3&Radius.unit=arcmin&submit=submit+query&%20id&output.format=votable&output.max=1&output.params=main_id,ids,otype(V),coo(d;ICRS),coo(d;GAL),coo(d;ECL),pmra,pmdec,rv_value,z_value,sp,mt,dim_majaxis,dim_minaxis,dim_angle,flux(B),flux(V),angdist' 2>/dev/null
         # get DSS cutout 
         wget -O $ogif 'archive.eso.org/dss/dss/image?ra='$RA'&dec='$DEC'&x=10&y=10&Sky-Survey=DSS1&mime-type=download-gif' &>/dev/null 2>&1
         if [ $convt -eq 1 ]; then
            convert $ogif $ojpg &>/dev/null 2>&1
            echo '\includegraphics[width=7cm]{'"$ojpg"'}\\\\' >>$ofile
         else
            echo "* WARNING convert not found in your PATH, please install imagemagick or modify your PATH variable" >> $ofile
            echo "%\includegraphics[width=7cm]{""$ogif""}\\\\ %don't work with gif files" >>$ofile
         fi
         echo "\begin{Verbatim}[breaklines=true]"  >> $ofile
         echo 'OBJECT =' $objname >> $ofile
         #echo 'RA     =' $RA      >> $ofile
         #echo 'DEC    =' $DEC     >> $ofile
         echo                     >> $ofile
         echo                     >> $ofile
         #cat $simb 2>&1 |sed -e 's/\#/*/g' -e 's/_/\_/g' | grep . |awk '{printf"%s \n",$0}'>>$ofile
         grep '<TR>' $simb | sed -e 's/<TR><TD>//g' -e 's/<\/TD><\/TR>//g' -e 's/<\/TD><TD>/,/g' | awk -F, '{printf"Main ID: %s\nother IDs: %s\ndescr: %s\nRA: %s deg\nDEC: %s deg \nGLON: %s deg \nGLAT: %s deg\nELON: %s deg \nELAT: %s deg \npmRA: %s mas/yr\npmDEC: %s mas/yr\nRV: %s km/s\nredshift: %s\nspectraltype: %s\nmorphologicaltype: %s\nsemimajorax: %s arcsec\nsemiminorax: %s arcsec\nposangle: %s deg\nBmag: %s\nVmag: %s\nangdist: %s arcsec\n",$1,$2,$3,$4,$5,$9,$10,$14,$15,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30 }'>>$ofile
      else
         echo "* WARNING wget not found in your PATH, please install it or modify your PATH variable" >> $ofile
      fi
      echo "\end{Verbatim}"  >> $ofile
      #less $simb
      else
         echo "* WARNING coordinates not defined !!!"
      fi
      echo "\clearpage"  >> $ofile
   fi

   # get info on the file
   if [ "$tipo" == 'image' ] || [ "$tipo" == 'imef' ] || [ "$tipo" == 'mefi' ] || [ "$tipo" == 'iflux' ] || [ "$tipo" == 'visi' ]; then
      echo "\subsection{File content}" >>  $ofile
      if [ $convt -eq 1 ]; then
         oojpg="$f".jpg
         convert -equalize "$f" "$oojpg" &>/dev/null 2>&1
         echo '\includegraphics[width=7cm]{'"$oojpg"'}\\\\' >>$ofile
      fi
      echo "\clearpage" >>  $ofile
   fi


   # get info on the file
   if [ $fileinfo = true ]; then
		echo "\subsection{File description}" >>  $ofile
		echo "\begin{verbatim}"  >> $ofile
		rpath=$(readlink -f "$f"|sed -e 's/_/\_/g')
		lsla=$(ls -la "$f"|sed -e 's/_/\_/g')
		lsatt=$(lsattr "$f"|sed -e 's/_/\_/g')
		dush=$(du -sh "$f"|sed -e 's/_/\_/g')
		filet=$(file "$f"|sed -e 's/_/\_/g')
		md5s=$(md5sum "$f"|sed -e 's/_/\_/g')
   
		echo "hostname"  >> $ofile
		hostname -a      >> $ofile
		echo             >> $ofile
		echo "date"      >> $ofile
		date             >> $ofile
		echo             >> $ofile
		echo "bash version">> $ofile
		echo $BASH_VERSION >> $ofile
		echo             >> $ofile
		echo "file path" >> $ofile
		echo $rpath      >> $ofile
		echo             >> $ofile
		echo "ls -la   " >> $ofile
		echo $lsla       >> $ofile
		echo             >> $ofile
		echo "lsattr   " >> $ofile
		echo $lsatt      >> $ofile
		echo             >> $ofile
		echo "du -sh   " >> $ofile
		echo $dush       >> $ofile
		echo             >> $ofile
		echo "filetype " >> $ofile
		echo $filet      >> $ofile
		echo             >> $ofile
		echo "md5sum   " >> $ofile
		echo $md5s       >> $ofile
		echo             >> $ofile
		echo "stat     " >> $ofile
		stat "$f"|sed -e 's/_/\\_/g' | awk '{printf"%s \n",$0}' >> $ofile
		echo "\end{verbatim}"  >> $ofile
		echo "\clearpage"  >> $ofile
	fi

   # run fitsverify
   echo "\clearpage"  >> $ofile
   echo "\subsection{Fitsverify}" >>  $ofile
   echo "\begin{verbatim}">> $ofile
   $fitsver "$f" 2>&1 |sed -e 's/\#/*/g' -e 's/_/\_/g' | awk '{printf"%s \n",$0}'>>$ofile
   echo "\end{verbatim}"  >> $ofile
   echo "\clearpage"  >> $ofile


   # get basic header information
   echo "\subsection{Basic headers informations}" >>  $ofile
   # get the number of extensions
   #ne=$(hexdump -e '80/1 "%_p""\n"' "$f" | grep -c  ^XTENSION) # works well if the files are small
   ne=$(grep -c ^XTENSION $hdr)
   echo "\begin{verbatim}">> $ofile
   echo                   >> $ofile
   echo "number of extensions found :" $ne >> $ofile
   echo                   >> $ofile
   basic_ckw
   echo "\end{verbatim}"  >> $ofile
   echo "\clearpage"  >> $ofile
   if [ $proginfo = true ]; then
		echo "\subsection{Program info}" >>  $ofile
		echo $title"\\\\" >> $ofile
		echo "\\\\" >> $ofile
		echo $pi"\\\\" >> $ofile
		echo "\\\\" >> $ofile
		echo $abstract"\\\\" >> $ofile
		echo "\\\\" >> $ofile
		echo nfiles $nfiles"\\\\">>$ofile
		echo "\\\\" >> $ofile
		echo nbib $nbib"\\\\">>$ofile
		echo "\\\\" >> $ofile
		echo "\clearpage"  >> $ofile
   fi
   
   #check primary header
   echo "\subsection{Primary header: required keywords}" >>  $ofile
   $mydfits -x -1 "$f" > $hdr
   nh=1 # define the header number
   naxis0=$(grep -c "NAXIS   =                    0" $hdr)
   #echo naxis0 $naxis0
   if [ $naxis0 -eq 1 ]; then
   if   [ "$tipo" == 'image'  ]; then
        true
        echo "ERROR NO data in the primary hdu !!!" 
   elif [ "$tipo" == 'iflux'  ]; then
        true
        echo "ERROR NO data in the primary hdu !!!" 
   fi
   fi
   if [ $naxis0 -ne 1 ]; then
   if   [ "$tipo" == 'imef'   ]; then
        true
        echo "ERROR data in the primary hdu !!!" 
   elif [ "$tipo" == 'mefi'   ]; then
        true
        echo "ERROR data in the primary hdu !!!" 
   elif [ "$tipo" == 'spec'   ]; then
        true
        echo "ERROR data in the primary hdu !!!" 
   elif [ "$tipo" == 'cube'   ]; then
        true
        echo "ERROR data in the primary hdu !!!" 
   elif [ "$tipo" == 'visi'   ]; then
        true
        echo "ERROR data in the primary hdu !!!" 
   elif [ "$tipo" == 'srctbl' ]; then
        true
        echo "ERROR data in the primary hdu !!!" 
   elif [ "$tipo" == 'mcat'   ]; then
        true
        echo "ERROR data in the primary hdu !!!" 
   elif [ "$tipo" == 'ctile'  ]; then
        true
        echo "ERROR data in the primary hdu !!!" 
   elif [ "$tipo" == 'catal'  ]; then
        true
        echo "ERROR data in the primary hdu !!!" 
   fi
   fi
   
   #exit 1
   #if   [ "$tipo" == 'image'  ]; then
        #true
   #elif [ "$tipo" == 'imef'   ]; then
        #true
   #elif [ "$tipo" == 'mefi'   ]; then
        #true
   #elif [ "$tipo" == 'iflux'  ]; then
        #true
   #elif [ "$tipo" == 'spec'   ]; then
        #true
   #elif [ "$tipo" == 'cube'   ]; then
        #true
   #elif [ "$tipo" == 'visi'   ]; then
        #true
   #elif [ "$tipo" == 'srctbl' ]; then
        #true
   #elif [ "$tipo" == 'mcat'   ]; then
        #true
   #elif [ "$tipo" == 'ctile'  ]; then
        #true
   #elif [ "$tipo" == 'catal'  ]; then
        #true
   #fi

	ckw

   # other extensions
   if [ $ne -gt 0 ]; then
      nh=2
      seq 1 $ne | while read enum
      do
         echo "\subsection{Extension "$enum": required keywords}" >>  $ofile
         $mydfits -x $enum "$f" > $hdr
         extype=$(grep XTENSION $hdr | cut -d= -f2 | cut -d/ -f1 | sed -e "s/'//g" -e 's/ //g')
         if   [ "$tipo" == 'image' ] && [ -n $extype ] ; then
              echo "EXTENSIONS NOT CHECKED" >>  $ofile
         elif [ "$tipo" == 'iflux' ] && [ -n $extype ] ; then
              echo "EXTENSIONS NOT CHECKED" >>  $ofile
         fi
         if [ $p3hdr = "true" ]; then
            echo >> $p3hdrfile
            echo extension $ne >> $p3hdrfile
            echo >> $p3hdrfile
         fi
         ckw
      done
   fi
   echo "\clearpage"  >> $ofile

   #print all headers
   if [ $allheaders == true ]; then 
      echo "\subsection{All headers}" >>  $ofile
      echo "\begin{verbatim}">> $ofile
      $mydfits -x 0 "$f" 2>&1 |sed -e 's/\#/*/g' -e 's/_/\_/g' | grep . |awk '{printf"%s \n",$0}'>>$ofile
      echo "\end{verbatim}"  >> $ofile
      echo "\clearpage"  >> $ofile
   fi
   echo "\cleardoublepage"  >> $ofile
   
done


echo "\end{document}"      >>$ofile
rm -rf $hdr $simb &>/dev/null 2>&1

if [ $latex_ = true ]; then

	#compile LaTeX
	pdflatex -synctex=1 -interaction=nonstopmode $ofile &>/dev/null 2>&1
	pdflatex -synctex=1 -interaction=nonstopmode $ofile &>/dev/null 2>&1
	pdflatex -synctex=1 -interaction=nonstopmode $ofile &>/dev/null 2>&1

	if [ -e $opdf ]; then
		pdfview_set=0
		for pdfview in evince okular xpdf gv
		do
			which $pdfview  &>/dev/null 2>&1
			ret_code=$?
			if (( $ret_code == 0 )) && (( $pdfview_set == 0 )); then
				pdfview_set=1
				#echo $edname $edname_set
				$pdfview $opdf &>/dev/null 2>&1&
				exit 1
			fi
		done
		#gnome-open $opdf&>/dev/null 2>&1
		echo 'report available in ' $ofile "and " $opdf
		echo
	else
		echo "error on compilation LaTeX report"
		echo "required LaTeX packages: graphicx,fancyvrb,fvextra,xcolor,hyperref"
		echo
		echo 'report available in ' $ofile
		echo
	fi
fi

echo 'done'

