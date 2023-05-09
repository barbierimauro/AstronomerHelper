mpfile="$(mktemp -p /dev/shm/)"
dfits -x 0 $1>$tmpfile

TTYPE=$(grep -c ^TTYPE $tmpfile)
TFORM=$(grep -c ^TFORM $tmpfile)
TCOMM=$(grep -c ^TCOMM $tmpfile)
TUCD=$( grep -c ^TUCD  $tmpfile)
TUNIT=$(grep -c ^TUNIT $tmpfile)

echo 'TTYPE' $TTYPE
echo 'TFORM' $TFORM
echo 'TCOMM' $TCOMM
echo 'TUCD ' $TUCD
echo 'TUNIT' $TUNIT


grep ^TTYPE $tmpfile | cut -d= -f1 |cut -c6-10| sort> ttype.list
grep ^TFORM $tmpfile | cut -d= -f1 |cut -c6-10| sort> tform.list
grep ^TCOMM $tmpfile | cut -d= -f1 |cut -c6-10| sort> tcomm.list
grep ^TUCD  $tmpfile | cut -d= -f1 |cut -c6-10| sort> tucd.list
grep ^TUNIT $tmpfile | cut -d= -f1 |cut -c6-10| sort> tunit.list


#comm --total -23 ttype.list tunit.list
#da i ttype che sono unici e non presenti in tunit

#comm -23 ttype.list tunit.list| awk '{printf"exthdr\['\''TUNIT%i'\''] = '\''        '\''\n",$1}'
#aggiunge le righe mancanti con python
