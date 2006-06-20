#/bin/bash

 
#   This program takes on stdin a file containing lines in the format
#
#   [NMPZDU](("0x"[0123456789ABCDEFabcdef]){8}" "){4}"#".* "\n"
#
#   e.g.
#
#   N 0xaffe1234 0x5678cafe 0x3FF921FB 0x54442D18 # RN(arccos(0xaffe12345678cafe))
#
#   Let be s the string obtained by concatenating the second and third word without the leading "0x".
#   No space character is taken into account.
#   Let be t be the string s in which all lower-case letters (a-f) have been replaced by upper-case 
#   letters.
#
#   The program prints out on stdout a file containing for each t the md5sum of the file containing the string t.
#
#


while read line 
do

xval=$(echo $line | sed -e "y/abcdef/ABCDEF/;s/#.*//;s/\([NUDMPZ]\) 0x\([012345678ABCDEF]\{8\}\) 0x\([012345678ABCDEF]\{8\}\) 0x\([012345678ABCDEF]\{8\}\) 0x\([012345678ABCDEF]\{8\}\)/\2\3/;s/[^0123456789ABCDEF]//g;/^$/d")

if test -n "$xval"; then
echo $xval | md5sum | sed -e "s/ -//"
fi

done