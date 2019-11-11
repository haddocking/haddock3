#!/bin/csh
#
./protocols/cns1 <<_Eod_ |grep CNS-OK
set message=off echo=off end
display CNS-OK
stop
_Eod_
