#!/bin/bash

project=$1 #project=PAS0471

grep -h "/fs/ess/$project" /users/reporting/storage/quota/*_quota.txt |
    sed -E 's/.*userid ([[:alnum:]]+).*used ([0-9]+).*/\1\t\2/' |
    sort -k2,2nr
    