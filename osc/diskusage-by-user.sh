#!/bin/bash

while read -r user; do
  echo -en "$user\t"
  find . -user "$user" -type f -print0 |
    du --files0-from=- --total -sh |
    tail -n 1 |
    awk '{print $1}'
done < <(find . -printf "%u\n"| sort -u)

# OR, check:
# grep -h "/fs/ess/PAS2880" /users/reporting/storage/quota/*_quota.txt