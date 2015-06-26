#!/bin/sh
date >> /home/samio/log/rsyncMaia
/usr/bin/rsync -avzr --exclude '*RESTART*' --exclude '*core.*' --exclude 'LATEST' --exclude 'GEOMETRY' --exclude 'KPTS_GENERATION' -e "ssh -i /home/samio/.ssh/id_rsa_maia" /home/samio/bin /home/samio/src/quantumChemistry /home/samio/Works/PhD changk@131.152.22.89:./save/R700_backup >> /home/samio/log/rsyncMaia

#/usr/bin/rsync -avzr --exclude '*RESTART*' --exclude '*core.*' --exclude 'LATEST' --exclude 'GEOMETRY' --exclude 'KPTS_GENERATION' --delete -e "ssh -i /home/samio/.ssh/id_rsa_maia" /home/samio/bin changk@131.152.22.89:./save/R700_backup >> /home/samio/log/rsyncMaia
echo >> /home/samio/log/rsyncMaia
