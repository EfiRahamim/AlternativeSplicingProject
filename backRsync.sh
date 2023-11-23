#!/bin/bash

source_dir=/private10/Projects/Efi/AML_thirdbatch/
destiny_dir=/private8/Projects/Efi/AML_backup/ThirdBatch_SFmutationsOnly/

#mkdir $destiny_dir

while [ 1 ]
do
    rsync -avX --timeout=60 --partial $source_dir $destiny_dir
    if [ "$?" = "0" ] ; then
        echo "rsync completed normally - ${source_dir}"
        break
    else
        echo "Rsync failure. Backing off and retrying..."
        sleep 180
    fi
done

echo "All Done!"