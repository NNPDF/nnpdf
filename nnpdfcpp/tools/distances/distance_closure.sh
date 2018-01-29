#!/bin/bash

make clean
make


rm -rf *.res
touch davtl.res

# Level 0, all data, LevelZeroTL=(1 5 10 20 40 80)
#for fit in "130918-r1178-010-jr" "130918-r1178-008-jr" "130918-r1178-006-jr" "130918-r1178-001-jr"  "130918-r1178-003-jr" "130913-r1176-001-ag" 

# Level 2, all data, LevelTwoTL=(1 5 10 20 80)
#for fit in  "130918-r1178-011-jr" "130918-r1178-009-jr" "130918-r1178-007-jr" "130918-r1178-002-jr" "130913-r1176-002-ag"

# Level 0, 50% of the data, LevelZeroTL=(1 5 10 20 40)
#for fit in "131005-r1200-001-jr" "131005-r1200-002-jr" "131005-r1200-003-jr" "131005-r1200-004-jr" "131005-r1200-005-jr"

# Level 2, 50% of the data, LevelZeroTL=(1 5 10 20 40)
#for fit in "131005-r1200-006-jr" "131005-r1200-007-jr" "131005-r1200-008-jr" "131005-r1200-009-jr" "131005-r1200-010-jr"

# For compare PDF errors in level 0 and level 2 fits, 50% data
#for fit in "131005-r1200-001-jr" "131005-r1200-002-jr" "131005-r1200-003-jr" "131005-r1200-004-jr" "131005-r1200-005-jr" "131005-r1200-006-jr" "131005-r1200-007-jr" "131005-r1200-008-jr" "131005-r1200-009-jr" "131005-r1200-010-jr"

# Level 2, 25% of the data, LevelZeroTL=(1 5 10 20 40)
#for fit in "131007-r1200-007-jr" "131007-r1200-008-jr" "131007-r1200-009-jr" "131007-r1200-010-jr" "131007-r1200-011-jr"

# Level 2, 10% of the data, LevelZeroTL=(1 5 10 20 40)
#for fit in "131007-r1200-001-jr" "131007-r1200-003-jr" "131007-r1200-004-jr" "131007-r1200-005-jr" 

# Level0 , weight penalty fits
for fit in "131008-r1203-002-ag" "131008-r1203-003-ag" "131008-r1203-001-ag"

do

    rm -rf dav.res
    ./distance_closure $fit MSTW2008nlo68cl
    cat dav.res >> davtl.res

done

echo "davtl.res computed"

# Now for the plots
# Need to update the training lenght by hand
