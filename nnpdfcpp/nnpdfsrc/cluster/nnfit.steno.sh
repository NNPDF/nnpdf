#!/bin/bash

for I in {1..110}
do
   export REP=$I
   llsubmit nnfit.steno.ll.cmd
done
