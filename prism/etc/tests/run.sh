#!/bin/sh

set -e # Abort if one of the commands fail
set -x # Print commands as they are executed

# Tests to check key libraries/components are installed

# Can be run manually from this directory, with "prism-test ." (or "prism-auto -tm .")

# This script is intended to be run from from ../..
# (the below is generated from ../.. by: "prism-test -e etc/tests -p bin/prism")

bin/prism etc/tests/dtmc_pctl.prism etc/tests/dtmc_pctl.prism.props -ex -test
bin/prism etc/tests/dtmc_pctl.prism etc/tests/dtmc_pctl.prism.props -h -test
bin/prism etc/tests/test_ppl_smg.prism etc/tests/test_ppl_smg.prism.props -test
bin/prism etc/tests/test_yices_csg.prism etc/tests/test_yices_csg.prism.props -smtsolver yices -test
bin/prism etc/tests/test_z3_csg.prism etc/tests/test_z3_csg.prism.props -smtsolver z3 -test
