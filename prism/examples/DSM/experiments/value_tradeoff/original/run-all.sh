#!/bin/sh

../../../../../bin/prism mdsm2324.prism mdsm2324.props > mdsm2324.res
../../../../../bin/prism mdsm3334.prism mdsm3334.props > mdsm3334.res
../../../../../bin/prism mdsm4344.prism mdsm4344.props > mdsm4344.res
PRISM_JAVAMAXMEM=8g ../../../../../bin/prism mdsm5354.prism mdsm5354.props > mdsm5354.res
PRISM_JAVAMAXMEM=16g ../../../../../bin/prism mdsm6364.prism mdsm6364.props > mdsm6364.res
PRISM_JAVAMAXMEM=32g ../../../../../bin/prism mdsm7374.prism mdsm7374.props > mdsm7374.res