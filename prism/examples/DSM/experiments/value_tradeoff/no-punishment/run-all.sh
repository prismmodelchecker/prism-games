#!/bin/sh

../../../../../bin/prism mdsm2304.prism mdsm2304.props > mdsm2304.res
../../../../../bin/prism mdsm3304.prism mdsm3304.props > mdsm3304.res
../../../../../bin/prism mdsm4304.prism mdsm4304.props > mdsm4304.res
PRISM_JAVAMAXMEM=8g ../../../../../bin/prism mdsm5304.prism mdsm5304.props > mdsm5304.res
PRISM_JAVAMAXMEM=16g ../../../../../bin/prism mdsm6304.prism mdsm6304.props > mdsm6304.res
PRISM_JAVAMAXMEM=32g ../../../../../bin/prism mdsm7304.prism mdsm7304.props > mdsm7304.res