#!/bin/sh

PRISM_JAVAMAXMEM=8g ../../../bin/prism team-form-offline-fc.prism team-form-fc.props > res-fc-offline.txt;
PRISM_JAVAMAXMEM=8g ../../../bin/prism team-form-offline-r.prism team-form-r.props > res-r-offline.txt;
PRISM_JAVAMAXMEM=8g ../../../bin/prism team-form-offline-s.prism team-form-s.props > res-s-offline.txt;
PRISM_JAVAMAXMEM=8g ../../../bin/prism team-form-offline-ia.prism team-form-ia.props > res-ia-offline.txt;

PRISM_JAVAMAXMEM=8g ../../../bin/prism team-form-online-fc.prism team-form-fc.props > res-fc-online.txt;
PRISM_JAVAMAXMEM=8g ../../../bin/prism team-form-online-r.prism team-form-r.props > res-r-online.txt;
PRISM_JAVAMAXMEM=8g ../../../bin/prism team-form-online-s.prism team-form-s.props > res-s-online.txt;
PRISM_JAVAMAXMEM=8g ../../../bin/prism team-form-online-ia.prism team-form-ia.props > res-ia-online.txt;