#!/bin/sh

PRISM_JAVAMAXMEM=8g ../../../bin/prism team-form-offline-fc.prism team-form-fc.props > team-form-offline-fc.res;
PRISM_JAVAMAXMEM=8g ../../../bin/prism team-form-offline-r.prism team-form-r.props > team-form-offline-r.res;
PRISM_JAVAMAXMEM=8g ../../../bin/prism team-form-offline-s.prism team-form-s.props > team-form-offline-s.res;
PRISM_JAVAMAXMEM=8g ../../../bin/prism team-form-offline-ia.prism team-form-ia.props > team-form-offline-ia.res;

PRISM_JAVAMAXMEM=8g ../../../bin/prism team-form-online-fc.prism team-form-fc.props > team-form-online-fc.res;
PRISM_JAVAMAXMEM=8g ../../../bin/prism team-form-online-r.prism team-form-r.props > team-form-online-r.res;
PRISM_JAVAMAXMEM=8g ../../../bin/prism team-form-online-s.prism team-form-s.props > team-form-online-s.res;
PRISM_JAVAMAXMEM=8g ../../../bin/prism team-form-online-ia.prism team-form-ia.props > team-form-online-ia.res;