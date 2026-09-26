//==============================================================================
//	JNI wrapper for SoPlex (floating point), part of PRISM. GPL v2 or later.
//	Authors: Gabriel Santos (University of Oxford)
//==============================================================================

#include <jni.h>
#include <cmath>
#include <vector>
#include "soplex.h"

using namespace soplex;

static inline SoPlex *S(jlong p) { return reinterpret_cast<SoPlex *>(p); }

extern "C" {

JNIEXPORT jlong JNICALL Java_soplex_SoPlex_create(JNIEnv *, jclass)
{
	try {
		SoPlex *s = new SoPlex();
		s->setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);
		// defaults tuned for many tiny LPs (matrix games at each state): presolving and timing cost more than they save.
		// Note: with TIMER_OFF a time limit has no effect; setTimeLimit() switches the timer back on.
		s->setIntParam(SoPlex::SIMPLIFIER, SoPlex::SIMPLIFIER_OFF);
		s->setIntParam(SoPlex::TIMER, SoPlex::TIMER_OFF);
		return reinterpret_cast<jlong>(s);
	} catch (...) {
		return 0;
	}
}

JNIEXPORT void JNICALL Java_soplex_SoPlex_free(JNIEnv *, jclass, jlong p) { delete S(p); }

JNIEXPORT void JNICALL Java_soplex_SoPlex_clear(JNIEnv *, jclass, jlong p) { S(p)->clearLPReal(); }

JNIEXPORT void JNICALL Java_soplex_SoPlex_setObjSense(JNIEnv *, jclass, jlong p, jboolean max)
{
	S(p)->setIntParam(SoPlex::OBJSENSE, max ? SoPlex::OBJSENSE_MAXIMIZE : SoPlex::OBJSENSE_MINIMIZE);
}

JNIEXPORT void JNICALL Java_soplex_SoPlex_setVerbosity(JNIEnv *, jclass, jlong p, jint level) { S(p)->setIntParam(SoPlex::VERBOSITY, level); }

JNIEXPORT void JNICALL Java_soplex_SoPlex_setTimeLimit(JNIEnv *, jclass, jlong p, jdouble sec)
{
	S(p)->setIntParam(SoPlex::TIMER, SoPlex::TIMER_CPU);
	S(p)->setRealParam(SoPlex::TIMELIMIT, sec);
}

JNIEXPORT void JNICALL Java_soplex_SoPlex_setScaler(JNIEnv *, jclass, jlong p, jint scaler) { S(p)->setIntParam(SoPlex::SCALER, scaler); }

JNIEXPORT void JNICALL Java_soplex_SoPlex_setIntParam(JNIEnv *, jclass, jlong p, jint code, jint value)
{
	S(p)->setIntParam((SoPlex::IntParam) code, value);
}

JNIEXPORT void JNICALL Java_soplex_SoPlex_setRealParam(JNIEnv *, jclass, jlong p, jint code, jdouble value)
{
	S(p)->setRealParam((SoPlex::RealParam) code, value);
}

JNIEXPORT void JNICALL Java_soplex_SoPlex_setBoolParam(JNIEnv *, jclass, jlong p, jint code, jboolean value)
{
	S(p)->setBoolParam((SoPlex::BoolParam) code, value);
}

JNIEXPORT jint JNICALL Java_soplex_SoPlex_addCol(JNIEnv *, jclass, jlong p, jdouble obj, jdouble lb, jdouble ub)
{
	DSVector empty(0);
	S(p)->addColReal(LPCol(obj, empty, ub, lb));
	return S(p)->numColsReal() - 1;
}

JNIEXPORT jint JNICALL Java_soplex_SoPlex_addRow(JNIEnv *env, jclass, jlong p, jintArray idx, jdoubleArray vals, jdouble lhs, jdouble rhs)
{
	jsize k = env->GetArrayLength(idx);
	jint *ix = env->GetIntArrayElements(idx, nullptr);
	jdouble *vx = env->GetDoubleArrayElements(vals, nullptr);
	DSVector row(k);
	for (jsize i = 0; i < k; i++)
		if (vx[i] != 0.0)
			row.add(ix[i], vx[i]);
	env->ReleaseIntArrayElements(idx, ix, JNI_ABORT);
	env->ReleaseDoubleArrayElements(vals, vx, JNI_ABORT);
	S(p)->addRowReal(LPRow(lhs, row, rhs));
	return S(p)->numRowsReal() - 1;
}

/* Append columns in one call (all with empty column vectors; rows added afterwards). Returns the index of the first new column. */
JNIEXPORT jint JNICALL Java_soplex_SoPlex_addCols(JNIEnv *env, jclass, jlong p, jdoubleArray objj, jdoubleArray lbj, jdoubleArray ubj)
{
	jsize k = env->GetArrayLength(objj);
	jdouble *obj = env->GetDoubleArrayElements(objj, nullptr);
	jdouble *lb = env->GetDoubleArrayElements(lbj, nullptr);
	jdouble *ub = env->GetDoubleArrayElements(ubj, nullptr);
	int first = S(p)->numColsReal();
	DSVector empty(0);
	LPColSetReal cols(k);
	for (jsize i = 0; i < k; i++)
		cols.add(LPCol(obj[i], empty, ub[i], lb[i]));
	env->ReleaseDoubleArrayElements(objj, obj, JNI_ABORT);
	env->ReleaseDoubleArrayElements(lbj, lb, JNI_ABORT);
	env->ReleaseDoubleArrayElements(ubj, ub, JNI_ABORT);
	S(p)->addColsReal(cols);
	return first;
}

/* Replace the whole objective vector (length = number of columns); the current basis is kept (warm start). */
JNIEXPORT void JNICALL Java_soplex_SoPlex_changeObj(JNIEnv *env, jclass, jlong p, jdoubleArray objj)
{
	jsize n = env->GetArrayLength(objj);
	jdouble *obj = env->GetDoubleArrayElements(objj, nullptr);
	VectorBase<Real> v(n);
	for (jsize i = 0; i < n; i++)
		v[i] = obj[i];
	env->ReleaseDoubleArrayElements(objj, obj, JNI_ABORT);
	S(p)->changeObjReal(v);
}

JNIEXPORT void JNICALL Java_soplex_SoPlex_changeRowRange(JNIEnv *, jclass, jlong p, jint i, jdouble lhs, jdouble rhs)
{
	S(p)->changeRangeReal(i, lhs, rhs);
}

JNIEXPORT void JNICALL Java_soplex_SoPlex_changeColBounds(JNIEnv *, jclass, jlong p, jint i, jdouble lb, jdouble ub)
{
	S(p)->changeBoundsReal(i, lb, ub);
}

JNIEXPORT jint JNICALL Java_soplex_SoPlex_optimize(JNIEnv *, jclass, jlong p) { return (jint) S(p)->optimize(); }

JNIEXPORT jint JNICALL Java_soplex_SoPlex_getStatus(JNIEnv *, jclass, jlong p) { return (jint) S(p)->status(); }

JNIEXPORT jdouble JNICALL Java_soplex_SoPlex_getObjValue(JNIEnv *, jclass, jlong p) { return S(p)->objValueReal(); }

JNIEXPORT jint JNICALL Java_soplex_SoPlex_numRows(JNIEnv *, jclass, jlong p) { return S(p)->numRowsReal(); }

JNIEXPORT jint JNICALL Java_soplex_SoPlex_numCols(JNIEnv *, jclass, jlong p) { return S(p)->numColsReal(); }

JNIEXPORT void JNICALL Java_soplex_SoPlex_getPrimal(JNIEnv *env, jclass, jlong p, jdoubleArray x)
{
	jsize n = env->GetArrayLength(x);
	std::vector<double> buf(n);
	S(p)->getPrimalReal(buf.data(), n);
	env->SetDoubleArrayRegion(x, 0, n, buf.data());
}

JNIEXPORT void JNICALL Java_soplex_SoPlex_getDual(JNIEnv *env, jclass, jlong p, jdoubleArray y)
{
	jsize n = env->GetArrayLength(y);
	std::vector<double> buf(n);
	S(p)->getDualReal(buf.data(), n);
	env->SetDoubleArrayRegion(y, 0, n, buf.data());
}

JNIEXPORT void JNICALL Java_soplex_SoPlex_getRedCost(JNIEnv *env, jclass, jlong p, jdoubleArray r)
{
	jsize n = env->GetArrayLength(r);
	std::vector<double> buf(n);
	S(p)->getRedCostReal(buf.data(), n);
	env->SetDoubleArrayRegion(r, 0, n, buf.data());
}

JNIEXPORT jdouble JNICALL Java_soplex_SoPlex_getSolvingTime(JNIEnv *, jclass, jlong p) { return S(p)->solveTime(); }

/*
 * Matrix game in one call. Row player (min=false): columns x_0..x_{m-1} in [0,1], v in [vlb,vub];
 * rows: sum_i A[i][j] x_i - v >= 0 for each j; sum_i x_i = 1; maximise v.
 * Column player (min=true): columns y_0..y_{n-1}, w; rows: sum_j A[i][j] y_j - w <= 0 for each i; sum_j y_j = 1; minimise w.
 */
JNIEXPORT jdouble JNICALL Java_soplex_SoPlex_matrixGame(JNIEnv *env, jclass, jlong p, jdoubleArray Aj, jint m, jint n, jboolean min,
		jdouble vlb, jdouble vub, jdoubleArray stratOut)
{
	SoPlex *s = S(p);
	s->clearLPReal();
	s->setIntParam(SoPlex::OBJSENSE, min ? SoPlex::OBJSENSE_MINIMIZE : SoPlex::OBJSENSE_MAXIMIZE);
	jdouble *A = env->GetDoubleArrayElements(Aj, nullptr);
	int k = min ? n : m;    // strategy variables
	int h = min ? m : n;    // constraints (opponent's pure actions)
	DSVector empty(0);
	LPColSetReal cols(k + 1);
	for (int i = 0; i < k; i++)
		cols.add(LPCol(0.0, empty, 1.0, 0.0));
	cols.add(LPCol(1.0, empty, vub, vlb));
	s->addColsReal(cols);
	LPRowSetReal rows(h + 1);
	DSVector row(k + 1);
	for (int o = 0; o < h; o++) {
		row.clear();
		for (int i = 0; i < k; i++) {
			double a = min ? A[o * n + i] : A[i * n + o];
			if (a != 0.0)
				row.add(i, a);
		}
		row.add(k, -1.0);
		if (min)
			rows.add(LPRow(-infinity, row, 0.0));
		else
			rows.add(LPRow(0.0, row, infinity));
	}
	row.clear();
	for (int i = 0; i < k; i++)
		row.add(i, 1.0);
	rows.add(LPRow(1.0, row, 1.0));
	s->addRowsReal(rows);
	env->ReleaseDoubleArrayElements(Aj, A, JNI_ABORT);
	SPxSolver::Status st = s->optimize();
	if (st != SPxSolver::OPTIMAL)
		return NAN;
	if (stratOut != nullptr) {
		std::vector<double> x(k + 1);
		s->getPrimalReal(x.data(), k + 1);
		env->SetDoubleArrayRegion(stratOut, 0, k, x.data());
	}
	return s->objValueReal();
}

} // extern "C"
