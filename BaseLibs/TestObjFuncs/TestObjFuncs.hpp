#ifndef DEF_TESTOBJ_DEF00
#define DEF_TESTOBJ_DEF00

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
using namespace std;

#include <BaseDefines.hpp>

/*
//=======================================================
// ■ TestObjFuncs
//=======================================================
// ベンチマーク用関数の管理クラス
//=======================================================*/
class TestObjFuncs{
public:
	static double SingleObj(const int type, const int s, const double *x);												/* 単目的ベンチマーク関数値を返す */
	static double SingleObj(const string& func_name, const int s, const double *x);										/* 単目的ベンチマーク関数値を返す */
	static void MultiObj(const int obj_num, double* objF, const int type, const int s, const double *x);				/* 複数目的ベンチマーク関数値を返す */
	static void MultiObj(const int obj_num, double* objF, const string& func_name, const int s, const double *x);		/* 複数目的ベンチマーク関数値を返す */
																												/**/
	static double Sphere(const int s, const double *x);					/* Sphere関数 */
	static double Ellipsoid(const int s, const double *x);				/* Ellipsoid関数 */
	static double Tablet(const int s, const double *x);					/* Tablet関数 */
	static double RosenBrock_c(const int s, const double *x);			/* RosenBrock_c関数 */
	static double Bohachevsky(const int s, const double *x);			/* Bohachevsky関数 */
	static double Schaffer(const int s, const double *x);				/* Schaffer関数 */
	static double Rastrigin(const int s, const double *x);				/* Rastrigin関数 */
	/**/
};


#endif
