#ifndef DEF_MY_HEADER_1
#define DEF_MY_HEADER_1

#include <cmath>
#include <iostream>
using namespace std;

/* 簡略化の記号 */


/* 専用共有ライブラリ用名前空間 */
namespace CommonLibs{

	/*
	//=======================================================
	// ■ CommonFuncs:共有用関数たち管理クラス
	//=======================================================
	// 共有的に使う関数を提供
	//=======================================================*/
	class CommonFuncs{
	public:
		static int round(double m);
		inline static double sigmoid(double x, double constK){return( 1.0/(1.0 + exp(-1.0*constK*x))  );};	/* シグモイド */
		static void setGauss(double* tg, double* wg, int numP );																	/* ガウス積分点[-1, 1] */
		static void setGaussHEX(double* tg_x, double* tg_y, double* tg_z, double* wg, int numP );									/* ガウス積分点(六面体用) */
		static void setGaussPRI(double* tg_x, double* tg_y, double* tg_z, double* wg_xy, double* wg_z, int numP_xy, int numP_z);	/* ガウス積分点(プリズム用) */
	};

}

#endif