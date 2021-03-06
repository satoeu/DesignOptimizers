#include "CommonFuncs.hpp"


namespace CommonLibs{


/*
//--------------------------------------------------------------------------
// ● 四捨五入
//--------------------------------------------------------------------------*/
int CommonFuncs::round(double m){
	double x;
	if(m >= 0){
		x = m + 0.5;
	}else{
		x = m - 0.5;
	}
	return ( (int)(x) );
}




/* ********************************************************************
// ガウスの数値積分用の定数データをセット
********************************************************************* */
void CommonFuncs::setGauss(double* tg, double* wg, int numP ){

/* input:
//    numP:積分点数
//
// output;
//    tg  :積分点の局所座標
//    wg  :積分の重み
//
//  エラー時には例外ErrorFEMを投げる。戻り値はなし
//  tg wgはともにtg[numP], wg[numP]だけの領域を
//  あらかじめ確保しておく必要がある。
*/
	switch(numP){
		case 2:
			tg[0] = 5.77350269189626e-1;
			tg[1] = -tg[0];
			
			wg[0] = 1.0;
			wg[1] = wg[0];
			break;
		case 3:
			tg[0] = 0.0;
			tg[1] = 7.74596669241483e-1;
			tg[2] = -tg[1];
				
			wg[0] = 8.88888888888889e-1;
			wg[1] = 5.55555555555556e-1;
			wg[2] = wg[1];
			break;
		case 4:
			tg[0] = 8.61136311594053e-1;
			tg[1] = -tg[0];
			tg[2] = 3.39981043584856e-1;
			tg[3] = -tg[2];

			wg[0] = 3.47854845137456e-1;
			wg[1] = wg[0];
			wg[2] = 6.52145154862546e-1;
			wg[3] = wg[2];
			break;
		case 6:
			tg[0] = 9.32469514203152e-1;
			tg[1] = -tg[0];
			tg[2] = 6.61209386466265e-1;
			tg[3] = -tg[2];
			tg[4] = 2.38619186083197e-1;
			tg[5] = -tg[4];

			wg[0] = 1.71324492379170e-1;
			wg[1] = wg[0];
			wg[2] = 3.60761573048139e-1;
			wg[3] = wg[2];
			wg[4] = 4.67913934572691e-1;
			wg[5] = wg[4];
			break;
		case 8:
			tg[0] = 9.60289856497536e-1;
			tg[1] = -tg[0];
			tg[2] = 7.96666477413627e-1;
			tg[3] = -tg[2];
			tg[4] = 5.25532409916329e-1;
			tg[5] = -tg[4];
			tg[6] = 1.83434642495650e-1;
			tg[7] = -tg[6];

			wg[0] = 1.01228536290376e-1;
			wg[1] = wg[0];
			wg[2] = 2.22381034453374e-1;
			wg[3] = wg[2];
			wg[4] = 3.13706645877887e-1;
			wg[5] = wg[4];
			wg[6] = 3.62683783378362e-1;
			wg[7] = wg[6];
			break;
		case 10:
			tg[0] = 9.73906528517172e-1;
			tg[1] = -tg[0];
			tg[2] = 8.65063366688985e-1;
			tg[3] = -tg[2];
			tg[4] = 6.79409568299024e-1;
			tg[5] = -tg[4];
			tg[6] = 4.33395394129247e-1;
			tg[7] = -tg[6];
			tg[8] = 1.48874338981631e-1;
			tg[9] = -tg[8];

			wg[0] = 6.6671344308688e-2;
			wg[1] = wg[0];
			wg[2] = 1.49451349150581e-1;
			wg[3] = wg[2];
			wg[4] = 2.19086362515982e-1;
			wg[5] = wg[4];
			wg[6] = 2.69266719309996e-1;
			wg[7] = wg[6];
			wg[8] = 2.95524224714753e-1;
			wg[9] = wg[8];
			break;
		case 12:
			tg[0] = 9.81560634246719e-1;
			tg[1] = -tg[0];
			tg[2] = 9.04117256370475e-1;
			tg[3] = -tg[2];
			tg[4] = 7.69902674194305e-1;
			tg[5] = -tg[4];
			tg[6] = 5.87317954286617e-1;
			tg[7] = -tg[6];
			tg[8] = 3.67831498998180e-1;
			tg[9] = -tg[8];
			tg[10]= 1.25233408511469e-1;
			tg[11]= -tg[10];

			wg[0] = 4.7175336386512e-2;
			wg[1] = wg[0];
			wg[2] = 1.06939325995318e-1;
			wg[3] = wg[2];
			wg[4] = 1.60078328543346e-1;
			wg[5] = wg[4];
			wg[6] = 2.03167426723066e-1;
			wg[7] = wg[6];
			wg[8] = 2.33492536538355e-1;
			wg[9] = wg[8];
			wg[10]= 2.49147045813403e-1;
			wg[11]= wg[10];
			break;
			
		default:
			cout << "numP is wrong number (in setGauss)" << endl;
	}
	return;
}

/* ********************************************************************
// ガウスの数値積分用の定数データをセット（六面体・三次元点数ver）
********************************************************************* */
void CommonFuncs::setGaussHEX(double* tg_x, double* tg_y, double* tg_z, double* wg, int numP ){

/* input:
//    numP:積分点数
//
// output;
//    tg  :積分点の局所座標
//    wg  :積分の重み
//
//  エラー時には例外ErrorFEMを投げる。戻り値はなし
//  tg wgはともにtg[numP], wg[numP]だけの領域を
//  あらかじめ確保しておく必要がある。
*/
	switch(numP){
		case 4:
			tg_x[0] = 0.0;					tg_y[0] = sqrt(2.0/3.0);		tg_z[0] = -1.0/sqrt(3.0);
			tg_x[1] = 0.0;					tg_y[1] = -1.0*sqrt(2.0/3.0);	tg_z[1] = -1.0/sqrt(3.0);
			tg_x[2] = sqrt(2.0/3.0);		tg_y[2] = 0.0;					tg_z[2] = 1.0/sqrt(3.0);
			tg_x[3] = -1.0*sqrt(2.0/3.0);	tg_y[3] = 0.0;					tg_z[3] = 1.0/sqrt(3.0);
							
			wg[0] = 2.0; wg[1] = 2.0;
			wg[2] = 2.0; wg[3] = 2.0;
			break;
		case 6:
			tg_x[0] = 1.0/sqrt(6.0);		tg_y[0] = 1.0/sqrt(2.0);	tg_z[0] = -1.0/sqrt(3.0);
			tg_x[1] = 1.0/sqrt(6.0);		tg_y[1] = -1.0/sqrt(2.0);	tg_z[1] = -1.0/sqrt(3.0);
			tg_x[2] = -1.0/sqrt(6.0);		tg_y[2] = 1.0/sqrt(2.0);	tg_z[2] = 1.0/sqrt(3.0);
			tg_x[3] = -1.0/sqrt(6.0);		tg_y[3] = -1.0/sqrt(2.0);	tg_z[3] = 1.0/sqrt(3.0);
			tg_x[4] = -1.0*sqrt(2.0/3.0);	tg_y[4] = 0.0;				tg_z[4] = -1.0/sqrt(3.0);
			tg_x[5] = sqrt(2.0/3.0);		tg_y[5] = 0.0;				tg_z[5] = 1.0/sqrt(3.0);

			wg[0] = 4.0/3.0; wg[1] = 4.0/3.0;
			wg[2] = 4.0/3.0; wg[3] = 4.0/3.0;
			wg[4] = 4.0/3.0; wg[5] = 4.0/3.0;
			break;
		default:
			cout << "numP is wrong number (in setGauss)" << endl;
	}
	return;
}

/* ********************************************************************
// ガウスの数値積分用の定数データをセット（プリズム・三次元点数ver）
********************************************************************* */
void CommonFuncs::setGaussPRI(double* tg_x, double* tg_y, double* tg_z, double* wg_xy, double* wg_z, int numP_xy, int numP_z){

/* input:
//    numP:積分点数
//
// output;
//    tg  :積分点の局所座標
//    wg  :積分の重み
//
//  エラー時には例外ErrorFEMを投げる。戻り値はなし
//  tg wgはともにtg[numP], wg[numP]だけの領域を
//  あらかじめ確保しておく必要がある。
*/
	switch(numP_xy){
		case 1: 
			tg_x[0] = 1.0/3.0;				tg_y[0] = 1.0/3.0;
							
			wg_xy[0] = 1.0;
			break;
		case 3: 
			tg_x[0] = 0.5;					tg_y[0] = 0.5;
			tg_x[1] = 0.0;					tg_y[1] = 0.5;
			tg_x[2] = 0.5;					tg_y[2] = 0.0;
							
			wg_xy[0] = 1.0/3.0; wg_xy[1] = 1.0/3.0; wg_xy[2] = 1.0/3.0;
			break;
		default:
			cout << "numPxy is wrong number (in setGauss)" << endl;
	}
	switch(numP_z){
		case 1: 
			tg_z[0] = 0.0;							
			wg_z[0] = 2.0;
			break;
		case 2: 
			tg_z[0] = 1.0/sqrt(3.0); tg_z[1] = -1.0/sqrt(3.0);
			wg_z[0] = 1.0; wg_z[1] = 1.0;
			break;
		case 3: 
			tg_z[0] = 0.0; tg_z[1] = -1.0*sqrt(3.0/5.0);  tg_z[2] = 1.0*sqrt(3.0/5.0);
							
			wg_z[0] = 8.0/9.0; wg_z[1] = 5.0/9.0; wg_z[2] = 5.0/9.0;
			break;
		default:
			cout << "numPz is wrong number (in setGauss)" << endl;
	}
	return;
}

}
