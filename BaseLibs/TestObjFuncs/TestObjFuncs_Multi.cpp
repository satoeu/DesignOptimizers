
#include "TestObjFuncs.hpp"

/*
//=======================================================
// ■ TestObjFuncs
//=======================================================
// ベンチマーク用関数の管理クラス
//=======================================================*/


/*//===========================================================
// ● 複数目的ベンチマーク関数値を返す 
//=========================================================== */
void TestObjFuncs::MultiObj(const int obj_num, double* objF, const string& func_name, const int s, const double *x){
	if(func_name == "ZDT1"){
		MultiObj(obj_num, objF, 0, s, x);
		return;
	}
	if(func_name == "ZDT2"){
		MultiObj(obj_num, objF, 1, s, x);
		return;
	}
	if(func_name == "ZDT3"){
		MultiObj(obj_num, objF, 2, s, x);
		return;
	}
	if(func_name == "ZDT4"){
		MultiObj(obj_num, objF, 3, s, x);
		return;
	}
	if(func_name == "DTLZ1"){
		MultiObj(obj_num, objF, 4, s, x);
		return;
	}
	if(func_name == "DTLZ2"){
		MultiObj(obj_num, objF, 5, s, x);
		return;
	}
	if(func_name == "DTLZ3"){
		MultiObj(obj_num, objF, 6, s, x);
		return;
	}
}

/*//===========================================================
// ● 複数目的ベンチマーク関数値を返す
//=========================================================== */
void TestObjFuncs::MultiObj(const int obj_num, double* objF, const int type, const int s, const double *x){
	/* ZDT1 */
	if(type == 0){
		objF[0] = x[0];
		double temp = 0;
		for(int i = 1 ; i < s ; i++){
			temp += x[i];
		}
		temp /= s-1.0;
		double g = 1.0 + 9.0 * temp;
		double h = 1.0 - sqrt(objF[0]/g);
		objF[1] = g*h;
	/* ZDT2 */
	}else if(type == 1){
		objF[0] =x[0];
		double temp = 0;
		for(int i = 1 ; i < s ; i++){
			temp += x[i];
		}
		temp /= s-1.0;
		double g = 1.0 + 9.0 * temp;
		double h = 1.0 - objF[0]*objF[0]/g/g;
		objF[1] = g*h;
	/* ZDT3 */
	}else if(type == 2){
		objF[0] =x[0];
		double temp = 0;
		for(int i = 1 ; i < s ; i++){
			temp += x[i];
		}
		temp /= s-1.0;
		double g = 1.0 + 9.0 * temp;
		double h = 1.0 - sqrt(objF[0]/g) - (objF[0]/g)*sin(10.0*3.141592*objF[0]);
		objF[1] = g*h;
	/* ZDT4 */
	}else if(type == 3){
		objF[0] = x[0];
		double temp = 0;
		for(int i = 1 ; i < s ; i++){
			temp += x[i]*x[i] - 10.0*cos(4.0*3.141592*x[0]);
		}
		double g = 1.0 + 10.0*(s-1.0) + temp;
		double h = 1.0 - sqrt(objF[0]/g);
		objF[1] = g*h;
	/* DTLZ1 */
	}else if(type == 4){
		double k = s - obj_num + 1;
		double g = 0.0;

		for (int i = s - k; i < s; i++){
			g += (x[i] - 0.5) * (x[i] - 0.5) - cos(20.0 * 3.141592 * (x[i] - 0.5));
		}
		g = 100.0 * (k + g);
		for (int i = 0; i < obj_num; i++){
			objF[i] = (1.0 + g) * 0.5;
		}

		for(int i = 0; i < obj_num; i++){
			for(int j = 0; j < obj_num - (i + 1); j++){
				objF[i] *= x[j];
			}
			if(i != 0){
				int aux = obj_num - (i + 1);
				objF[i] *= 1 - x[aux];
			}
		}
	/* DTLZ2 */
	}else if(type == 5){
		double k = s - obj_num + 1;
		double g = 0.0;

		for (int i = s - k; i < s; i++){
			g += (x[i] - 0.5) * (x[i] - 0.5);
		}
		for (int i = 0; i < obj_num; i++){
			objF[i] = 1.0 + g;
		}
		for (int i = 0; i < obj_num; i++) {
			for (int j = 0; j < obj_num - (i + 1); j++){
				objF[i] *= cos(x[j] * 0.5 * 3.141592);
			}
			if (i != 0) {
				int aux = obj_num - (i + 1);
				objF[i] *= sin(x[aux] * 0.5 * 3.141592);
			}
		}
	/* DTLZ3 */
	}else if(type == 6){
		double temp = 0;
		for(int i = obj_num-1 ; i < s ; i++){
			temp += ( (x[i]-0.5)*(x[i]-0.5) - cos(20.0*3.141592*(x[i]-0.5)) );
		}
		double g = 100.0*( s - obj_num + 1.0 + temp );
		double f1 = 1.0;
		for(int i = 0 ; i < obj_num-1 ; i++){
			f1 *= cos(x[i]*3.141592*0.5);
		}
		objF[0] = (1.0 + g) * f1;

		for(int mm = 1 ; mm < obj_num-1 ; mm++){
			double ff=1.0;
			for(int i = 0 ; i < obj_num-mm ; i++){
				ff *= cos(x[i]*3.141592*0.5);
			}
			objF[mm] = (1.0 + g) * ff * sin(3.141592*0.5*x[obj_num-mm]);
		}
		objF[obj_num-1] = (1.0 + g) * sin(3.141592*0.5*x[0]);
	}
}

