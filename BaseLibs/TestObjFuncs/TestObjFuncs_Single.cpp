
#include "TestObjFuncs.hpp"

/*
//=======================================================
// ■ TestObjFuncs
//=======================================================
// ベンチマーク用関数の管理クラス
//=======================================================*/

/*//===========================================================
// ● 単目的ベンチマーク関数値を返す
//=========================================================== */
double TestObjFuncs::SingleObj(const int type, const int s, const double *x){
	double tmp;
	switch(type){
	case 0: default:
		tmp = TestObjFuncs::Sphere(s, x);
		break;
	case 1:
		tmp = TestObjFuncs::Ellipsoid(s, x);
		break;
	case 2:
		tmp = TestObjFuncs::Tablet(s, x);
		break;
	case 3:
		tmp = TestObjFuncs::RosenBrock_c(s, x);
		break;
	case 4:
		tmp = TestObjFuncs::Bohachevsky(s, x);
		break;
	case 5:
		tmp = TestObjFuncs::Schaffer(s, x);
		break;
	case 6:
		tmp = TestObjFuncs::Rastrigin(s, x);
		break;
	}
	return tmp;
}

/*//===========================================================
// ● 単目的ベンチマーク関数値を返す
//=========================================================== */
double TestObjFuncs::SingleObj(const string& func_name, const int s, const double *x){
	double tmp=0;
	if(func_name == "Sphere"){
		tmp = TestObjFuncs::Sphere(s, x);
		return tmp;
	}
	if(func_name == "Ellipsoid"){
		tmp = TestObjFuncs::Ellipsoid(s, x);
		return tmp;
	}
	if(func_name == "Tablet"){
		tmp = TestObjFuncs::Tablet(s, x);
		return tmp;
	}
	if(func_name == "RosenBrock_c"){
		tmp = TestObjFuncs::RosenBrock_c(s, x);
		return tmp;
	}
	if(func_name == "Bohachevsky"){
		tmp = TestObjFuncs::Bohachevsky(s, x);
		return tmp;
	}
	if(func_name == "Schaffer"){
		tmp = TestObjFuncs::Schaffer(s, x);
		return tmp;
	}
	if(func_name == "Rastrigin"){
		tmp = TestObjFuncs::Rastrigin(s, x);
		return tmp;
	}
	cout <<"Function name erorr !" << endl;
	exit(1);
	return 1.0;
}

/*//===========================================================
// ● Sphere
//=========================================================== */
double TestObjFuncs::Sphere(const int s, const double *x){
	double temp = 0;
	for(int i = 0 ; i < s ; i++){
		temp += x[i]*x[i];
	}
	return temp;
}
/*//===========================================================
// ● Ellipsoid
//=========================================================== */
double TestObjFuncs::Ellipsoid(const int s, const double *x){	
	double temp = 0;
	double A;
	for(int i = 0 ; i < s ; i++){
		A = pow( 1000.0, (i-1.0)/(s-1.0) );
		temp += A*A*x[i]*x[i];
	}
	return temp;
}
/*//===========================================================
// ● Tablet
//=========================================================== */
double TestObjFuncs::Tablet(const int s, const double *x){	
	double temp = 0;
	const int k = s / 4;
	for(int i = 0 ; i < k ; i++){
		temp += x[i]*x[i];
	}
	for(int i = k ; i < s ; i++){
		temp += 100.0*100.0*x[i]*x[i];
	}
	return temp;
}
/*//===========================================================
// ● RosenBrock_c
//=========================================================== */
double TestObjFuncs::RosenBrock_c(const int s, const double *x){	
	double temp = 0;
	double a1,a2;
	for(int i = 0 ; i < s-1 ; i++){
		a1 = 100.0*(x[i+1]-x[i]*x[i])*(x[i+1]-x[i]*x[i]);
		a2 = (1.0 - x[i])*(1.0 - x[i]);
		temp += a1 + a2;
	}
	return temp;
}
/*//===========================================================
// ● Bohachevsky
//=========================================================== */
double TestObjFuncs::Bohachevsky(const int s, const double *x){	
	double temp = 0;
	double a1,a2,a3,a4;
	for(int i = 0 ; i < s-1 ; i++){
		a1 = x[i]*x[i];
		a2 = 2.0 * x[i+1]*x[i+1];
		a3 = 0.3 * cos(3.0*3.1415*x[i]);
		a4 = 0.4 * cos(4.0*3.1415*x[i+1]);
		temp += a1 + a2 - a3 - a4 + 0.7;
	}
	return temp;
}
/*//===========================================================
// ● Schaffer
//=========================================================== */
double TestObjFuncs::Schaffer(const int s, const double *x){	
	double temp = 0;
	double a1,a2,a3,a4;
	for(int i = 0 ; i < s-1 ; i++){
		a1 = x[i]*x[i] + x[i+1]*x[i+1];
		a2 = pow(a1, 0.25);
		a3 = pow(a1, 0.1);
		a4 = sin(50.0*a3+1.0) * sin(50.0*a3+1.0);
		temp += a2 * a4;
	}
	return temp;
}
/*//===========================================================
// ● Rastrigin
//=========================================================== */
double TestObjFuncs::Rastrigin(const int s, const double *x){	
	double temp = 10.0*s;//200;
	double a1,a2;
	for(int i = 0 ; i < s ; i++){
		a1 = (x[i]-1.0) * (x[i]-1.0);
		a2 = 10.0 * cos(2.0*3.1415 * (x[i]-1.0));
		temp += a1 - a2;
	}
	return temp;
}
