#include "Mt.hpp"


namespace CommonLibs{


/*
//=======================================================
// ■ メルセンヌ乱数使用したメソッド管理staticクラス
//=======================================================
//  メルセンヌ乱数を用いた処理を提供するクラス
//=======================================================*/

/* 初期化 */
#ifdef MT_TYPE_64
std::mt19937_64 Mt::mt_r;
#else
std::mt19937 Mt::mt_r;
#endif

/*//=======================================================
// ● シード初期化
//=======================================================*/
void Mt::init_rand() { 
#ifdef MT_TYPE_64
	std::random_device seed_gen;
	std::mt19937_64 temp(seed_gen());
	mt_r = std::move(temp);
#else
	std::random_device seed_gen;
	std::mt19937 temp(seed_gen());
	mt_r = std::move(temp);
#endif
} 
/*//=======================================================
// ● シード初期化
//=======================================================*/
void Mt::init_rand(int xx) { 
#ifdef MT_TYPE_64
	std::mt19937_64 temp(xx);
	mt_r = std::move(temp);
#else
	std::mt19937 temp(xx);
	mt_r = std::move(temp);
#endif
} 


/*//=======================================================
// ● 整数乱数デフォルト生成（0～最大まで）
//=======================================================*/
int Mt::genrand_int32(){
#ifdef MT_TYPE_64
	std::uint64_t result = mt_r();
#else
	std::uint32_t result = mt_r();
#endif
	return( result );
}

/*//=======================================================
// ● ｍ未満の整数を返す
//=======================================================*/
int Mt::mrand(int m){
	std::uniform_int_distribution<> dist(0, m-1);
	const int result = dist(mt_r);
	return result;
}
/*//=======================================================
// ● a<=x<=bの範囲の一様整数乱数を返す
//=======================================================*/
int Mt::mrand(int a, int b){
	std::uniform_int_distribution<> dist(a, b);
	const int result = dist(mt_r);
	return result;
}


/*//=======================================================
// ● a<=x<=bの範囲の一様乱数を返す
//=======================================================*/
double Mt::unif_rand(double a, double b){
	std::uniform_real_distribution<> dist(a, b);
	const double result = dist(mt_r);
	return result;
}


/*//=======================================================
// ● 正規分布生成(limit_sig=true: 絶対値４sigmaまでの範囲でしか出力しない) 
//=======================================================*/
double Mt::normal_rand(double u, double sigma, bool limit_sig){
	std::normal_distribution<> dist(u, sigma);
	if(limit_sig){
		while(true){
			const double result = dist(mt_r);
			if(u-4.0*sigma <= result && result <= u+4.0*sigma){
				return result;
			}
		}
	}else{
		const double result = dist(mt_r);
		return result;
	}
}

/*//=======================================================
// ● 対数正規分布生成
//=======================================================*/
double Mt::lognorm_rand(double u, double sigma){
	std::lognormal_distribution<> dist(u, sigma);
	const double result = dist(mt_r);
	return result;
}


/*//=======================================================
// ● 指数分布
//=======================================================*/
double Mt::exp_rand(double lam){
	std::exponential_distribution<> dist(lam);
	const double result = dist(mt_r);
	return result;
}

/*//=======================================================
// ● カイ2乗分布（n：自由度）
//=======================================================*/
double Mt::chai_rand(double n){
	std::chi_squared_distribution<> dist(n);
	const double result = dist(mt_r);
	return result;
}


/*//=======================================================
// ● フィッシャーのF分布（m,n：自由度）
//=======================================================*/
double Mt::fisher_rand(double m, double n){
	std::fisher_f_distribution<> dist(m, n);
	const double result = dist(mt_r);
	return result;
}

/*//=======================================================
// ● studentのｔ分布生成（n：自由度）
//=======================================================*/
double Mt::student_t_rand(double n){
	std::student_t_distribution<> dist(n);
	const double result = dist(mt_r);
	return result;
}


/*//=======================================================
// ● ポアソン分布生成 
//=======================================================*/
double Mt::poisson_rand(double u){
	std::poisson_distribution<> dist(u);
	const double result = dist(mt_r);
	return result;
}


/*//=======================================================
// ● ガンマ分布生成
//=======================================================*/
double Mt::gamma_rand(double a, double b){
	std::gamma_distribution<> dist(a, b);
	const double result = dist(mt_r);
	return result;
}

/*//=======================================================
// ● ワイブル分布生成
//=======================================================*/
double Mt::weibull_rand(double a, double b){
	std::weibull_distribution<> dist(a, b);
	const double result = dist(mt_r);
	return result;
}

/*//=======================================================
// ● ベルヌーイ分布生成(確率aでtrueを生成) 
//=======================================================*/
bool Mt::bernoulli_rand(double a){
	std::bernoulli_distribution dist(a);
	const bool result = dist(mt_r);
	return result;
}


/*//=======================================================
// ● 二項分布生成(成功確率a5の事象をn回施行し、成功した回数を返す)
//=======================================================*/
int Mt::binomial_rand(double a, int n){
	std::binomial_distribution<> dist(n, a);
	const int result = dist(mt_r);
	return result;
}

/* end of namespace */
}
