#ifndef DEF_RGENE
#define DEF_RGENE

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <BaseDefines.hpp>
#include <vector>
#include <string>
#include "01_OptSetting.hpp"


/*
//=============================================================
// ■ 遺伝子座のクラス　for　多目的最適化
//=============================================================
//	GAに使用する個体の設定をまとめたクラス(多目的用, Fは最小化でセット)
//============================================================= */
class RGeneS{
protected:
	static int size;									/* 遺伝子座のサイズ */
	static int feature_num;								/* 適応度計算に使う特性値数 */
/*------------------------------------------------------------------*/
	double *gene;										/* 遺伝子座 */
	int evaluate_correct_type;							/* 正しく評価されたか（-1：未評価。+1：正常評価。+2：解析せずにダメと評価 */

	double adapts;										/* 適応度 */
	double simple_adapts;								/* 適応度(残差を考慮していない純適応度) */
	double residuals;									/* 制約非達成時の残差関数値 */
	double *FeatureValues;								/* 適応度計算に使った特性値の数 */
/*------------------------------------------------------------------*/
	void setResidual(double xx){residuals=xx;};			/* 制約残差関数値セッタ */
	void setNormalAdapts(double x){adapts = x;};
	void set_ResidualPenerlyAdapt();					/* 残差を込みした適応度の計算 */
	void set_oracle_penarly();							/* オラクル制約処理で適応度補正 */
/*------------------------------------------------------------------*/
	void adjust();													/* 遺伝子補正 */
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
public:
	RGeneS();															/* コンストラクタ */
	~RGeneS();															/* デストラクタ */
	RGeneS(const RGeneS& RG);	
	RGeneS(RGeneS&& RG);	
	RGeneS& operator=(const RGeneS& RG);									/* 代入オペレータ（遺伝子と適応度を全てコピー） */
	RGeneS& operator=(RGeneS&& RG) noexcept;								/* 代入オペレータ（遺伝子と適応度を全てコピー） */
/*------------------------------------------------------------------*/
	static void setSize(int s){size=s;};								/* 遺伝子サイズ */
	static void setFeatureNum(int s) { feature_num = s; };				/* 特性値数 */
	static int getSize(){return size;};
	static int gettFeatureNum(){return feature_num;};
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
	void setGene(double *x);											/* 遺伝子座ｘの内容を遺伝子座にセット */
	void getGene(double *xx) const;										/* 遺伝子座の内容を x に渡す */
	void setGene(int pos, double x){gene[pos]=x;};						/* 遺伝子座ｘの内容を遺伝子座にセット */
	double getGene(int i) const{return gene[i];};						/* 遺伝子座 x[i] の内容を返す */
	int getEvaluateCorrect() const{return evaluate_correct_type;};

	double getAdapts() const{return adapts;};							/* 適応度ゲット */
	double getSimpleAdapts() const{return simple_adapts;};				/* 修正適応度ゲット */
	void setAdaptNormal(double xx){adapts=xx;};							/* 適応度単純セット */
	void setSimpleAdapt(double xx){simple_adapts=xx;};					/* 修正適応度単純セット */
	void setFeatureVal(int p, double xx){FeatureValues[p]=xx;};			/* 特性値セッタ */
	double getFeatureVal(int p) const{return FeatureValues[p];};		/* 特性値ゲッタ */
	double getResidual() const{return residuals;};						/* 制約残差関数値ゲッタ */
/*------------------------------------------------------------------*/
	void makeGene();													/* 遺伝子座初期化 */
	bool check_simple_violate_init();									/* 単純にチェックできる制約を評価（初期個体生成時） */
	bool check_simple_violate_child();									/* 単純にチェックできる制約を評価（子個体生成時） */
	void setAdapts(double* oracle_results, int conv_type);				/* 適応度セット */
/*------------------------------------------------------------------*/
	void PrintGene(const string& filename);								/* 遺伝子情報をファイルへ書き出し */
	void PrintGeneJSON(const string& filename);							/* 遺伝子情報をファイルへ書き出し(JSON形式) */
/*------------------------------------------------------------------*/
};

#endif

