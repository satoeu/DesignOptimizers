
#include "RGeneS.hpp"
#include <fstream>
#include <iomanip> //std::setw(int w)

int RGeneS::size;
int RGeneS::feature_num;


/*//===========================================================
  // ● コンストラクタ
  //=========================================================== */
RGeneS::RGeneS(){
	/* 各種初期化 */
	gene = new double[size];	
	adapts = 1.0e+12;
	simple_adapts = 1.0e+12;
	residuals = 1.0e+12;
	FeatureValues = new double[feature_num];
	evaluate_correct_type = -1;
}

RGeneS::~RGeneS(){
	delete[] gene;
	delete[] FeatureValues;
}



/*//===========================================================
  // ● コンストラクタ（代入）
  //=========================================================== */
RGeneS::RGeneS(const RGeneS& RG){
	gene = new double[size];	
	FeatureValues = new double[feature_num];
	*this = RG;
}

/*//===========================================================
  // ● コンストラクタ（move）
  //=========================================================== */
RGeneS::RGeneS(RGeneS&& RG){
	*this = std::move(RG);
}

/*//===========================================================
  // ● 代入オペレータ（遺伝子と適応度を全てコピー）
  //=========================================================== */
RGeneS& RGeneS::operator=(const RGeneS& RG){
	/* 遺伝子をコピー */
	for(int i = 0 ; i < size; i++){
		this->gene[i] = RG.gene[i];
	}
	/* 適応度を渡す */
	this->adapts = RG.adapts;
	this->simple_adapts = RG.simple_adapts;
	this->residuals = RG.residuals;
	/* 特性値を渡す */
	for(int i = 0 ; i < feature_num; i++){
		this->FeatureValues[i] = RG.FeatureValues[i];
	}
	this->evaluate_correct_type = RG.evaluate_correct_type;
	return *this;
}
/*//===========================================================
  // ● 代入オペレータ（遺伝子と適応度を全てコピー）
  //=========================================================== */
RGeneS& RGeneS::operator=(RGeneS&& RG) noexcept{
	this->gene = RG.gene;
	RG.gene = nullptr;
	this->adapts = RG.adapts;
	this->simple_adapts = RG.simple_adapts;
	this->FeatureValues = RG.FeatureValues;
	RG.FeatureValues = nullptr;
	this->residuals = RG.residuals;
	this->evaluate_correct_type= RG.evaluate_correct_type;
	return *this;
}



/*//===========================================================
  // ● 遺伝子座ｘの内容を遺伝子座にセット
  //=========================================================== */
void RGeneS::setGene(double *x){
	for(int i=0 ; i< size ; i++){
		gene[i] = x[i];
	}
	this->adjust();
}
/*//===========================================================
  // ● 遺伝子座の内容を x に渡す
  //=========================================================== */
void RGeneS::getGene(double *xx) const{
	for(int i=0 ; i< size ; i++){
		xx[i] = gene[i];
	}
}

/*//===========================================================
// ● 適合度セット(conv_type 0:発散、1：収束、-10：未評価、全然ダメ)
//=========================================================== */
void RGeneS::setAdapts(double* oracle_results, int conv_type){
	/* 制約違反残差初期化 */
	this->residuals = 0.0;

	/* ちゃんと収束しているとき＝セット */
	if(conv_type == 1){
		/**/
		/* 正評価フラグをonに */
		evaluate_correct_type = 1;

		simple_adapts = oracle_results[0];
		residuals = oracle_results[1];
		for(int i=0 ; i < feature_num ; i++){
			FeatureValues[i] = oracle_results[i+2];
		}		

		/* 残差を込みした適応度の計算 */
		this->set_ResidualPenerlyAdapt();

	/**/
	/* 未収束のとき・・・ */
	}else{
		/* ここはダメ評価フラグをonに */
		evaluate_correct_type = 2;

		//cout << "Adapt error .... "<< object_num << ", " << feature_num << endl;
		adapts = 1.0e+6;
		residuals = 1.0e+6;
		//set_oracle_penarly();
		for(int i = 0 ; i < feature_num; i ++){
			FeatureValues[i] = 1000;
		}
	}
}


/*//===========================================================
// ● 残差を込みした適応度の計算
//=========================================================== */
void RGeneS::set_ResidualPenerlyAdapt(){
	if(OptSettings::USE_ORACLE){
		set_oracle_penarly();
	}else{
		adapts = simple_adapts + residuals;
	}
}

/*//===========================================================
// ● オラクル制約処理で適応度補正
//=========================================================== */
void RGeneS::set_oracle_penarly(){
	/* 全ての目的関数へオラクル補正 */
	double temp_adapt;
	double objF = simple_adapts;
	double ORACLE_PARAMETER = OptSettings::ORACLE_PARAMETERS;
	/* オラクル処理 */
	double absDiff = fabs(objF - ORACLE_PARAMETER);
	if(objF <= ORACLE_PARAMETER && residuals < 1.0e-3){
		temp_adapt = -1.0 * absDiff;
	}else{
		double thres1 = absDiff / 3.0;
		double alpha=0;
		if(objF > ORACLE_PARAMETER && residuals < thres1){
			double xx = absDiff*(6.0*sqrt(3.0)-2.0)/6.0/sqrt(3.0) - residuals;
			alpha = xx / (absDiff-residuals);
		}else if(objF > ORACLE_PARAMETER && thres1 <= residuals && residuals <= absDiff){
			alpha = 1.0 - 1.0 / (2.0*sqrt(absDiff/residuals));
		}else if(objF > ORACLE_PARAMETER && residuals >= absDiff){
			alpha = sqrt(absDiff/residuals) / 2.0;
		}
		temp_adapt = alpha*absDiff + (1.0-alpha)*residuals;
	}
	//(min->F0)
	adapts = temp_adapt;
}

/*//===========================================================
  // ● 遺伝子情報をファイルへ書き出し
  //=========================================================== */
void RGeneS::PrintGene(const string& filename){
	fstream fpx(filename, std::ios::out);
	fpx << "adapts " << adapts << endl;	
	fpx << "simple_adapts " << simple_adapts << endl;	
	for(int i = 0 ; i < feature_num-1 ; i++){
		fpx << OptSettings::FEATURES_NAME[i] << ", " ;
	}
	fpx << OptSettings::FEATURES_NAME[feature_num-1] <<  endl;
	fpx << "residual " << residuals << endl;
	fpx << "Gene" << endl;
	for(int i = 0 ; i < size ; i++){
		fpx << scientific << gene[i] << endl;
	}
	fpx.close();
}
/*//===========================================================
// ● 遺伝子情報をファイルへ書き出し(JSON形式)
//=========================================================== */
void RGeneS::PrintGeneJSON(const string& filename){
	fstream fpx(filename, std::ios::out);
	fpx << "{"<<endl;
	/**/
	fpx << " \"adapts\" : " << adapts << ","<<endl;;
	fpx << " \"simple_adapts\" : " << simple_adapts << ","<<endl;
	/**/
	for(int i = 0 ; i < feature_num-1 ; i++){
		fpx << " \"" << OptSettings::FEATURES_NAME[i] << "\" : " << setw(15) << scientific << FeatureValues[i] << ", "<< endl;
	}
	fpx << " \"" << OptSettings::FEATURES_NAME[feature_num-1] << "\" : " << setw(15) << scientific << FeatureValues[feature_num-1] << ", "<< endl;
	/**/
	fpx << " \"residuals\" : " << residuals << "," << endl;
	/**/
	fpx << " \"Gene\" : {";
	for(int i = 0 ; i < size-1 ; i++){
		fpx << scientific << gene[i] << ", ";
	}
	fpx << scientific << gene[size-1] << "}" << endl;
	/**/
	fpx << "}"<<endl;

	fpx.close();
}
