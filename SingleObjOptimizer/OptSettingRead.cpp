#include "01_OptSetting.hpp"
#include <fstream>

/*
//=======================================================
// ■ OptSettings
//=======================================================
// 最適化の定数系管理staticクラス
//=======================================================*/

int OptSettings::NUM_FEATURE;
int OptSettings::GENE_SIZE;
int OptSettings::PARA_OMP_CALC_NUM;
int OptSettings::TOTAL_SIZE;
int OptSettings::GENERATION_SIZE;
bool OptSettings::USE_ORACLE;
double OptSettings::ORACLE_PARAMETERS;
vector<string> OptSettings::FEATURES_NAME;
/*//===========================================================//*/
int OptSettings::REX_Settings::NUM_CHILD;
int OptSettings::REX_Settings::PARENT_NUM;



/*//===========================================================
// ● 最適化設定読み取り
//=========================================================== */
void OptSettings::opt_set_read(){
	fstream fp("001_BaseOptSettings.dat", std::ios::in);
	string temp_str;
    /* 基本read */
	fp >> temp_str >> GENE_SIZE;
	fp >> temp_str >> NUM_FEATURE;
	fp >> temp_str >> PARA_OMP_CALC_NUM;
	fp >> temp_str >> TOTAL_SIZE;
	fp >> temp_str >> GENERATION_SIZE;

    fp >> temp_str;

    /* Oracle読み取り */
    int temp_value;
	fp >> temp_str >> temp_value >> ORACLE_PARAMETERS;
    USE_ORACLE = (temp_value==1);
	/**/
	/**/
	/* 特性値のお名前取得 */
	fp >> temp_str;
	fp >> temp_str;
	fp >> temp_str;
	for(int i=0 ; i < NUM_FEATURE; i++){
		string str00;
		fp >> str00;
		FEATURES_NAME.push_back(str00);
	}
	/**/
	fp.close();


	/**/
	/* 読み終わったら、RGeneのstaticをセットする */
	OptSettings::setRGeneStaticParas();
}


/*//===========================================================
// ● REXstarの設定読み取り
//=========================================================== */
void OptSettings::REX_Settings::opt_set_read_rex(){
	fstream fp("021_REXgeneratorSettings.dat", std::ios::in);
	string temp_str;
	/* 基本read */
	fp >> temp_str >> NUM_CHILD;
	fp >> temp_str >> PARENT_NUM;

	fp.close();
}




/*//===========================================================
// ● RGeneのStatic変数をセットする
//=========================================================== */
void OptSettings::setRGeneStaticParas(){
	/* Gene定数static初期化 */
	RGeneS::setSize(OptSettings::GENE_SIZE);
	RGeneS::setFeatureNum(OptSettings::NUM_FEATURE);
}
