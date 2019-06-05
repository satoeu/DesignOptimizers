#include "FunctionalOracle.hpp"

#include <TestObjFuncs/TestObjFuncs.hpp>


/*
//=============================================================
// ■ FunctionalOracle
//=============================================================
//	個体評価実施クラス（関数値オラクル）
//============================================================= */


/*//===========================================================
// ● コンストラクタ 
//=========================================================== */
FunctionalOracle::FunctionalOracle(){
	inputs = new double[OptSettings::GENE_SIZE];
}
/*//===========================================================
// ● デストラクタ 
//=========================================================== */
FunctionalOracle::~FunctionalOracle(){
    delete[] inputs;
}   

/*//===========================================================
// ● 遺伝子情報をこのオラクルにセット
//=========================================================== */
void FunctionalOracle::setGeneParas(double* gene_values){
	for(int i = 0 ; i < OptSettings::GENE_SIZE ; i++){
		inputs[i] = gene_values[i];
	}
}

/*//===========================================================
// ● 評価後にオラクル内の情報リセット用メソッド
//=========================================================== */
void FunctionalOracle::reset_oracle_data(){
    
}

/*//===========================================================
// ●  評価値計算
//=========================================================== */
int FunctionalOracle::run_oracle(double* result_values){
	//TestObjFuncs::MultiObj(2, result_values, 0, 2, inputs);
	result_values[0] = TestObjFuncs::SingleObj("Rastrigin", OptSettings::GENE_SIZE, inputs);
	
	/* 残差 */
	result_values[1] = 0;
	/* 特徴量 */
	for(int i = 0 ; i < OptSettings::GENE_SIZE ; i++){
		result_values[i+2] = inputs[i];
	}
	return 1;
}