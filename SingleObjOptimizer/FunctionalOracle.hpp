#ifndef DEF_ORACLE_FUNCTION
#define DEF_ORACLE_FUNCTION

#include <BaseDefines.hpp>
#include "01_OptSetting.hpp"
#include "FunctionalOracle.hpp"

/*
//=============================================================
// ■ FunctionalOracle
//=============================================================
//	個体評価実施クラス（関数値オラクル）
//============================================================= */
class FunctionalOracle{
private:
    double* inputs;
public:
    FunctionalOracle();                                 /* コンストラクタ */
    ~FunctionalOracle();                                /* デストラクタ */
    void setGeneParas(double* gene_values);             /* 遺伝子情報をこのオラクルにセット */
    void reset_oracle_data();                           /* 評価後にオラクル内の情報リセット用メソッド */
    int run_oracle(double* result_values);				/* 評価値計算 */
};


#endif
