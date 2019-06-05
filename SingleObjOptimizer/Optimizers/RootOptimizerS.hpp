#ifndef DEF_ROOT_OPTIMIZER_MODE
#define DEF_ROOT_OPTIMIZER_MODE

#include <BaseDefines.hpp>
#include "../01_OptSetting.hpp"
#include "../RPopulationS.hpp"


/*
//=============================================================
// ■ RootOptimizer
//=============================================================
//	最適化における解探索ソルバのスーパークラス
//============================================================= */
class RootOptimizerS{
protected:
    int population_size;	                                /* 個体集団の数 */									
public:
    RootOptimizerS(){return;};
    RootOptimizerS(int p_num){
        population_size = p_num;
    };                                                      /* コンストラクタ */
    ~RootOptimizerS(){;};                                    /* デストラクタ */
    /* 生存選択 */
    virtual void selection(int select_num, vector<RGeneS>& ParentPop, vector<RGeneS>& Children){
            return;
    };    
};

#endif
