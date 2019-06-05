#ifndef DEF_ROOT_CGENERATOR_MODE
#define DEF_ROOT_CGENERATOR_MODE

#include <BaseDefines.hpp>
#include "../01_OptSetting.hpp"
#include "../RPopulationS.hpp"


/*
//=============================================================
// ■ RootCGenerator
//=============================================================
//	最適化における解探索ソルバのスーパークラス
//============================================================= */
class RootCGeneratorS{
protected:
    int child_num;                                          /* 子個体生成数 */
public:
    RootCGeneratorS(){return;};
    RootCGeneratorS(int c_num){child_num=c_num;};             /* コンストラクタ */
    ~RootCGeneratorS(){;};                                    /* デストラクタ */
    /* 子個体生成 */
    virtual void create_children(const vector<RGeneS>& populations, const vector<RGeneS>& elite_data, vector<RGeneS>& children){
            return;
    };    
};

#endif
