#ifndef DEF_REXs_CROSS_DEF
#define DEF_REXs_CROSS_DEF

#include <BaseFuncs/Mt.hpp>
#include "../RGeneS.hpp"
#include "RootCGeneratorS.hpp"


/*
//=============================================================
// ■ AREX
//=============================================================
//	AREX交叉による子個体生成（ただし範囲拡張などはない簡単化ver）
//============================================================= */
class SimpleAREX : public RootCGeneratorS{
protected:
    int parent_num;                                             /* 親個体の数 */
	/*--------------------------------------------------------------------*/
	void selectNewParents(int* parents, int total_size);																				/* 親個体のランダム選択 */
    void crossREX(vector<RGeneS>& children, const vector<RGeneS>& populations, int* parents);											/* REX交叉 */
    void child_generationREX(RGeneS& child, const vector<RGeneS>& populations, int* parents, double *par_grav, double *best_genes);		/* 交叉実態部分 */
	/*--------------------------------------------------------------------*/
public:
	vector<int> selected_parent;
    /* コンストラクタ */
	SimpleAREX(int c_num, int p_num) : RootCGeneratorS(c_num){        
        child_num=c_num;
        parent_num = p_num;
    };
	~SimpleAREX(){;};																	/* デストラクタ */
	/* 子個体生成 */
    void create_children(const vector<RGeneS>& populations, vector<RGeneS>& children);
};





#endif
