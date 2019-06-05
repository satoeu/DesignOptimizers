#ifndef DEF_REXSTAR00_CROSS_DEF
#define DEF_REXSTAR00_CROSS_DEF

#include <BaseFuncs/Mt.hpp>
#include "../RGeneS.hpp"
#include "RootCGeneratorS.hpp"


/*
//=============================================================
// ■ REXstar
//=============================================================
//	REXstar交叉による子個体生成（ただし範囲拡張などはない簡単化ver）
//============================================================= */
class REXstar : public RootCGeneratorS{
protected:
	double rex_t;												/* 降下方向ステップサイズ */
	int parent_num;                                             /* 親個体の数 */
/*--------------------------------------------------------------------*/
	void selectNewParents(int* parents, int total_size);																				/* 親個体のランダム選択 */
	void crossREX(vector<RGeneS>& children, const vector<RGeneS>& populations, int* parents);											/* REX交叉 */
	void child_generationREX(RGeneS& child, const vector<RGeneS>& populations, int* parents, double *par_grav, double *best_genes);		/* 交叉実態部分 */
/*--------------------------------------------------------------------*/
public:
	vector<int> selected_parent;
	/* コンストラクタ */
	REXstar(int c_num, int p_num, double tt) : RootCGeneratorS(c_num){        
		child_num=c_num;
		parent_num = p_num;
		rex_t = tt;
	};
	~REXstar(){;};																	/* デストラクタ */
																						/* 子個体生成 */
	void create_children(const vector<RGeneS>& populations, vector<RGeneS>& children);
};


#endif

